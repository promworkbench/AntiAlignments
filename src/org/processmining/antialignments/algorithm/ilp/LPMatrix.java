package org.processmining.antialignments.algorithm.ilp;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

import java.util.Arrays;

import lpsolve.LpSolve;
import lpsolve.LpSolveException;

public abstract class LPMatrix<S> {

	public static final int INFEASIBLE = 1;
	public static final int OPTIMAL = 0;
	public static final int OTHER = 1;

	public static int GE = LpSolve.GE;
	public static int LE = LpSolve.LE;
	public static int EQ = LpSolve.EQ;
	public static double INFINITE = 1.0e30;

	public static class LPMatrixException extends Exception {

		/**
		 * 
		 */
		private static final long serialVersionUID = 798019891801690375L;

		public LPMatrixException(Throwable cause) {
			super(cause);
		}
	}

	public static class LPSOLVE extends LPMatrix<LpSolve> {

		public LPSOLVE(int rows, int columns) {
			super(rows, columns, 1);
		}

		public LpSolve toSolver() throws LPMatrixException {
			LpSolve lp;
			try {
				lp = LpSolve.makeLp(0, getNcolumns());
				lp.setAddRowmode(true);
				lp.setObjFn(obj);

				for (int r = 0; r < matrix.length; r++) {
					lp.addConstraint(matrix[r], types[r], rhs[r]);
				}
				lp.setAddRowmode(false);

				for (int r = 0; r < matrix.length; r++) {
					if (rowNames[r] != null) {
						lp.setRowName(r + 1, rowNames[r]);
					}
				}
				for (int c = 1; c <= getNcolumns(); c++) {
					lp.setInt(c, integer[c]);
					lp.setLowbo(c, lowBo[c]);
					lp.setUpbo(c, upBo[c]);
					if (colNames[c] != null) {
						lp.setColName(c, colNames[c]);
					}
				}

				if (isMinimizing) {
					lp.setMinim();
				} else {
					lp.setMaxim();
				}

				lp.setVerbose(1);

				lp.setScaling(LpSolve.SCALE_GEOMETRIC | LpSolve.SCALE_EQUILIBRATE | LpSolve.SCALE_INTEGERS);
				lp.setScalelimit(5);
				lp.setPivoting(LpSolve.PRICER_DEVEX | LpSolve.PRICE_ADAPTIVE);
				lp.setMaxpivot(250);
				lp.setBbFloorfirst(LpSolve.BRANCH_AUTOMATIC);
				lp.setBbRule(LpSolve.NODE_PSEUDONONINTSELECT | LpSolve.NODE_GREEDYMODE | LpSolve.NODE_DYNAMICMODE
						| LpSolve.NODE_RCOSTFIXING);
				lp.setBbDepthlimit(-50);
				lp.setAntiDegen(LpSolve.ANTIDEGEN_FIXEDVARS | LpSolve.ANTIDEGEN_STALLING);
				lp.setImprove(LpSolve.IMPROVE_DUALFEAS | LpSolve.IMPROVE_THETAGAP);
				lp.setBasiscrash(LpSolve.CRASH_NOTHING);
				lp.setSimplextype(LpSolve.SIMPLEX_DUAL_PRIMAL);

				//				lp.printLp();

				return lp;
			} catch (LpSolveException e) {
				throw new LPMatrixException(e);
			}

		}

		public int solve(double[] vars) throws LPMatrixException {

			LpSolve lp = toSolver();

			try {
				solverResult = lp.solve();

				if (solverResult == lp.INFEASIBLE) {
					return INFEASIBLE;
				} else if (solverResult == OPTIMAL) {
					//					double[] varTmp = new double[vars.length + 1];
					lp.getVariables(vars);
					//					System.arraycopy(varTmp, offset, vars, 0, vars.length);
					return OPTIMAL;
				} else {
					return OTHER;
				}

			} catch (LpSolveException e) {
				throw new LPMatrixException(e);

			} finally {

				lp.deleteAndRemoveLp();
			}
		}

	}

	public static class GUROBI extends LPMatrix<GRBModel> {

		private final GRBEnv env;

		public GUROBI(GRBEnv env, int rows, int columns) {
			super(rows, columns, 0);
			this.env = env;
		}

		public GRBModel toSolver() throws LPMatrixException {
			try {
				GRBModel model = new GRBModel(env);

				// Add variables to the model
				char[] vType = new char[integer.length];
				for (int v = 0; v < vType.length; v++) {
					vType[v] = integer[v] ? GRB.INTEGER : GRB.CONTINUOUS;
				}

				double[] minObj = new double[obj.length];
				for (int i = 0; i < obj.length; i++) {
					if (isMinimizing) {
						minObj[i] = obj[i];
					} else {
						minObj[i] = -obj[i];
					}
				}
				GRBVar[] vars = model.addVars(lowBo, upBo, minObj, vType, colNames);
				model.update();

				// Populate A matrix
				char sense;
				for (int i = 0; i < getNrows(); i++) {
					GRBLinExpr expr = new GRBLinExpr();
					for (int j = 0; j < getNcolumns(); j++)
						if (matrix[i][j] != 0)
							expr.addTerm(matrix[i][j], vars[j]);
					sense = types[i] == EQ ? GRB.EQUAL : types[i] == LE ? GRB.LESS_EQUAL : GRB.GREATER_EQUAL;
					model.addConstr(expr, sense, rhs[i], rowNames[i]);
				}

				return model;
			} catch (GRBException e) {
				throw new LPMatrixException(e);
			}
		}

		public int solve(double[] vars) throws LPMatrixException {
			GRBModel grbModel = toSolver();

			try {
				grbModel.optimize();
				solverResult = grbModel.get(GRB.IntAttr.Status);

				if (grbModel.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {

					for (int j = 0; j < vars.length; j++)
						vars[j] = grbModel.getVars()[j].get(GRB.DoubleAttr.X);
					return OPTIMAL;
				} else if (grbModel.get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE
						|| grbModel.get(GRB.IntAttr.Status) == GRB.Status.INF_OR_UNBD) {
					return INFEASIBLE;

				} else {
					return OTHER;
				}

			} catch (GRBException e) {
				throw new LPMatrixException(e);

			} finally {
				grbModel.dispose();
			}
		}

	}

	protected final int offset;
	protected final double[][] matrix;
	protected final double[] rhs;
	protected final int[] types;
	protected final double[] obj;
	protected final double[] upBo;
	protected final double[] lowBo;
	protected final boolean[] integer;
	protected boolean isMinimizing = true;
	protected final String[] colNames;
	protected final String[] rowNames;
	protected int solverResult;

	protected LPMatrix(int rows, int columns, int offset) {
		this.offset = offset;
		matrix = new double[rows][columns + offset];
		rhs = new double[rows];
		types = new int[rows];
		Arrays.fill(types, LPMatrix.EQ);
		upBo = new double[columns + offset];
		Arrays.fill(upBo, INFINITE);
		lowBo = new double[columns + offset];
		obj = new double[columns + offset];
		integer = new boolean[columns + offset];
		colNames = new String[columns + offset];
		rowNames = new String[rows];
	}

	public abstract S toSolver() throws LPMatrixException;

	public abstract int solve(double[] vars) throws LPMatrixException;

	public int getLastSolverResult() {

		return solverResult;

	}

	public double getMat(int row, int column) {
		return matrix[row][column + offset];
	}

	public void setMat(int row, int column, double value) {
		matrix[row][column + offset] = value;
	}

	public void setInt(int column, boolean isInt) {
		integer[column + offset] = isInt;
	}

	public void setUpbo(int column, double value) {
		upBo[column + offset] = value;
	}

	public void setLowbo(int column, double value) {
		lowBo[column + offset] = value;
	}

	public void setColName(int column, String name) {
		colNames[column + offset] = name;
	}

	public void setRowName(int row, String name) {
		rowNames[row] = name;
	}

	public void setConstrType(int row, int type) {
		types[row] = type;
	}

	public void setMaxim() {
		isMinimizing = false;
	}

	public void setMinim() {
		isMinimizing = true;
	}

	public int getNrows() {
		return matrix.length;
	}

	public void setRhVec(double[] rhs) {
		System.arraycopy(rhs, 0, this.rhs, 0, rhs.length);
	}

	public int getNcolumns() {
		return integer.length - offset;
	}

	public void setRh(int row, double value) {
		rhs[row] = value;
	}

	public void setObjective(int column, double value) {
		obj[column + offset] = value;
	}

	public double getRh(int row) {
		return rhs[row];
	}

	public String getColName(int column) {
		return colNames[column + offset];
	}

	public double getObjective(int column) {
		return obj[column + offset];
	}

	public String[] getColNames() {
		return colNames;
	}

	//	public void printLp() {
	//		LpSolve lp = null;
	//		try {
	//			lp = toLpSolve();
	//			lp.printLp();
	//		} catch (LpSolveException e) {
	//			// TODO Auto-generated catch block
	//			e.printStackTrace();
	//		} finally {
	//			lp.deleteAndRemoveLp();
	//		}
	//
	//	}

}
