package org.processmining.antialignments.algorithm;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

import java.util.Arrays;

import lpsolve.LpSolve;
import lpsolve.LpSolveException;

public class LPMatrix {

	public static final int INFEASIBLE = 0;
	public static final int OPTIMAL = 0;
	public static int GE = LpSolve.GE;
	public static int LE = LpSolve.LE;
	public static int EQ = LpSolve.EQ;
	public static double INFINITE = 1.0e10;

	private final double[][] matrix;
	private final double[] rhs;
	private final int[] types;
	private final double[] obj;
	private final double[] upBo;
	private final double[] lowBo;
	private final boolean[] integer;
	private boolean isMinimizing = true;
	private final String[] colNames;
	private final String[] rowNames;

	public LPMatrix(int rows, int columns) {
		matrix = new double[rows][columns];
		rhs = new double[rows];
		types = new int[rows];
		upBo = new double[columns];
		Arrays.fill(upBo, INFINITE);
		lowBo = new double[columns];
		obj = new double[columns];
		integer = new boolean[columns];
		colNames = new String[columns];
		rowNames = new String[rows];
	}

	public double getMat(int row, int column) {
		return matrix[row][column];
	}

	public void setMat(int row, int column, double value) {
		matrix[row][column] = value;
	}

	public void setInt(int column, boolean isInt) {
		integer[column] = isInt;
	}

	public void setUpbo(int column, double value) {
		upBo[column] = value;
	}

	public void setLowbo(int column, double value) {
		lowBo[column] = value;
	}

	public void setColName(int column, String name) {
		colNames[column] = name;
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
		System.arraycopy(rhs, 0, this.rhs, 0, this.rhs.length);
	}

	public int getNcolumns() {
		return integer.length;
	}

	public void setRh(int row, double value) {
		rhs[row] = value;
	}

	public void setObjective(int column, double value) {
		obj[column] = value;
	}

	public double getRh(int row) {
		return rhs[row];
	}

	public LpSolve toLpSolve() throws LpSolveException {
		LpSolve lp = LpSolve.makeLp(0, getNcolumns());
		lp.setAddRowmode(true);
		double[] row = new double[getNcolumns() + 1];
		System.arraycopy(obj, 0, row, 1, obj.length);
		lp.setObjFn(row);

		for (int r = 0; r < matrix.length; r++) {
			System.arraycopy(matrix[r], 0, row, 1, matrix[r].length);
			lp.addConstraint(row, types[r], rhs[r]);
		}
		lp.setAddRowmode(false);

		for (int r = 0; r < matrix.length; r++) {
			if (rowNames[r] != null) {
				lp.setRowName(r + 1, rowNames[r]);
			}
		}
		for (int c = 0; c < getNcolumns(); c++) {
			lp.setInt(c + 1, integer[c]);
			lp.setLowbo(c + 1, lowBo[c]);
			lp.setUpbo(c + 1, upBo[c]);
			if (colNames[c] != null) {
				lp.setColName(c + 1, colNames[c]);
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

		return lp;

	}

	public GRBModel toGurobi(GRBEnv env) throws GRBException {
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
	}

	public String getColName(int column) {
		return colNames[column];
	}

	public void printLp() {
		LpSolve lp = null;
		try {
			lp = toLpSolve();
			lp.printLp();
		} catch (LpSolveException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			lp.deleteAndRemoveLp();
		}

	}

}
