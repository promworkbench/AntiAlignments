package org.processmining.antialignments.alignments;

import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;

import java.util.Arrays;
import java.util.Stack;

import lpsolve.LpSolve;

import org.processmining.antialignments.algorithm.AbstractILPCalculator;
import org.processmining.antialignments.algorithm.ilp.LPMatrix;
import org.processmining.antialignments.algorithm.ilp.LPMatrix.LPMatrixException;
import org.processmining.framework.util.Pair;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetEdge;
import org.processmining.models.graphbased.directed.petrinet.PetrinetGraph;
import org.processmining.models.graphbased.directed.petrinet.PetrinetNode;
import org.processmining.models.graphbased.directed.petrinet.elements.Arc;
import org.processmining.models.graphbased.directed.petrinet.elements.InhibitorArc;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class AlignmentILPCalculator extends AbstractILPCalculator {

	private static final double EPSILON = 0.1;
	private int spCols;
	private int spRows;
	private int synchronousTransitions;
	private short[] syncTransitionMap;
	private short[] syncLabelMap;
	private short[] logMoveMap;
	private short[] syncEventMap;

	public AlignmentILPCalculator(PetrinetGraph net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short, TShortObjectMap<String> short2label, short[][] log) {
		super(net, initialMarking, finalMarking, label2short, short2label, log);
		mode = MODE_LPSOLVE;
		//		mode = MODE_LPSOLVE;
		//		try {
		//			gbEnv = new GRBEnv();
		//			gbEnv.set(GRB.IntParam.OutputFlag, 0);
		//			mode = MODE_GUROBI;
		//		} catch (GRBException _) {
		//			mode = MODE_LPSOLVE;
		//		}
		//
	}

	public void getAlignments(Marking initialMarking, Marking finalMarking) throws LPMatrixException {

		mode = MODE_LPSOLVE;

		LPMatrix<?> matrix = setupLpForHybrid(3, true, initialMarking, finalMarking, 0, 0);
		((LpSolve) matrix.toSolver()).printLp();

	}

	protected LPMatrix<?> setupLpForHybrid(int maxLengthX, boolean integerVariables, Marking initialMarking,
			Marking finalMarking, int traceToConsider, int startTraceAt) {

		// Find the transitions in the net labeled with an event in 
		// log[traceToConsider][startTraceAt] .. log[traceToConsider][startTraceAt+lengthX-1]

		// A synchronous transition is needed for each (t,e) combination where e is an event

		short[] label2pos = new short[labels];
		{
			Arrays.fill(label2pos, (short) -1);

			logMoveMap = new short[Math.min(log[traceToConsider].length - startTraceAt, maxLengthX)];

			TShortList trs = new TShortArrayList();
			TShortList evts = new TShortArrayList();
			TShortList indices = new TShortArrayList();
			short pos = 0;
			for (short e = (short) startTraceAt; e < log[traceToConsider].length && e < startTraceAt + maxLengthX; e++) {
				for (short t = invisibleTransitions; t < transitions; t++) {
					if (equalLabel(t, log[traceToConsider][e])) {
						if (label2pos[log[traceToConsider][e]] < 0) {
							label2pos[log[traceToConsider][e]] = pos;
							pos++;
						}
						trs.add(t);
						evts.add(log[traceToConsider][e]);
						indices.add((short) (e - startTraceAt));
					}
				}
				logMoveMap[e - startTraceAt] = log[traceToConsider][e];
			}
			syncTransitionMap = trs.toArray();
			syncLabelMap = evts.toArray();
			syncEventMap = indices.toArray();

			for (short l = 0; l < labels; l++) {
				if (label2pos[l] < 0) {
					label2pos[l] = pos;
					pos++;
				}
			}
		}

		synchronousTransitions = syncLabelMap.length;
		spCols = transitions + synchronousTransitions + logMoveMap.length;
		spRows = places + logMoveMap.length + 1;

		//-  
		//-                            AB   >= m0
		//-                            AAB  >= m0
		//-                            AAAA == mf        
		LPMatrix<?> lp = setupMatrix(spRows * maxLengthX + (places + labels) + //
				// Xi - Xi+1<=0        Xi <= 1  Y <= ly*X_lx
				(maxLengthX - 1) + maxLengthX + 1,//
				//  X part is small           Y part is full
				spCols * maxLengthX + 2 * transitions - invisibleTransitions + labels);

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short p, t;
			int dir;
			int type = LpSolve.GE;
			Transition trans;
			if (e instanceof Arc) {
				if (e.getSource() instanceof Place) {
					p = place2int.get(e.getSource());
					trans = (Transition) e.getTarget();
					t = trans2int.get(trans);
					dir = -((Arc) e).getWeight();
				} else {
					trans = (Transition) e.getSource();
					t = trans2int.get(trans);
					p = place2int.get(e.getTarget());
					dir = ((Arc) e).getWeight();
				}
			} else if (e instanceof InhibitorArc) {
				p = place2int.get(e.getSource());
				trans = (Transition) e.getTarget();
				t = trans2int.get(trans);
				dir = 0;
				type = LpSolve.EQ;
			} else {
				continue;
			}

			short tLabel = trans2label[t];

			int labelRow;// label is the relative position of the corresponding label row.
			if (tLabel >= 0) {
				labelRow = places - p + label2pos[tLabel];
			} else {
				labelRow = 0;
			}

			for (int block = 0; block <= maxLengthX; block++) {
				// First the whole A matrices.
				int r = block * spRows + p;

				lp.setConstrType(r, type);
				lp.setRowName(r, "A" + block + "_" + p);
				for (int c = t; c < block * spCols; c += spCols) {
					// update the SP A matrix
					lp.setMat(r, c, lp.getMat(r, c) + dir);
					// update the SP A' matrix for mapped visible transitions

					// Synchronous transitions (SYNCMOVE)
					for (int sm = 0; sm < synchronousTransitions; sm++) {
						if (syncTransitionMap[sm] == t) {
							int pos = c - t + transitions + sm;
							lp.setMat(r, pos, lp.getMat(r, pos) + dir);
							if (block < maxLengthX) {
								// determine consumption row for event class
								int ecr = syncEventMap[sm];
								// token flow
								lp.setMat(r - p + places + ecr, pos, -1);
								lp.setMat(r - p + places + ecr + 1, pos, 1);
							} else {
								// count
								lp.setMat(r + labelRow, pos, 1);
							}
						}
					}
				}

				// Then, the  A- matrix.
				int c = block * spCols + t;
				if ((dir < 0 || trans.isInvisible()) && block < maxLengthX) {
					lp.setColName(c, "M_" + trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
					lp.setObjective(c, getCostForModelMove(trans));
					lp.setInt(c, integerVariables);
					if (!trans.isInvisible()) {
						// set upper bound of 1 for visible transitions
						lp.setUpbo(c, 1.0);
					}
					// Synchronous transitions (SYNCMOVE)
					for (int sm = 0; sm < synchronousTransitions; sm++) {
						if (syncTransitionMap[sm] == t) {
							int pos = c - t + transitions + sm;

							lp.setColName(pos, "S_" + trans.getLabel() + sm + "-" + block);
							lp.setObjective(pos, getCostForSync(trans, syncLabelMap[sm]));

							lp.setInt(pos, integerVariables);
							lp.setUpbo(pos, 1.0);
							lp.setMat(r, pos, lp.getMat(r, pos) + dir);

							// determine consumption row for event class
							int ecr = syncEventMap[sm];
							lp.setMat(r - p + places + ecr, pos, -1);
							// lp.setMat(r + label + 1, c + pos, 1);
						}
					} // update the A matrix only for consumption and for invisible transitions
						// or in the last block

					// update the SP A matrix for consumption in case of visible,
					// or always in case of invisible transition
					lp.setMat(r, c, lp.getMat(r, c) + dir);
				}

				// Then, in the last block
				if (block == maxLengthX) {

					// In this Block, the positions of the labels are now different,
					// as are the positions of the synchronous transitions.

					// The relative position of the synchronous counterpart is at:
					int pos = transitions - invisibleTransitions;

					if (finalMarking != null) {
						lp.setConstrType(r, LpSolve.EQ);
					}
					// update the SP A matrix
					lp.setMat(r, c, lp.getMat(r, c) + dir);
					lp.setColName(c, "M_" + trans.getLabel().replace("\\n$invisible$", ""));
					lp.setObjective(c, getCostForModelMove(trans) + EPSILON);

					lp.setInt(c, integerVariables);

					lp.setColName(c + pos, "S_" + trans.getLabel() + "-" + block);
					lp.setObjective(c + pos, getCostForSync(trans, tLabel) + EPSILON);

					lp.setInt(c + pos, integerVariables);
					lp.setMat(r, c + pos, lp.getMat(r, c + pos) + dir);
					if (labelRow > 0) {
						lp.setMat(r + labelRow, c + pos, 1);
					}

				}

			}
		}

		for (short l = 0; l < logMoveMap.length; l++) {
			for (int block = 0; block <= maxLengthX; block++) {
				// First the whole A matrices.
				int r = block * spRows + places + l;
				lp.setConstrType(r, LPMatrix.GE);
				lp.setRowName(r, "B" + block + "_" + l);
				if (l == logMoveMap.length - 1) {
					lp.setConstrType(r + 1, LPMatrix.GE);
					lp.setRowName(r + 1, "B" + block + "_" + (l + 1));
				}
				for (int c = transitions + synchronousTransitions + l; c < block * spCols; c += spCols) {
					// update the SP B matrix
					if (block < maxLengthX) {
						// token flow
						lp.setMat(r, c, -1);
						lp.setMat(r + 1, c, 1);
					} else {
						// count
						lp.setMat(r, c, 1);
					}
					lp.setUpbo(c, 1);
					lp.setInt(c, integerVariables);
				}

				// Then, the  A- matrix.
				int c = block * spCols + transitions + synchronousTransitions + l;
				if (block < maxLengthX) {
					lp.setColName(c, "L_" + short2label.get(logMoveMap[l]) + l + "-" + block);
					lp.setObjective(c, getCostForLogMove(logMoveMap[l]));
					lp.setInt(c, integerVariables);
					lp.setUpbo(c, 1.0);
					lp.setMat(r, c, -1);
				}

			}
		}
		for (short l = 0; l < label2pos.length; l++) {
			int r = maxLengthX * spRows + places + label2pos[l];
			lp.setRowName(r, "C_" + short2label.get(l));

			int c = maxLengthX * spCols + 2 * transitions - invisibleTransitions + label2pos[l];

			// Then, in the last block
			lp.setConstrType(r, LpSolve.EQ);
			// count the number of B's
			lp.setMat(r, c, 1);
			lp.setColName(c, "L_" + short2label.get(l));
			lp.setObjective(c, getCostForLogMove(l) + EPSILON);
			lp.setInt(c, integerVariables);

		}

		int row = spRows * maxLengthX + places + labels;
		// set up the comparisons and sums
		for (int i = 0; i < maxLengthX - 1; i++) {
			for (int t = i * spCols + invisibleTransitions; t < (i + 1) * spCols; t++) {
				lp.setMat(row + i, t, 1);
			}
			for (int t = (i + 1) * spCols + invisibleTransitions; t < (i + 2) * spCols; t++) {
				lp.setMat(row + i, t, -1);
			}
			lp.setRowName(row + i, "X" + i + "~X" + (i + 1));
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// Y < lm*l_lx
		row += maxLengthX - 1;
		for (int t = (maxLengthX - 1) * spCols + invisibleTransitions; t < maxLengthX * spCols; t++) {
			//TODO: Find some constant!
			lp.setMat(row, t, transitions * transitions);
		}
		for (int t = maxLengthX * spCols + invisibleTransitions; t < lp.getNcolumns(); t++) {
			lp.setMat(row, t, -1);
		}
		lp.setRowName(row, "Y~X" + maxLengthX);
		lp.setConstrType(row, LpSolve.GE);

		// sos constraints, sum max 1
		row++;
		for (int i = 0; i < maxLengthX; i++) {
			for (int t = i * spCols + invisibleTransitions; t < (i + 1) * spCols; t++) {
				lp.setMat(row + i, t, 1);
			}
			lp.setRowName(row + i, "X" + i + ".1");
			lp.setConstrType(row + i, LpSolve.LE);
		}

		lp.setMinim();

		double[] rhs = new double[lp.getNrows()];

		// first -initial Marking
		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row < spRows * maxLengthX + places; row += spRows) {
				rhs[row] = -initialMarking.occurrences(p);
			}
			if (log[traceToConsider].length > startTraceAt) {
				for (row = places; row < spRows * maxLengthX; row += spRows) {
					rhs[row] = -1;
				}
			}
		}
		// then add final marking
		if (finalMarking != null) {
			for (Place p : finalMarking.baseSet()) {
				row = maxLengthX * spRows + place2int.get(p);
				rhs[row] += finalMarking.occurrences(p);
			}
		}
		// then, the parikh vector of the log
		row = maxLengthX * spRows + places;
		// Count the number of times each label is in the trace from position startTraceAt
		for (int e = startTraceAt; e < log[traceToConsider].length; e++) {
			rhs[row + label2pos[log[traceToConsider][e]]]++;
		}

		// then X.i - X.i+1 >= 0
		// then Y   - ly * X.lx >= 0

		// then X.1 <= 1
		for (row = maxLengthX * spRows + places + labels + maxLengthX; row < rhs.length; row++) {
			rhs[row] = 1;
		}

		lp.setRhVec(rhs);

		return lp;
	}

	private double getCostForModelMove(Transition trans) {
		// TODO Auto-generated method stub
		return trans.isInvisible() ? 0 : 1;
	}

	private double getCostForLogMove(short s) {
		// TODO Auto-generated method stub
		return 1;
	}

	private double getCostForSync(Transition trans, short s) {
		// TODO Auto-generated method stub
		return 0;
	}

	protected void solveSequential(Marking initialMarking, Marking finalMarking, final int traceToConsider)
			throws LPMatrixException {

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		Marking marking = initialMarking;
		int startTracesAt = 0;
		LPMatrix<?> matrix;
		Stack<Pair<Transition, Short>> moves = new Stack<>();

		cutOffLength = 3;

		do {
						if (VERBOSE) {
							System.out.println("Trying to get from " + initialMarking + " to " + finalMarking + //
									" starting with  " + cutOffLength + " exact steps.");
						}

			matrix = setupLpForHybrid(cutOffLength, true, marking, finalMarking, traceToConsider, startTracesAt);
			//			((LpSolve) matrix.toSolver()).printLp();

			// Compute the new marking
			double[] vars = new double[matrix.getNcolumns()];
			int result = matrix.solve(vars);

			if (result == LPMatrix.OPTIMAL) {

				// Get the intermediate Marking reached after the cutOffLength
				marking = getIntermediateMarking(marking, matrix, vars);

				// Compute the number of explained events
				int l = updateListOfMoves(vars, log[traceToConsider], startTracesAt, moves);
				startTracesAt += l;

				if (VERBOSE) {
					System.out.println("Reached " + marking + " now.");
					System.out.println("Partial alignment :" + moves.toString());

				}
			} else {
				// BACKTRACKING NEEDED
				break;
			}
		} while (!marking.equals(finalMarking));

		if (VERBOSE) {
			System.out.println();
			System.out.println("Alignment done:");
			System.out.println(moves);
			System.out.println();

		}
		// translate vars into log, sync and model moves
		//		System.out.println(Arrays.toString(matrix.getColNames()));
		//		System.out.println(Arrays.toString(vars));

	}

	protected Marking getIntermediateMarking(Marking marking, LPMatrix<?> matrix, double[] vars) {
		Marking newMarking = new Marking(marking);
		for (short p = 0; p < places; p++) {
			// compute the effect of the x vectors on place p.
			double v = matrix.product(vars, 0, cutOffLength * spCols, cutOffLength * spRows + p);
			while (v < -.5) {
				newMarking.remove(short2place[p]);
				v += 1.0;
			}
			if (v > 0) {
				newMarking.add(short2place[p], (int) (v + .5));
			}
		}
		return newMarking;
	}

	private int updateListOfMoves(double[] vars, short[] trace, int startTraceAt, Stack<Pair<Transition, Short>> moves) {
		int l = 0;
		for (int c = 0; c < cutOffLength * spCols; c++) {
			if (vars[c] > 0.5 || c % spCols < invisibleTransitions) {
				if (c % spCols < transitions) {
					int v = (int) (vars[c] + 0.5);
					while (v > 0) {
						moves.push(new Pair<>(short2trans[c % spCols], Short.MIN_VALUE));
						v--;
					}
				} else if (c % spCols < transitions + synchronousTransitions) {
					assert (int) (vars[c] + 0.5) == 1;
					int t = syncTransitionMap[(c % spCols) - transitions];
					moves.push(new Pair<>(short2trans[t], syncLabelMap[(c % spCols) - transitions]));
					assert trace[startTraceAt + l] == syncLabelMap[(c % spCols) - transitions];
					l++;
				} else {
					moves.push(new Pair<>((Transition) null, logMoveMap[(c % spCols) - transitions
							- synchronousTransitions]));
					assert trace[startTraceAt + l] == logMoveMap[(c % spCols) - transitions - synchronousTransitions];
					l++;
				}
			}
		}
		return l;
	}
}
