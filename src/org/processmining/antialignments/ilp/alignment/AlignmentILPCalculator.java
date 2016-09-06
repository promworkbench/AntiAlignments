package org.processmining.antialignments.ilp.alignment;

import gnu.trove.list.TIntList;
import gnu.trove.list.TShortList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import lpsolve.LpSolve;
import nl.tue.astar.util.LPMatrix;
import nl.tue.astar.util.LPMatrix.LPMatrixException;

import org.deckfour.xes.classification.XEventClass;
import org.processmining.antialignments.ilp.AbstractILPCalculator;
import org.processmining.framework.plugin.Progress;
import org.processmining.framework.util.Pair;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;

public class AlignmentILPCalculator extends AbstractILPCalculator {

	static double EPSILON = 0.001;

	private static final int NOMOVE = 0x0000FFFF;

	// number of columns in synchronous product matrix
	private int spCols;
	// number of rows in synchronous product matrix
	private int spRows;
	// number of synchronous transitions in synchronous product
	private int synchronousTransitions;

	// maps from [0..synchronousTransitions> to the actual transition
	// that is synchronous at that position
	private short[] syncTransitionMap;
	// maps from [0..synchronousTransitions> to the label of the event
	// that is synchronous at that position
	private short[] syncLabelMap;
	// maps from [0..synchronousTransitions> to the actual event
	// that is synchronous at that position
	private short[] syncEventMap;

	// represents the window of the considered trace
	private short[] traceWindow;

	// presents a sorting on the labels from [0..labels>, such that the synchronous
	// labels are first, then the rest
	private final short[] label2pos;

	private int minEvents = 2;

	private int solve;

	private int setup;

	private int alignmentCosts;

	protected int steps;

	private final int[] modelMoveCost;

	private final int[] logMoveCost;

	private final int[] syncMoveCost;

	public AlignmentILPCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<XEventClass> label2short, TShortObjectMap<XEventClass> short2label,
			TransEvClassMapping mapping, short[][] log, Map<Transition, Integer> mapTrans2Cost,
			Map<XEventClass, Integer> mapEvClass2Cost, Map<Transition, Integer> mapSync2Cost) {
		super(net, initialMarking, finalMarking, label2short, short2label, mapping, log);

		label2pos = new short[labels];

		modelMoveCost = new int[transitions];
		syncMoveCost = new int[transitions];
		for (short i = 0; i < transitions; i++) {
			modelMoveCost[i] = mapTrans2Cost.get(short2trans[i]);
			Integer sc = mapSync2Cost.get(short2trans[i]);
			syncMoveCost[i] = sc == null ? 0 : sc;
		}
		logMoveCost = new int[labels];
		for (short i = 0; i < labels; i++) {
			logMoveCost[i] = mapEvClass2Cost.get(short2label.get(i));
		}
		//		try {
		//			gbEnv = new GRBEnv();
		//			gbEnv.set(GRB.IntParam.OutputFlag, 0);
		//			mode = MODE_GUROBI;
		//		} catch (GRBException _) {
		//			mode = MODE_LPSOLVE;
		//		}
		//
	}

	public TIntList getAlignment(Progress progress, Marking initialMarking, Marking finalMarking, int trace)
			throws LPMatrixException {

		return solveSequential(progress, initialMarking, finalMarking, log[trace]);

	}

	public TIntList getAlignmentWithoutTrace(Progress progress, Marking initialMarking, Marking finalMarking)
			throws LPMatrixException {

		return solveSequential(progress, initialMarking, finalMarking, new short[0]);

	}

	public Pair<Transition, Short> toPair(int move) {
		short t = (short) (move >>> 16);
		short l = (short) (move & 0x0000FFFF);
		Transition trans = t == (short) NOMOVE ? null : short2trans[t];
		Short label = l == (short) NOMOVE ? null : l;

		return new Pair<>(trans, label);
	}

	private void prepareSyncProduct(short[] trace, int maxLengthX, int startTraceAt) {
		Arrays.fill(label2pos, (short) -1);

		traceWindow = new short[Math.min(trace.length - startTraceAt, maxLengthX)];

		TShortList trs = new TShortArrayList();
		TShortList evts = new TShortArrayList();
		TShortList indices = new TShortArrayList();
		short pos = 0;
		for (short e = (short) startTraceAt; e < trace.length && e < startTraceAt + maxLengthX; e++) {
			for (short t = invisibleTransitions; t < transitions; t++) {
				if (equalLabel(t, trace[e])) {
					if (label2pos[trace[e]] < 0) {
						label2pos[trace[e]] = pos;
						pos++;
					}
					trs.add(t);
					evts.add(trace[e]);
					indices.add((short) (e - startTraceAt));
				}
			}
			traceWindow[e - startTraceAt] = trace[e];
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

		synchronousTransitions = syncLabelMap.length;
		spCols = transitions + synchronousTransitions + traceWindow.length;
		spRows = places + traceWindow.length;

	}

	protected LPMatrix<?> setupLpForHybrid(int maxLengthX, int minEvent, boolean integerVariables,
			Marking initialMarking, Marking finalMarking, short[] trace, int startTraceAt) {

		// Find the transitions in the net labeled with an event in 
		// log[traceToConsider][startTraceAt] .. log[traceToConsider][startTraceAt+lengthX-1]

		// A synchronous transition is needed for each (t,e) combination where e is an event

		prepareSyncProduct(trace, maxLengthX, startTraceAt);

		LPMatrix<?> lp;
		//-  
		//-                  AB   >= m0
		//-                  AAB  >= m0
		//-                  AAAA == mf        
		int rows = spRows * maxLengthX + (places + labels) + //
				// Xi - Xi+1<=0        Xi <= 1  
				(maxLengthX - 1) + maxLengthX + 1 + 1 + 1 + 1;
		//  X part is small           Y part is full
		int cols = spCols * maxLengthX + 2 * transitions - invisibleTransitions + labels;

		lp = setupMatrix(rows, cols);

		// row for the alignment costs

		int acRow = lp.getNrows() - 1;
		int progressRow = lp.getNrows() - 2;

		lp.setConstrType(acRow, LPMatrix.GE);

		int labelCountRowOffset = maxLengthX * spRows + places;

		// Fist set up all the sub-matrices for the exact part
		for (int block = 0; block < maxLengthX; block++) {
			int firstRow = block * spRows;
			int firstEventRow = block * spRows + places;
			// Setup the constraint types
			for (int p = 0; p < places; p++) {
				if (NAMES) {
					lp.setRowName(firstRow + p, "r" + block + "P" + p);
				}
				lp.setConstrType(firstRow + p, LPMatrix.GE);
			}
			for (int e = 0; e < traceWindow.length; e++) {
				if (NAMES) {
					lp.setRowName(firstEventRow + e, "r" + block + "T" + traceWindow[e]);
				}
				lp.setConstrType(firstEventRow + e, LPMatrix.GE);
			}
			if (NAMES) {
				lp.setRowName(firstEventRow + traceWindow.length, "r" + block + "Tf");
			}
			lp.setConstrType(firstEventRow + traceWindow.length, LPMatrix.GE);
			for (int min = 0; min < block; min++) {
				// First the whole incidence matrix
				matrixA.copyIntoMatrix(lp, firstRow, min * spCols);
				// then the SyncMove part
				for (int sm = 0; sm < synchronousTransitions; sm++) {
					int col = min * spCols + transitions + sm;
					matrixA.copyColumnIntoMatrix(syncTransitionMap[sm], lp, block * spRows, col);
					// token flow on event
					lp.setMat(firstEventRow + syncEventMap[sm], col, -1);
					if (syncEventMap[sm] < synchronousTransitions - 1) {
						lp.setMat(firstEventRow + syncEventMap[sm] + 1, col, 1);
					}
				}
				// then the LogMove part
				for (int e = 0; e < traceWindow.length; e++) {
					int col = min * spCols + transitions + synchronousTransitions + e;
					// token flow on event
					lp.setMat(firstEventRow + e, col, -1);
					if (e < traceWindow.length - 1) {
						lp.setMat(firstEventRow + e + 1, col, 1);
					}

				}
			}
			// Second, Aminus
			matrixAMin.copyIntoMatrix(lp, block * spRows, block * spCols);
			matrixA.copyIntoMatrix(lp, maxLengthX * spRows, block * spCols);
			for (short t = 0; t < transitions; t++) {
				int col = block * spCols + t;
				// Set Objective
				lp.setObjective(col, getCostForModelMove(t) + (minEvent == 0 ? 0 : EPSILON));
				lp.setMat(progressRow, col, 1);

				lp.setBinary(col, integerVariables);
				lp.setUpbo(col, 1);
				lp.setMat(acRow, col, getCostForModelMove(t));
				// Add labels
				if (NAMES) {
					lp.setColName(col, "c" + block + "M" + t);
				}
			}

			// then the SyncMove part
			for (int sm = 0; sm < synchronousTransitions; sm++) {
				int col = block * spCols + transitions + sm;
				matrixAMin.copyColumnIntoMatrix(syncTransitionMap[sm], lp, block * spRows, col);
				matrixA.copyColumnIntoMatrix(syncTransitionMap[sm], lp, maxLengthX * spRows, col);

				// token flow on event
				lp.setMat(firstEventRow + syncEventMap[sm], col, -1);

				// label count for syncmove in X
				lp.setMat(labelCountRowOffset + label2pos[syncLabelMap[sm]], col, 1);

				// Set Objective
				lp.setObjective(col, getCostForSync(syncTransitionMap[sm], syncLabelMap[sm]));
				lp.setMat(acRow, col, getCostForSync(syncTransitionMap[sm], syncLabelMap[sm]));
				lp.setMat(progressRow, col, 1);
				lp.setBinary(col, integerVariables);
				lp.setUpbo(col, 1);

				// Add labels
				if (NAMES) {
					lp.setColName(col, "c" + block + "S" + syncTransitionMap[sm]);
				}
			}
			// then the LogMove part
			for (int e = 0; e < traceWindow.length; e++) {
				int col = block * spCols + transitions + synchronousTransitions + e;
				// token flow on event
				lp.setMat(firstEventRow + e, col, -1);

				// label count for logmove in X
				lp.setMat(labelCountRowOffset + label2pos[traceWindow[e]], col, 1);

				// Set Objective
				lp.setObjective(col, getCostForLogMove(traceWindow[e]));
				lp.setMat(acRow, col, getCostForLogMove(traceWindow[e]));
				lp.setMat(progressRow, col, 1);
				lp.setBinary(col, integerVariables);
				lp.setUpbo(col, 1);

				// Add labels
				if (NAMES) {
					lp.setColName(col, "c" + block + "L" + traceWindow[e]);
				}
			}
		}

		// That concludes the X part, i.e. the sync product part of the incidence matrix
		// Now the Y part
		// First the Model Moves
		matrixA.copyIntoMatrix(lp, maxLengthX * spRows, maxLengthX * spCols);
		for (int p = 0; p < places; p++) {
			if (NAMES) {
				lp.setRowName(maxLengthX * spRows + p, "rP" + p);
			}
			// set Constraint type
			lp.setConstrType(maxLengthX * spRows + p, LPMatrix.EQ);
		}
		for (short t = 0; t < transitions; t++) {
			int col = maxLengthX * spCols + t;
			// Set Objective
			if (t < invisibleTransitions) {
				lp.setObjective(col, getCostForModelMove(t) + 0.001);
			} else {
				lp.setObjective(col, getCostForModelMove(t) + (minEvent == 0 ? EPSILON : 0));
			}

			lp.setMat(acRow, col, getCostForModelMove(t));

			// count progress
			lp.setMat(progressRow, col, 1);

			lp.setInt(col, integerVariables);
			// Add labels
			if (NAMES) {
				lp.setColName(col, "cM" + t);
			}
		}

		// Second the synchronous moves
		matrixA.copyIntoMatrixFromColumn(invisibleTransitions, lp, maxLengthX * spRows, maxLengthX * spCols
				+ transitions);
		for (short t = invisibleTransitions; t < transitions; t++) {
			int col = maxLengthX * spCols + transitions - invisibleTransitions + t;
			// Set Objective
			lp.setObjective(col, getCostForSync(t, trans2label[t]) + 2 * EPSILON);
			lp.setMat(acRow, col, getCostForSync(t, trans2label[t]));

			// count progress
			lp.setMat(progressRow, col, 1);
			if (trans2label[t] >= 0) {
				lp.setInt(col, integerVariables);
			}

			// Add labels
			if (NAMES) {
				lp.setColName(col, "cS" + t);
			}

			// label count for logmove in X
			lp.setMat(labelCountRowOffset + label2pos[trans2label[t]], col, 1);
		}

		// Third the log moves
		for (short l = 0; l < labels; l++) {
			int col = maxLengthX * spCols + 2 * transitions - invisibleTransitions + label2pos[l];
			int row = labelCountRowOffset + label2pos[l];

			// label count for logmove in X
			lp.setMat(row, col, 1);

			// Set Objective
			lp.setObjective(col, getCostForLogMove(l) + EPSILON);
			lp.setMat(acRow, col, getCostForLogMove(l));

			// count progress
			lp.setMat(progressRow, col, 1);
			lp.setInt(col, integerVariables);

			// set Constraint type
			lp.setConstrType(row, LPMatrix.EQ);

			// Add labels
			if (NAMES) {
				lp.setRowName(row, "rL" + l);
				lp.setColName(col, "cL" + l);
			}
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
			if (NAMES) {
				lp.setRowName(row + i, "X" + i + "_X" + (i + 1));
			}
			lp.setConstrType(row + i, LpSolve.GE);
		}
		row += maxLengthX - 1;

		int constant = Math.max(trace.length - startTraceAt, transitions); //TODO:??
		// C.lastX - Y >= 0
		for (int t = (maxLengthX - 1) * spCols + invisibleTransitions; t < maxLengthX * spCols; t++) {
			lp.setMat(row, t, constant);
		}
		for (int t = maxLengthX * spCols + invisibleTransitions; t < maxLengthX * spCols + 2 * transitions
				- invisibleTransitions; t++) {
			lp.setMat(row, t, -1);
		}
		if (NAMES) {
			lp.setRowName(row, "X" + (maxLengthX - 1) + "_Y");
		}
		lp.setConstrType(row, LpSolve.GE);

		row++;

		//	Events progress required
		for (int t = 0; t < maxLengthX * spCols; t++) {
			if ((t % spCols >= transitions)) {
				// add up sync of log moves
				lp.setMat(row, t, 1);
			}
		}
		if (NAMES) {
			lp.setRowName(row, "EVTS_" + minEvent);
		}
		lp.setConstrType(row, LpSolve.GE);
		row++;

		// sos constraints, sum max 1
		for (int i = 0; i < maxLengthX; i++) {
			for (int t = i * spCols + invisibleTransitions; t < (i + 1) * spCols; t++) {
				lp.setMat(row + i, t, 1);
			}
			if (NAMES) {
				lp.setRowName(row + i, "X" + i + ".1");
			}
			lp.setConstrType(row + i, LpSolve.LE);
		}

		//		row += maxLengthX;
		//		for (int t = 0; t < lp.getNcolumns(); t++) {
		//			lp.setMat(row, t, 1);
		//		}
		if (NAMES) {
			lp.setRowName(progressRow, "SUM");
		}
		lp.setConstrType(progressRow, LPMatrix.GE);

		lp.setMinim();
		//		lp.setMaxim();

		double[] rhs = new double[lp.getNrows()];

		// first -initial Marking
		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row < spRows * maxLengthX + places; row += spRows) {
				rhs[row] = -initialMarking.occurrences(p);
			}
			if (trace.length > startTraceAt) {
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
		for (int e = startTraceAt; e < trace.length; e++) {
			rhs[row + label2pos[trace[e]]]++;
		}

		// then X.i - X.i+1 >= 0
		// then ensure progress >= 0

		row = maxLengthX * spRows + places + labels + maxLengthX;
		rhs[row] = minEvent;

		// then X.1 <= 1
		for (row = row + 1; row < rhs.length - 2; row++) {
			rhs[row] = 1;
		}

		// do something!
		rhs[rhs.length - 2] = 1;

		lp.setRhVec(rhs);

		return lp;
	}

	private double getCostForModelMove(short t) {
		return modelMoveCost[t];
	}

	private double getCostForLogMove(short label) {
		return logMoveCost[label];
	}

	private double getCostForSync(short t, short s) {
		return syncMoveCost[t];
	}

	public void doExperiment(Marking initialMarking, Marking finalMarking) throws LPMatrixException {
		VERBOSE = false;
		NAMES = false;

		for (mode = 1; mode < 3; mode++) {
			FileWriter writer;
			try {
				writer = new FileWriter("D:/temp/antialignment/testResults"
						+ (mode == MODE_GUROBI ? "Gurobi" : "LpSolve") + ".csv");
				writer.write("cutOffLength;minEvents;trace;setuptime;solvetime;totalTime;alignmentCosts;isFS;isTrace\r\n");
				writer.flush();
			} catch (IOException e1) {
				return;
			}

			for (int traceToConsider = 0; traceToConsider < log.length; traceToConsider++) {
				for (cutOffLength = 1; cutOffLength < log[traceToConsider].length; cutOffLength++) {
					for (minEvents = 1; minEvents <= cutOffLength; minEvents++) {
						long start = System.currentTimeMillis();
						TIntList moves = solveSequential(null, initialMarking, finalMarking, log[traceToConsider]);
						long end = System.currentTimeMillis();

						try {
							writer.write(cutOffLength + ";" + minEvents + ";" + traceToConsider + ";" + setup + ";"
									+ solve + ";" + (end - start) + ";" + alignmentCosts + ";"
									+ checkAndReorderFiringSequence(moves, initialMarking, finalMarking, false) + ";"
									+ checkTrace(moves, log[traceToConsider], true) + "\r\n");
							writer.flush();

						} catch (IOException e) {
						}

					}
				}
			}
			try {
				writer.close();
			} catch (IOException e) {

			}
		}

	}

	protected TIntList solveSequential(Progress progress, Marking initialMarking, Marking finalMarking, short[] trace)
			throws LPMatrixException {

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		LPMatrix<?> matrix;

		int startTraceAt = 0;
		solve = 0;
		setup = 0;
		alignmentCosts = 0;
		steps = 0;

		TIntList moves = new TIntArrayList();
		Marking marking = initialMarking;

		int maxLengthX = cutOffLength;
		if (VERBOSE) {
			System.out.println("---------------------------- Starting trace -----------------------------------------");
		}
		do {
			if (VERBOSE) {
				System.out.println("Trying to get to " + finalMarking
						+ //
						" starting with  " + maxLengthX + " exact steps with [" + minEvents + ".." + maxLengthX
						+ "] events explained, starting from index " + startTraceAt + "/" + trace.length + ".");
			}

			long start = System.currentTimeMillis();
			matrix = setupLpForHybrid(maxLengthX, Math.min(minEvents, trace.length - startTraceAt), true, marking,
					finalMarking, trace, startTraceAt);

			if (VERBOSE) {
				System.out
						.println("ILP size: " + matrix.getNrows() + " rows and " + matrix.getNcolumns() + " columns.");
			}

			// Compute the new marking
			double[] vars = new double[matrix.getNcolumns()];
			long mid = System.currentTimeMillis();
			int result = matrix.solve(vars);
			steps++;

			if (result == LPMatrix.OPTIMAL) {

				// Update the alignment costs
				double costX = matrix.product(vars, 0, spCols * maxLengthX, matrix.getNrows() - 1);
				//				double costY = matrix.product(vars, spCols * maxLengthX, vars.length, matrix.getNrows() - 1);
				alignmentCosts += costX;

				// Compute the number of explained events
				int l = updateListOfMoves(maxLengthX, vars, trace, startTraceAt, moves, false);
				startTraceAt += l;

				//				// Remove trailing model moves if progress was made
				//				int removed = 0;
				//				int move = moves.get(moves.size() - 1);
				//				int varsX = spCols * maxLengthX;
				//				while (l > 0 && ((move & NOMOVE) == NOMOVE)) {
				//					// trailing model move
				//					short t = (short) (move >>> 16);
				//					if (VERBOSE) {
				//						System.out.println("Removing: " + toPair(move));
				//					}
				//					if (trans2label[t] >= 0) {
				//						// visible transition:
				//						removed++;
				//					}
				//					alignmentCosts -= getCostForModelMove(t);
				//					moves.removeAt(moves.size() - 1);
				//					move = moves.get(moves.size() - 1);
				//					for (int i = varsX; i-- > varsX - spCols;) {
				//						vars[i] = 0;
				//					}
				//					vars[spCols * (maxLengthX - 1 - removed) + t] = 0;
				//				}
				//				if (VERBOSE && removed > 0) {
				//					System.out.println(removed + " trailing model moves removed.");
				//				}

				// Get the intermediate Marking reached after the cutOffLength
				Marking reachedMarkingX = getIntermediateMarking(marking, matrix, vars, false);
				if (VERBOSE) {
					System.out.println("T: " + printTrace(trace));
					System.out.print("A: ");
					printMoves(moves, trace);
					System.out.println();
					assert checkAndReorderFiringSequence(moves, initialMarking, reachedMarkingX, false);
					assert checkTrace(moves, trace, false);
				}
				if (l == 0 && reachedMarkingX.equals(marking)) {
					// reachedMarkingX is equal to start marking. COuld be a loop, but maybe termination is
					// reached with invisible steps in Y.
					Marking reachedMarkingY = getIntermediateMarking(reachedMarkingX, matrix, vars, true);
					if (!reachedMarkingY.equals(finalMarking)) {
						// Stuck in a model move loop. Force the solver out
						// Terminate the model...
						maxLengthX++;
						if (VERBOSE) {
							System.out.println("Local loop reached, increasing maxLengthX to " + maxLengthX);
						}
					} else {
						// add invisible steps to the list of moves and update alignmentCostsF
						updateListOfMoves(maxLengthX, vars, trace, startTraceAt, moves, true);
						alignmentCosts += matrix.product(vars, spCols * maxLengthX, spCols * maxLengthX
								+ invisibleTransitions, matrix.getNrows() - 1);

						reachedMarkingX = reachedMarkingY;
						if (VERBOSE) {
							System.out.println("Finishing the model with tau-steps from Y.");
						}
					}
				} else {
					maxLengthX = cutOffLength;
				}
				if (VERBOSE) {
					assert checkAndReorderFiringSequence(moves, initialMarking, reachedMarkingX, false);
				}
				marking = reachedMarkingX;

				if (marking.equals(finalMarking)) {
					// Copy the remaining log Moves.
					for (int e = startTraceAt; e < trace.length; e++) {
						moves.add(NOMOVE << 16 | trace[e]);
					}

				}

				if (VERBOSE) {
					System.out.println("Alignment costs so far: " + alignmentCosts);
					System.out.println("Reached intermediate marking with objective with costs " + costX);
					System.out.println("Marking: " + marking);
					System.out.println("Vars: " + Arrays.toString(vars));
					System.out.print("Partial alignment :");
					printMoves(moves, trace);

				}
			} else {
				// BACKTRACKING NEEDED
				// So, the current marking has an infeasible answer, i.e. we cannot actually get 
				// to the final marking. 
				//
				// Backtrack and increase minEvents at the first step. This will
				// ensure that we don't loop.
				System.err.println("Infeasible model, waiting for input...");
				try {
					System.in.read();
					return null;
				} catch (IOException e) {
				}

				FileWriter writer;
				try {
					//					((LpSolve) matrix.toSolver()).writeLp("D:/temp/antialignment/debugLP-Alignment.lp");
					//					((LpSolve) matrix.toSolver()).writeMps("D:/temp/antialignment/debugLP-Alignment.mps");
					writer = new FileWriter("D:/temp/antialignment/debugLP-Alignment.csv");
					matrix.printLp(writer, ";");
					writer.close();

				} catch (Exception e1) {
					e1.printStackTrace();
					return null;
				}
				break;
			}
			long end = System.currentTimeMillis();
			setup += mid - start;
			solve += end - mid;
		} while (!marking.equals(finalMarking) && (progress == null || !progress.isCancelled()));

		if (VERBOSE) {
			System.out.println();
			System.out.println("Alignment done:");
			printMoves(moves, trace);
			System.out.println("Trace: " + printTrace(trace));
			assert checkAndReorderFiringSequence(moves, initialMarking, finalMarking, false);
			assert checkTrace(moves, trace, true);
			System.out.println("Is a firing sequence: "
					+ checkAndReorderFiringSequence(moves, initialMarking, finalMarking, false));
			System.out.println("Is the trace: " + checkTrace(moves, trace, true));
			System.out.println("Alignment costs: " + alignmentCosts);
			System.out.println("Setup time: " + setup);
			System.out.println("Solve time: " + solve);
			System.out.println();

		}
		return moves;
	}

	private String printTrace(short[] trace) {
		if (trace == null)
			return "null";
		int iMax = trace.length - 1;
		if (iMax == -1)
			return "[]";

		StringBuilder b = new StringBuilder();
		b.append('[');
		for (int i = 0;; i++) {
			b.append(short2label.get(trace[i]));
			if (i == iMax)
				return b.append(']').toString();
			b.append(", ");
		}
	}

	protected Marking getIntermediateMarking(Marking marking, LPMatrix<?> matrix, double[] vars, boolean invisibleY) {
		Marking newMarking = new Marking(marking);
		for (short p = 0; p < places; p++) {
			// compute the effect of the x vectors on place p.
			double v = matrix.product(vars, invisibleY ? cutOffLength * spCols : 0, cutOffLength * spCols
					+ (invisibleY ? invisibleTransitions : 0), cutOffLength * spRows + p);
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

	public void printMoves(TIntList moves, short[] trace) {
		double cost = 0;
		int mmSeq = 0;
		int lmSeq = 0;
		int smSeq = 0;
		int mmSeqMax = 0;
		int lmSeqMax = 0;
		int smSeqMax = 0;

		for (int i = 0; i < moves.size(); i++) {
			if (mmSeq > mmSeqMax) {
				mmSeqMax = mmSeq;
			}
			if (lmSeq > lmSeqMax) {
				lmSeqMax = lmSeq;
			}
			if (smSeq > smSeqMax) {
				smSeqMax = smSeq;
			}
			short t = (short) (moves.get(i) >>> 16);
			short l = (short) (moves.get(i) & 0x0000FFFF);

			System.out.print("[M:");
			System.out.print(t == (short) NOMOVE ? ">>" : short2trans[t].getLabel());
			System.out.print(",L:");
			System.out.print(l == (short) NOMOVE ? ">>" : short2label.get(l));
			System.out.print("],");
			if (t == (short) NOMOVE) {
				cost += getCostForLogMove(l);
				mmSeq = 0;
				lmSeq++;
				smSeq = 0;
			} else if (l == (short) NOMOVE) {
				cost += getCostForModelMove(t);
				mmSeq++;
				lmSeq = 0;
				smSeq = 0;
			} else {
				cost += getCostForSync(t, l);
				mmSeq = 0;
				lmSeq = 0;
				smSeq++;
			}
		}
		if (mmSeq > mmSeqMax) {
			mmSeqMax = mmSeq;
		}
		if (lmSeq > lmSeqMax) {
			lmSeqMax = lmSeq;
		}
		if (smSeq > smSeqMax) {
			smSeqMax = smSeq;
		}
		System.out.println();
		System.out.println("Cost of alignment: " + cost);
		System.out.println("Max sequences: MM: " + mmSeqMax + ", LM: " + lmSeqMax + ", SM: " + smSeqMax);
	}

	private int updateListOfMoves(int maxLengthX, double[] vars, short[] trace, int startTraceAt, TIntList moves,
			boolean invisibleY) {
		int l = 0;
		for (int c = invisibleY ? maxLengthX * spCols : 0; c < maxLengthX * spCols
				+ (invisibleY ? invisibleTransitions : 0); c++) {
			if (vars[c] > 0.5 || c % spCols < invisibleTransitions) {
				if (c % spCols < transitions) {
					int v = (int) (vars[c] + 0.5);
					while (v > 0) {
						//						moves.push(new Pair<>(short2trans[c % spCols], Short.MIN_VALUE));
						moves.add(((c % spCols) << 16) | NOMOVE);
						v--;
					}
				} else if (c % spCols < transitions + synchronousTransitions) {
					assert (int) (vars[c] + 0.5) == 1;
					int t = syncTransitionMap[(c % spCols) - transitions];
					//					moves.push(new Pair<>(short2trans[t], syncLabelMap[(c % spCols) - transitions]));

					moves.add((t << 16) | syncLabelMap[(c % spCols) - transitions]);
					assert trace[startTraceAt + l] == syncLabelMap[(c % spCols) - transitions];
					l++;
				} else {
					//					moves.push(new Pair<>((Transition) null, traceWindow[(c % spCols) - transitions
					//							- synchronousTransitions]));
					moves.add(NOMOVE << 16 | traceWindow[(c % spCols) - transitions - synchronousTransitions]);
					assert trace[startTraceAt + l] == traceWindow[(c % spCols) - transitions - synchronousTransitions];
					l++;
				}
			}
		}
		return l;
	}

	public boolean checkAndReorderFiringSequence(TIntList moves, Marking initialMarking, Marking finalMarking,
			boolean mergeMoves) {

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		short[] modelMoveStack = new short[moves.size()];
		short[] modelMoveLocationStack = new short[moves.size()];
		int mPos = -1;

		short[] logMoveStack = new short[moves.size()];
		short[] logMoveLocationStack = new short[moves.size()];
		int lPos = -1;

		int checked = -1;
		boolean ok = true;
		for (short t_i = 0; ok && t_i < moves.size(); t_i++) {
			int move = moves.get(t_i);
			short e = (short) (move & 0x0000FFFF);
			short t = (short) (move >>> 16);
			if (t == (short) NOMOVE) {
				if (mPos >= 0 && modelMoveStack[mPos] == e) {
					// last model move can be merged with this log move
					int m = moves.get(modelMoveLocationStack[mPos]) & 0xFFFF0000;
					m = m | (move & 0x0000FFFF);
					moves.set(t_i, m);
					System.out.println("Moving model move: " + modelMoveLocationStack[mPos] + " forward. Merges with "
							+ t_i + " Checked:" + checked);
					moves.removeAt(modelMoveLocationStack[mPos]);
					// checked--, since array was shortened by one
					if (modelMoveLocationStack[mPos] <= checked) {
						checked = t_i - 1;
					}
					// t-- means: go to next move
					// t -=2 means: process this move again
					// t -=3 means: process previous move.
					mPos--;
					t_i -= 3;
				} else if (mergeMoves) {
					// push on logMoveStack
					lPos++;
					logMoveStack[lPos] = e;
					logMoveLocationStack[lPos] = t_i;
				}
				// logMove
				continue;
			}

			// LogMove handled

			if (e == (short) NOMOVE && trans2label[t] >= 0) {
				// Visible model move
				if (lPos >= 0 && logMoveStack[lPos] == trans2label[t]) {
					// there is a log move on the stack with this label
					// merge them
					int m = moves.get(logMoveLocationStack[lPos]) & 0x0000FFFF;
					m = m | (move & 0xFFFF0000);
					moves.set(t_i, m);
					System.out.println("Moving log move: " + logMoveLocationStack[lPos] + " forward. Merges with "
							+ t_i + " Checked:" + checked);
					moves.removeAt(logMoveLocationStack[lPos]);
					checked = logMoveLocationStack[lPos] - 1;

					lPos--;
					t_i -= 3;
				} else if (mergeMoves) {
					mPos++;
					modelMoveStack[mPos] = trans2label[t];
					modelMoveLocationStack[mPos] = t_i;
				}
			} else if (trans2label[t] >= 0) {
				// empty the modelmove stack
				mPos = -1;
				lPos = -1;
			}

			if (t_i > checked) {
				checked = t_i;
				try {
					semantics.executeExecutableTransition(short2trans[t]);
				} catch (IllegalTransitionException _) {
					// so this transition was not enabled.
					//				assert (t.isInvisible());
					ok &= short2trans[t].isInvisible();
					assert ok;
					// push forward to first visible transition
					int j;
					for (j = t_i + 1; j < moves.size(); j++) {
						Transition tj = short2trans[moves.get(j) >>> 16];
						if (tj != null && tj.isInvisible()) {
							moves.set(j - 1, moves.get(j));
						} else {
							moves.set(j - 1, move);
							break;
						}
					}
					if (j == moves.size()) {
						moves.set(j - 1, move);
					}
					t_i--;
					continue;
				}
			}
		}

		return ok && semantics.getCurrentState().equals(finalMarking);
	}

	public boolean checkTrace(TIntList moves, short[] trace, boolean ensureFullTrace) {
		int i = 0;
		boolean ok = true;
		for (int t_i = 0; ok && t_i < moves.size() && i < trace.length; t_i++) {
			short e = (short) (moves.get(t_i) & 0x0000FFFF);
			if (e == (short) NOMOVE) {
				continue;
			}
			//			/assert trace[i] == e;
			ok &= trace[i] == e;
			i++;
		}

		return ok && (!ensureFullTrace || i == trace.length);
	}

	public double getCost(Transition t, XEventClass label) {
		if (t == null) {
			return getCostForLogMove(label2short.get(label));
		} else if (label == null) {
			return getCostForModelMove(trans2short.get(t));
		} else {
			return getCostForSync(trans2short.get(t), label2short.get(label));
		}
	}

	public double getCost(TIntList moves, short[] trace) {
		double cost = 0;

		for (int i = 0; i < moves.size(); i++) {

			short t = (short) (moves.get(i) >>> 16);
			short l = (short) (moves.get(i) & 0x0000FFFF);

			if (t == (short) NOMOVE) {
				cost += getCostForLogMove(l);
			} else if (l == (short) NOMOVE) {
				cost += getCostForModelMove(t);
			} else {
				cost += getCostForSync(t, l);
			}
		}

		return cost;
	}

	public void setMinEvents(int minEvents) {
		this.minEvents = minEvents;
	}

	public int getMinEvents() {
		return minEvents;
	}

}
