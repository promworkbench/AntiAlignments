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

import nl.tue.astar.util.ilp.LPMatrix;
import nl.tue.astar.util.ilp.LPMatrixException;

import org.deckfour.xes.classification.XEventClass;
import org.processmining.antialignments.ilp.AbstractILPCalculator;
import org.processmining.antialignments.ilp.util.HybridEquationResult;
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

	private static final short NOMOVE = -1;

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

	private int alignmentCostRow;

	private int progressRow;

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
		short t = splitT(move);
		short l = splitL(move);
		Transition trans = t == NOMOVE ? null : short2trans[t];
		Short label = l == NOMOVE ? null : l;

		return new Pair<>(trans, label);
	}

	private void prepareSyncProduct(short[] trace, int maxLengthX, int startTraceAt) {
		Arrays.fill(label2pos, NOMOVE);

		traceWindow = new short[Math.min(trace.length - startTraceAt, maxLengthX)];

		TShortList trs = new TShortArrayList();
		TShortList evts = new TShortArrayList();
		TShortList indices = new TShortArrayList();
		short pos = 0;
		for (short e = (short) startTraceAt; e < trace.length && e < startTraceAt + traceWindow.length; e++) {
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
			Marking initialMarking, Marking finalMarking, short[] trace, int startTraceAt) throws LPMatrixException {

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

		alignmentCostRow = lp.getNrows() - 1;
		progressRow = lp.getNrows() - 2;

		lp.setConstrType(alignmentCostRow, LPMatrix.GE);

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
					if (syncEventMap[sm] < traceWindow.length - 1) {
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
				lp.setObjective(col, getCostForModelMove(t));
				lp.setMat(alignmentCostRow, col, getCostForModelMove(t));

				lp.setMat(progressRow, col, 1);

				lp.setBinary(col, integerVariables);
				lp.setUpbo(col, 1);
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
				lp.setMat(alignmentCostRow, col, getCostForSync(syncTransitionMap[sm], syncLabelMap[sm]));
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
				lp.setMat(alignmentCostRow, col, getCostForLogMove(traceWindow[e]));
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
				lp.setObjective(col, getCostForModelMove(t));
			}

			lp.setMat(alignmentCostRow, col, getCostForModelMove(t));

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
			lp.setObjective(col, getCostForSync(t, trans2label[t]));
			lp.setMat(alignmentCostRow, col, getCostForSync(t, trans2label[t]));

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
			lp.setObjective(col, getCostForLogMove(l));
			lp.setMat(alignmentCostRow, col, getCostForLogMove(l));

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
			lp.setConstrType(row + i, LPMatrix.GE);
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
		lp.setConstrType(row, LPMatrix.GE);

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
		lp.setConstrType(row, LPMatrix.GE);
		row++;

		// sos constraints, sum max 1
		for (int i = 0; i < maxLengthX; i++) {
			for (int t = i * spCols + invisibleTransitions; t < (i + 1) * spCols; t++) {
				lp.setMat(row + i, t, 1);
			}
			if (NAMES) {
				lp.setRowName(row + i, "X" + i + ".1");
			}
			lp.setConstrType(row + i, LPMatrix.LE);
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

	private double getCostForSync(short t, short l) {
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
		HybridEquationResult intermediate = new HybridEquationResult(initialMarking, 0, Integer.MAX_VALUE);

		int startTraceAt = 0;
		int maxLengthX = cutOffLength;
		int backTrackDepth = 0;
		int minEventsLocal = minEvents;

		solve = 0;
		setup = 0;
		alignmentCosts = 0;
		steps = 0;

		TIntList moves = new TIntArrayList();

		if (VERBOSE) {
			System.out.println("---------------------------- Starting trace -----------------------------------------");
		}
		do {
			if (startTraceAt < 0) {
				// So we backtracked a bit to far :)
				System.err.println("HMM");
			}

			if (VERBOSE) {
				System.out.println("Trying to get from " + intermediate.getMarking() + " to " + finalMarking
						+ //
						" starting with  " + maxLengthX + " exact steps with [" + minEventsLocal + ".." + maxLengthX
						+ "] events explained, starting from index " + startTraceAt + "/" + (trace.length - 1) + ".");
			}

			long start = System.currentTimeMillis();
			matrix = setupLpForHybrid(maxLengthX, Math.min(minEventsLocal, trace.length - startTraceAt), true,
					intermediate.getMarking(), finalMarking, trace, startTraceAt);

			if (VERBOSE) {
				System.out
						.println("ILP size: " + matrix.getNrows() + " rows and " + matrix.getNcolumns() + " columns.");
			}

			// Compute the new marking
			double[] vars = new double[matrix.getNcolumns()];
			long mid = System.currentTimeMillis();

			HybridEquationResult nextResult = determineSplitMarkingForHybrid(matrix, vars, trace, maxLengthX,
					backTrackDepth == 0);
			//			int result = matrix.solve(vars);
			steps++;

			if (nextResult == null) {
				if (VERBOSE) {
					System.err.println("Infeasibility detected, exception...");
				}
				try {
					matrix.printLpToCSV("d:/temp/antialignment/infeasibleLP.csv");
				} catch (IOException e) {

				}
				//TODO: EXCEPTION
				assert false;

			} else if (backTrackDepth < maxBackTrackDepth //- We can still backtrack
					&& intermediate.getLengthY() > 0 //- expected optimal match. Not found, hence does not exist.
					&& nextResult.getLengthY() > 0 //- there's something in Y, so no optimal solution found
					&& (nextResult.getLengthX() + nextResult.getLengthY() > //- over backtracking threshold
					/*                            */intermediate.getLengthY() * backtrackThreshold && startTraceAt > 0)) {

				// BACKTRACKING NEEDED
				// nextResult did not live up to its expectations. 
				// go back one step.
				if (VERBOSE) {
					System.out.println("Backtracking from " + nextResult.getMarking() + " since "
							+ (nextResult.getLengthX() + nextResult.getLengthY()) + " > " + backtrackThreshold
							* intermediate.getLengthY());
					assert checkAndReorderFiringSequence(moves, initialMarking, intermediate.getMarking(), false);
					assert checkTrace(moves, trace, false);
				}
				// undo last synchronous step in the alignment.
				int m_i = moves.size() - 1;
				int move = moves.get(m_i);
				short t = splitT(move);
				short l = splitL(move);
				if (t == NOMOVE) {
					// logMove
					alignmentCosts -= getCostForLogMove(l);
					startTraceAt--;
				} else if (l == NOMOVE) {
					// modelMove
					alignmentCosts -= getCostForModelMove(t);
					unfire(short2trans[t], intermediate.getMarking());
				} else {
					alignmentCosts -= getCostForSync(t, l);
					// sync move
					unfire(short2trans[t], intermediate.getMarking());
					startTraceAt--;
				}
				moves.removeAt(m_i);

				if (VERBOSE) {
					System.out.println("T: " + printTrace(trace));
					System.out.print("A: ");
					printMoves(moves, trace);
					System.out.println();
					assert checkAndReorderFiringSequence(moves, initialMarking, intermediate.getMarking(), false);
					assert checkTrace(moves, trace, false);
				}

				// increase the maxLength of X by 1.
				// startTraceAt has been adjusted
				// TODO: How to guarantee progress, i.e. how to not get stuck in an infinite loop.
				maxLengthX++;
				minEventsLocal++;
				backTrackDepth++;

			} else {
				// accept the solution.

				// Update the alignment costs
				//				double costX = matrix.product(vars, 0, spCols * maxLengthX, alignmentCostRow);
				//				double costY = matrix.product(vars, spCols * maxLengthX, vars.length, matrix.getNrows() - 1);
				alignmentCosts += nextResult.getLengthX();
				// Compute the number of explained events
				int l = updateListOfMoves(maxLengthX, vars, trace, startTraceAt, moves, backTrackDepth == 0);

				// Get the intermediate Marking reached after the cutOffLength
				nextResult
						.setMarking(getIntermediateMarking(intermediate.getMarking(), maxLengthX, matrix, vars, false));

				startTraceAt += l;

				// Get the intermediate Marking reached after the cutOffLength
				Marking reachedMarkingX = nextResult.getMarking();//getIntermediateMarking(marking, matrix, vars, false);
				if (VERBOSE) {
					System.out.println("T: " + printTrace(trace));
					System.out.print("A: ");
					printMoves(moves, trace);
					System.out.println();
					assert checkAndReorderFiringSequence(moves, initialMarking, reachedMarkingX, false);
					assert checkTrace(moves, trace, false);
				}
				//				if (l == 0 && reachedMarkingX.equals(intermediate.getMarking())) {
				//					// reachedMarkingX is equal to start marking. COuld be a loop, but maybe termination is
				//					// reached with invisible steps in Y.
				//					Marking reachedMarkingY = getIntermediateMarking(reachedMarkingX, matrix, vars, true);
				//					if (!reachedMarkingY.equals(finalMarking)) {
				//						// Stuck in a model move loop. Force the solver out
				//						// Terminate the model...
				//						maxLengthX++;
				//						if (VERBOSE) {
				//							System.out.println("Local loop reached, increasing maxLengthX to " + maxLengthX);
				//						}
				//						
				//					} else {
				//						// add invisible steps to the list of moves and update alignmentCostsF
				//						updateMovesWithInvisiblesFromY(maxLengthX, vars, moves);
				//						alignmentCosts += matrix.product(vars, spCols * maxLengthX, spCols * maxLengthX
				//								+ invisibleTransitions, matrix.getNrows() - 1);
				//
				//						reachedMarkingX = reachedMarkingY;
				//						if (VERBOSE) {
				//							System.out.println("Finishing the model with tau-steps from Y.");
				//						}
				//					}
				//				} else {
				//				}
				if (VERBOSE) {
					assert checkAndReorderFiringSequence(moves, initialMarking, reachedMarkingX, false);
				}
				maxLengthX = cutOffLength;
				minEventsLocal = minEvents;
				intermediate = nextResult;
				//				marking = reachedMarkingX;

				if (intermediate.getMarking().equals(finalMarking)) {
					// Copy the remaining log Moves.
					for (int e = startTraceAt; e < trace.length; e++) {
						moves.add(NOMOVE << 16 | trace[e]);
						alignmentCosts += getCostForLogMove(trace[e]);
					}

				}

				if (VERBOSE) {
					System.out.println("Alignment costs so far: " + alignmentCosts);
					System.out.println("Reached intermediate marking with objective with costs "
							+ nextResult.getLengthX());
					System.out.println("Marking: " + intermediate.getMarking());
					System.out.println("Vars: " + Arrays.toString(vars));
					System.out.print("Partial alignment :");
					printMoves(moves, trace);

				}
				intermediate = nextResult;
				backTrackDepth = 0;
			}
			long end = System.currentTimeMillis();
			setup += mid - start;
			solve += end - mid;
		} while (!intermediate.getMarking().equals(finalMarking) && (progress == null || !progress.isCancelled()));

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

	protected HybridEquationResult determineSplitMarkingForHybrid(LPMatrix<?> matrix, double[] vars, short[] trace,
			int maxLengthX, boolean removeTrailingModelMoves) throws LPMatrixException {

		int result = matrix.solve(vars);

		if (result == LPMatrix.OPTIMAL) {

			int costX = 0;
			int costTrailModelMoves = 0;
			boolean onlyModelMove = true;
			// Update the alignment costs
			for (int c = 0; c < spCols * maxLengthX; c++) {
				if (vars[c] > 0) {
					int cost_c = (int) (vars[c] * matrix.getMat(alignmentCostRow, c) + 0.5);
					costX += cost_c;
					if (c % spCols < transitions) {
						// Trailing modelMove, to be removed later
						costTrailModelMoves += cost_c;
					} else {
						// log or sync move, reset trailing model move cost
						costTrailModelMoves = 0;
						onlyModelMove = false;
					}
				}
			}
			int costY = (int) (matrix.product(vars, spCols * maxLengthX, vars.length, alignmentCostRow) + 0.5);
			if (removeTrailingModelMoves && !onlyModelMove) {
				// remove the cost of trailing moves if there were other moves and move them from X to Y.
				costX -= costTrailModelMoves;
				costY += costTrailModelMoves;
			}

			return new HybridEquationResult(costX, costY);
		} else {
			return null;
		}
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

	protected Marking getIntermediateMarking(Marking marking, int maxLengthX, LPMatrix<?> matrix, double[] vars,
			boolean invisibleY) {
		Marking newMarking = new Marking(marking);
		int[] effect = new int[places];
		for (short p = 0; p < places; p++) {
			// compute the effect of the x vectors on place p.
			double v = matrix.product(vars, invisibleY ? maxLengthX * spCols : 0, maxLengthX * spCols
					+ (invisibleY ? invisibleTransitions : 0), maxLengthX * spRows + p);
			if (v < 0) {
				effect[p] += (int) (v - 0.5);
			} else if (v > 0) {
				newMarking.add(short2place[p], (int) (v + .5));
			}
		}
		for (short p = 0; p < places; p++) {
			while (effect[p] < 0) {
				newMarking.remove(short2place[p]);
				effect[p]++;
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
			short t = splitT(moves.get(i));
			short l = splitL(moves.get(i));

			System.out.print("[M:");
			System.out.print(t == NOMOVE ? ">>" : short2trans[t].getLabel());
			System.out.print(",L:");
			System.out.print(l == NOMOVE ? ">>" : short2label.get(l));
			System.out.print("],");
			if (t == NOMOVE) {
				cost += getCostForLogMove(l);
				mmSeq = 0;
				lmSeq++;
				smSeq = 0;
			} else if (l == NOMOVE) {
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
		assert cost == alignmentCosts;
		System.out.println("Max sequences: MM: " + mmSeqMax + ", LM: " + lmSeqMax + ", SM: " + smSeqMax);
	}

	private short[] splitTE(int move) {
		return new short[] { splitT(move), splitL(move) };

	}

	private short splitT(int move) {
		return (short) (move >>> 16);

	}

	private short splitL(int move) {
		return (short) (move & 0x0000FFFF);

	}

	private int updateListOfMoves(int maxLengthX, double[] vars, short[] trace, int startTraceAt, TIntList moves,
			boolean ignoreTrailingModelMoves) {
		int l = 0;
		TIntList tmpList = new TIntArrayList();
		TIntList tmpVars = new TIntArrayList();

		for (int b = 0; b < maxLengthX; b++) {
			int bl = 0;
			for (int c = b * spCols; c < (b + 1) * spCols; c++) {
				if (vars[c] > 0.5) {// || c % spCols < invisibleTransitions) {
					tmpVars.add(c);
					if (c % spCols < transitions) {
						int v = (int) (vars[c] + 0.5);
						while (v > 0) {
							//						moves.push(new Pair<>(short2trans[c % spCols], Short.MIN_VALUE));
							tmpList.add(((c % spCols) << 16) | 0x0000FFFF);
							v--;
						}
					} else if (c % spCols < transitions + synchronousTransitions) {
						assert (int) (vars[c] + 0.5) == 1;
						int t = syncTransitionMap[(c % spCols) - transitions];
						//					moves.push(new Pair<>(short2trans[t], syncLabelMap[(c % spCols) - transitions]));

						tmpList.add((t << 16) | syncLabelMap[(c % spCols) - transitions]);
						assert trace[startTraceAt + l] == syncLabelMap[(c % spCols) - transitions];
						bl++;
					} else {
						//					moves.push(new Pair<>((Transition) null, traceWindow[(c % spCols) - transitions
						//							- synchronousTransitions]));
						tmpList.add(0xFFFF0000 | traceWindow[(c % spCols) - transitions - synchronousTransitions]);
						assert trace[startTraceAt + l] == traceWindow[(c % spCols) - transitions
								- synchronousTransitions];
						bl++;
					}
				}
			}
			l += bl;
			if (!ignoreTrailingModelMoves || bl > 0) {
				// we are not ignoring trailing model moves and we found progress in the events
				// purge the tempList into Moves
				moves.addAll(tmpList);
				tmpList.clear();
				tmpVars.clear();
			}
		}
		if (l == 0) {
			// no progress at all, so purge the templist regardless of whether we ignore trailing
			// model moves
			moves.addAll(tmpList);
		} else if (ignoreTrailingModelMoves) {
			// remove trailing model moves by removing them from the variable set.
			if (VERBOSE) {
				System.out.println("Removing trailing model moves: " + tmpVars);
			}
			for (int i = 0; i < tmpVars.size(); i++) {
				vars[tmpVars.get(i)] = 0;
			}
		}
		return l;
	}

	private void updateMovesWithInvisiblesFromY(int maxLengthX, double[] vars, TIntList moves) {

		for (int c = (maxLengthX + 1) * spCols; c < (maxLengthX + 1) * spCols + invisibleTransitions; c++) {
			int v = (int) (vars[c] + 0.5);
			while (v > 0) {
				//						moves.push(new Pair<>(short2trans[c % spCols], Short.MIN_VALUE));
				moves.add(((c % spCols) << 16) | 0x0000FFFF);
				v--;
			}
		}
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
			short e = splitL(move);
			short t = splitT(move);
			if (t == NOMOVE) {
				if (mPos >= 0 && modelMoveStack[mPos] == e) {
					// last model move can be merged with this log move
					int m = moves.get(modelMoveLocationStack[mPos]) & 0xFFFF0000;
					m = m | (move & 0x0000FFFF);
					moves.set(t_i, m);
					if (VERBOSE) {
						System.out.println("Moving model move: " + modelMoveLocationStack[mPos]
								+ " forward. Merges with " + t_i + " Checked:" + checked);
					}
					moves.removeAt(modelMoveLocationStack[mPos]);
					// checked--, since array was shortened by one
					if (modelMoveLocationStack[mPos] <= checked) {
						if (checked <= t_i - 1) {
							// Set checked to the new location
							checked = t_i - 1;
						} else {
							// we remove one element
							checked--;
						}
					}
					// t-- means: go to next move
					// t -=2 means: process this move again
					// t -=3 means: process previous move.
					mPos--;
					if (t_i > 1) {
						t_i -= 3;
					} else {
						t_i = -1;
					}
				} else if (mergeMoves) {
					// push on logMoveStack
					lPos++;
					logMoveStack[lPos] = e;
					logMoveLocationStack[lPos] = t_i;
				}
				// logMove
				continue;
			}

			// There is a transition to fire. Check the firing before trying the merge.
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
						Transition tj = short2trans[splitT(moves.get(j))];
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
				}
			}

			// LogMove and firing handled
			if (e == NOMOVE && trans2label[t] >= 0) {
				// Visible model move
				if (lPos >= 0 && logMoveStack[lPos] == trans2label[t]) {
					// there is a log move on the stack with this label
					// merge them
					int m = moves.get(logMoveLocationStack[lPos]) & 0x0000FFFF;
					m = m | (move & 0xFFFF0000);
					moves.set(t_i, m);
					if (VERBOSE) {
						System.out.println("Moving log move: " + logMoveLocationStack[lPos] + " forward. Merges with "
								+ t_i + " Checked:" + checked);
					}
					moves.removeAt(logMoveLocationStack[lPos]);
					checked--;

					lPos--;
					if (t_i > 1) {
						t_i -= 3;
					} else {
						t_i = -1;
					}
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

		}

		return ok && semantics.getCurrentState().equals(finalMarking);
	}

	public boolean checkTrace(TIntList moves, short[] trace, boolean ensureFullTrace) {
		int i = 0;
		boolean ok = true;
		for (int t_i = 0; ok && t_i < moves.size() && i < trace.length; t_i++) {
			short e = (short) (moves.get(t_i) & 0x0000FFFF);
			if (e == NOMOVE) {
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

			short t = splitT(moves.get(i));
			short l = splitL(moves.get(i));

			if (t == NOMOVE) {
				cost += getCostForLogMove(l);
			} else if (l == NOMOVE) {
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

	public boolean isGurobi() {
		return mode == MODE_GUROBI;
	}

	public boolean isLpSolve() {
		return mode == MODE_LPSOLVE;
	}

}
