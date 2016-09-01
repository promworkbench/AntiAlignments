package org.processmining.antialignments.algorithm;

import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.util.Vector;

import lpsolve.LpSolve;

import org.processmining.antialignments.algorithm.LPMatrix.LPMatrixException;
import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetEdge;
import org.processmining.models.graphbased.directed.petrinet.PetrinetNode;
import org.processmining.models.graphbased.directed.petrinet.elements.Arc;
import org.processmining.models.graphbased.directed.petrinet.elements.InhibitorArc;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class AntiAlignmentILPCalculator {

	public class HybridEquationResult {
		private final Marking marking;
		private final int lengthX;
		private final int lengthY;
		private final TShortList antiAlignment;
		private final Stack<Transition> firingSequence;

		public HybridEquationResult(Marking marking, int lengthX, int lengthY, TShortList antiAlignment,
				Stack<Transition> firingSequence) {
			this.marking = marking;
			this.lengthX = lengthX;
			this.lengthY = lengthY;
			this.antiAlignment = antiAlignment;
			this.firingSequence = firingSequence;
		}

		public Marking getMarking() {
			return marking;
		}

		public int getLengthX() {
			return lengthX;
		}

		public int getLengthY() {
			return lengthY;
		}

		public boolean equals(Object o) {
			if (o != null && o instanceof HybridEquationResult) {
				HybridEquationResult h = (HybridEquationResult) o;
				return h.lengthX == lengthX && h.lengthY == lengthY && h.getMarking().equals(marking);
			}
			return false;
		}

		public int hashCode() {
			return marking.hashCode() + 31 * lengthX + 31 * 31 * lengthY;
		}

		public String toString() {
			return marking.toString() + " Lx: " + lengthX + " Ly: " + lengthY;
		}

		public TShortList getAntiAlignment() {
			return antiAlignment;
		}

		public Stack<Transition> getFiringSequence() {
			return firingSequence;
		}
	}

	public static boolean VERBOSE = true;

	private static final int MODE_LPSOLVE = 1;
	private static final int MODE_GUROBI = 2;
	private int mode;

	private final Petrinet net;
	//	private final Marking initialMarking;
	private final PetrinetSemantics semantics;
	//	private final Marking finalMarking;
	private final TObjectShortMap<String> label2short;
	private TShortObjectMap<String> short2label;
	private short transitions;
	private short places;
	private Transition[] short2trans;
	private Place[] short2place;

	private short[] trans2label;
	private final short[][] log;
	private final int maxFactor;

	private GRBEnv gbEnv;
	private TObjectShortHashMap<Object> trans2int;
	private TObjectShortHashMap<Object> place2int;
	private final int maxLength;
	private short invisibleTransitions;

	private int cutOffLength = 5;

	public AntiAlignmentILPCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short, TShortObjectMap<String> short2label, short[][] log, int maxLength,
			int maxFactor) {
		this.net = net;
		//		this.initialMarking = initialMarking;
		//		this.finalMarking = finalMarking;
		this.label2short = label2short;
		this.short2label = short2label;
		this.log = log;
		this.maxLength = maxLength;
		this.maxFactor = maxFactor;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

		setupDataStructures();

		mode = MODE_LPSOLVE;
		try {
			gbEnv = new GRBEnv();
			gbEnv.set(GRB.IntParam.OutputFlag, 0);
			mode = MODE_GUROBI;
		} catch (GRBException e) {
			mode = MODE_LPSOLVE;
		}

	}

	protected void solveSequential(int maxLength, Marking initialMarking, Marking finalMarking,
			final Vector<Transition> firingSequence, final TShortList antiAlignment, final int traceToIgnore)
			throws LPMatrixException {

		LPMatrix<?> matrix;
		HybridEquationResult intermediate = new HybridEquationResult(initialMarking, 0, 0, new TShortArrayList(1),
				new Stack<Transition>());
		int startTraceAt = 0;
		int maxLengthY = maxLength - cutOffLength;
		int maxLengthX = cutOffLength;
		if (maxLength < cutOffLength) {
			maxLengthY = 0;
			maxLengthX = maxLength;
		}
		List<HybridEquationResult> splits = new ArrayList<HybridEquationResult>();
		Marking marking = intermediate.getMarking();
		splits.add(intermediate);

		do {
			if (VERBOSE) {
				System.out.println("Trying to get from " + marking + " to " + finalMarking + " in " + maxLengthX
						+ " steps and " + maxLengthY + " approximate steps.");
			}

			if (startTraceAt < 0) {
				// So we backtracked a bit to far :)
				System.err.println("HMM");
			}

			matrix = setupLpForHybrid(maxLengthX, maxLengthY, true, marking, finalMarking, traceToIgnore, startTraceAt);
			// Determine an intermediate marking at lengthX visible steps
			HybridEquationResult nextResult = determineSplitMarkingForHybrid(matrix, maxLengthX, maxLengthY == 0);

			if (nextResult == null) {
				// The model could not be solved, i.e. there is no path from intermediateMarking to the finalMarking in
				// no more than maxLengthX+maxLengthY steps.
				// This is where the fun begins....
				//				((LpSolve) matrix.toSolver()).printLp();
				System.err.println("HMM");
			}
			if (VERBOSE) {
				System.out.println("Marking reached: " + nextResult.getMarking() + " in " + nextResult.getLengthX()
						+ " steps and estimating " + nextResult.getLengthY() + " remaining steps.");
			}
			if (nextResult.lengthX + nextResult.lengthY < intermediate.getLengthY()) {
				// nextResult did not live up to its expectations. 
				// go back one step.
				if (VERBOSE) {
					System.out.println("Backtracking since " + (nextResult.lengthX + nextResult.lengthY) + " < "
							+ intermediate.getLengthY());
				}
				// remove the last transition from the stack.
				Transition last;
				intermediate.getAntiAlignment().removeAt(intermediate.getAntiAlignment().size() - 1);

				do {
					last = intermediate.getFiringSequence().pop();

					for (PetrinetEdge<?, ?> e : net.getInEdges(last)) {
						if (e instanceof Arc) {
							Arc a = (Arc) e;
							Place p = (Place) a.getSource();
							int w = a.getWeight();
							intermediate.getMarking().add(p, w);
						}
					}
					for (PetrinetEdge<?, ?> e : net.getOutEdges(last)) {
						if (e instanceof Arc) {
							Arc a = (Arc) e;
							Place p = (Place) a.getTarget();
							int w = a.getWeight();
							while (w > 0) {
								intermediate.getMarking().remove(p);
								w--;
							}
						}
					}
					// and remove all invisible transitions before, as they belong to the last 
					// X transition.
				} while (!intermediate.getFiringSequence().isEmpty()
						&& intermediate.getFiringSequence().peek().isInvisible());

				// We've removed the last transition in the row, so allow for one more 
				// transition in Y
				maxLengthY++;
				// and include one additional event.
				startTraceAt--;
			} else {
				// Setup for next iteration
				intermediate = nextResult;
				splits.add(intermediate);
				marking = intermediate.getMarking();
				maxLengthY -= maxLengthX;
				startTraceAt += maxLengthX;
			}
		} while (maxLengthY >= 0 && !intermediate.getMarking().equals(finalMarking));

		if (!intermediate.getMarking().equals(finalMarking)) {
			splits.add(intermediate);
			matrix = setupLpForHybrid(maxLengthX + maxLengthY, 0, true, marking, finalMarking, traceToIgnore,
					startTraceAt);
			intermediate = determineSplitMarkingForHybrid(matrix, maxLengthX + maxLengthY, true);
			splits.add(intermediate);
		}

		assert intermediate.getMarking().equals(finalMarking);

		for (HybridEquationResult r : splits) {
			antiAlignment.addAll(r.getAntiAlignment());
			firingSequence.addAll(r.getFiringSequence());
		}

		if (VERBOSE) {
			System.out.println("Anti alignment: " + toString(antiAlignment.toArray()));
			System.out.println("Dlog: " + getMinimalHammingDistanceToLog(antiAlignment.toArray(), log, traceToIgnore));
			System.out.println("Firing Sequence: " + firingSequence);

			System.out.println(splits);
		}
		//		System.out.println("Test");

	}

	private boolean checkFiringSequence(Vector<Transition> firingSequence, Marking initialMarking, Marking finalMarking) {

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		for (int t_i = 0; t_i < firingSequence.size(); t_i++) {
			Transition t = firingSequence.get(t_i);
			try {
				semantics.executeExecutableTransition(t);
			} catch (IllegalTransitionException e) {
				// so this transition was not enabled.
				assert (t.isInvisible());
				// push forward to first visible transition
				int j;
				for (j = t_i + 1; j < firingSequence.size(); j++) {
					if (firingSequence.get(j).isInvisible()) {
						firingSequence.set(j - 1, firingSequence.get(j));
					} else {
						firingSequence.set(j - 1, t);
						break;
					}
				}
				if (j == firingSequence.size()) {
					firingSequence.set(j - 1, t);
				}
				t_i--;
				continue;
			}
		}

		return semantics.getCurrentState().equals(finalMarking);
	}

	protected void solveByDrillingDown(final int splitInParts, int maxLength, Marking initialMarking,
			Marking finalMarking, final Vector<Transition> firingSequence, final TShortList antiAlignment,
			final int traceToIgnore, int startTracesAt) throws LPMatrixException {

		if (VERBOSE) {
			System.out.println("Solving maxlength " + maxLength + " from marking " + initialMarking + " to "
					+ finalMarking);
		}
		if (maxLength > cutOffLength) {

			int lengthX = maxLength / splitInParts + (maxLength % splitInParts > 0 ? 1 : 0);
			int maxLengthY = maxLength - lengthX;

			LPMatrix<?> matrix = setupLpForSplit(lengthX, maxLengthY, true, initialMarking, finalMarking,
					traceToIgnore, startTracesAt);
			double[] vars = new double[matrix.getNcolumns()];

			// Determine an intermediate marking at lengthX visible steps
			Marking intermediate = determineSplitMarking(matrix, traceToIgnore, vars);

			if (intermediate != null && !initialMarking.equals(intermediate)) {
				// Success!
				if (VERBOSE) {
					System.out.println("Intermediate marking: " + intermediate);
				}
				solveByDrillingDown(splitInParts, lengthX, initialMarking, intermediate, firingSequence, antiAlignment,
						traceToIgnore, startTracesAt);
				if (!finalMarking.equals(intermediate)) {
					solveByDrillingDown(splitInParts, maxLengthY, intermediate, finalMarking, firingSequence,
							antiAlignment, traceToIgnore, startTracesAt + lengthX);
				}

			} else {
				// Initial Marking is kept intact, or marking not reachabel in lengthX VISIBLE steps.
				// Possible spurious trace, or transition invariant detected. Now do 1+ rest search.
				lengthX = 1;
				maxLengthY = maxLength - 1;
				matrix = setupLpForSplit(lengthX, maxLengthY, true, initialMarking, finalMarking, traceToIgnore,
						startTracesAt);
				vars = new double[matrix.getNcolumns()];
				intermediate = determineSplitMarking(matrix, traceToIgnore, vars);

				if (intermediate != null) {
					// Success. We found an intermediate marking at 1 VISIBLE step (this can be added to the firing sequence, including
					// possible invisible steps.
					if (VERBOSE) {
						System.out.println("Intermediate marking: " + intermediate);
					}
					int length = firingSequence.size();
					for (int t = 0; t < transitions; t++) {
						if (vars[t] > 0.5) {
							assert ((int) (vars[t] + 0.5)) == 1 || short2trans[t].isInvisible();
							do {
								Transition trans = short2trans[t];
								if (trans.isInvisible()) {
									firingSequence.add(length, trans);
								} else {
									antiAlignment.add(label2short.get(trans.getLabel()));
									firingSequence.add(trans);
								}
								if (VERBOSE) {
									System.out.print(trans);
									System.out.print(",");
								}
								vars[t] = vars[t] - 1;
							} while (vars[t] > 0.5);
						}

					}
					if (VERBOSE) {
						System.out.println();
					}
					// Continue with the search from intermediate marking
					if (!intermediate.equals(finalMarking)) {
						solveByDrillingDown(splitInParts, maxLength - lengthX, intermediate, finalMarking,
								firingSequence, antiAlignment, traceToIgnore, startTracesAt + 1);
					}

				} else {
					// final marking not reachable from initialmarking in lengthX VISIBLE step followed by at most maxLengthY other visible steps. 
					// try invisible steps :)
					matrix = setupLpForFinalInvisibleSteps(true, initialMarking, finalMarking);
					vars = new double[matrix.getNcolumns()];
					int result = matrix.solve(vars);
					if (result != LPMatrix.OPTIMAL) {
						// Cannot reach the final marking with invisible steps. Error :)
						System.out.println("Cannot reach the final marking " + finalMarking + " from " + initialMarking
								+ "!");
					}
					for (int t = 0; t < invisibleTransitions; t++) {
						if (vars[t] > 0.5) {
							assert ((int) (vars[t] + 0.5)) == 1 || short2trans[t].isInvisible();
							do {
								Transition trans = short2trans[t];
								assert (trans.isInvisible());
								firingSequence.add(trans);
								if (VERBOSE) {
									System.out.print(trans);
									System.out.print(",");
								}
								vars[t] = vars[t] - 1;
							} while (vars[t] > 0.5);
						}

					}
					if (VERBOSE) {
						System.out.println();
					}
				}
			}

		} else {
			LPMatrix<?> matrix = setupLpForFullSequence(maxLength, true, initialMarking, finalMarking, traceToIgnore,
					startTracesAt);
			//			matrix.printLp();
			int length = 0;
			double[] vars = new double[matrix.getNcolumns()];
			solveForFullSequence(matrix, -1, maxLength, null, vars);
			for (int t = 0; t < vars.length - 2; t++) {
				if (t % transitions == 0) {
					length = firingSequence.size();
				}
				if (vars[t] > 0.5) {
					assert ((int) (vars[t] + 0.5)) == 1;
					Transition trans = short2trans[t % transitions];
					if (trans.isInvisible()) {
						firingSequence.add(length, trans);
					} else {
						antiAlignment.add(label2short.get(trans.getLabel()));
						firingSequence.add(trans);
					}
					if (VERBOSE) {
						System.out.print(trans);
						System.out.print(",");
					}

				}

			}
			if (VERBOSE) {
				System.out.println();
			}
		}

	}

	protected Marking determineSplitMarking(LPMatrix<?> matrix, int traceToIgnore, double[] vars)
			throws LPMatrixException {
		//			matrix.toLpSolve().printLp();

		//			matrix.printLp();

		int result = matrix.solve(vars);

		if (result == LPMatrix.OPTIMAL) {
			// get the intermediate marking:
			Marking intermediate = new Marking();
			// Iterate over the X vector
			for (int p = 0; p < places; p++) {
				double effect = -matrix.getRh(p + places);
				for (int t = 0; t < transitions; t++) {
					effect += vars[t] * matrix.getMat(p, t);
				}
				assert effect >= 0;
				if (effect > 0) {
					intermediate.add(short2place[p], (int) (effect + 0.5));
				}
			}
			return intermediate;
		} else {
			return null;
		}
	}

	public boolean doAntiAlignmentTest(Marking initialMarking, Marking finalMarking) throws LPMatrixException {
		VERBOSE = true;
		mode = MODE_LPSOLVE;
		mode = MODE_GUROBI;

		String sep = "\t";
		double exp = 1;
		cutOffLength = 4;
		int maxCutOffLength = 10;//maxLength * maxFactor

		System.out.println("mode" + sep + "cutOffLength" + sep + "time(ms)" + sep + "Dlog");
		while (cutOffLength < maxCutOffLength) {
			cutOffLength += 1;
			long sum = 0;
			long dist = 0;

			Vector<Transition> firingSequence;
			TShortList antiAlignment = null;

			for (int i = 0; i < exp; i++) {
				firingSequence = new Vector<Transition>(maxLength * maxFactor);
				antiAlignment = new TShortArrayList(maxLength * maxFactor);
				long start = System.currentTimeMillis();
				//				solveByDrillingDown(2, maxLength * maxFactor, initialMarking, finalMarking, firingSequence,
				//						antiAlignment, -1, 0);
				solveSequential(maxLength * maxFactor, initialMarking, finalMarking, firingSequence, antiAlignment, -1);
				long end = System.currentTimeMillis();
				assert checkFiringSequence(firingSequence, initialMarking, finalMarking);
				sum += end - start;

				short[] aa = antiAlignment.toArray();
				int hd = getMinimalHammingDistanceToLog(aa, log, -1);

				dist += hd;
			}

			System.out.println((mode == MODE_GUROBI ? "Gurobi" : "LpSolve") + sep + cutOffLength + sep + (sum / exp)
					+ sep + (dist / exp));

		}
		System.out.println("done");
		return true;
		//		System.exit(0);
	}

	public AntiAlignments getAntiAlignments(Marking initialMarking, Marking finalMarking) throws LPMatrixException {
		if (doAntiAlignmentTest(initialMarking, finalMarking)) {
			return null;
		}

		if (VERBOSE) {
			System.out.println("Solving by drilling down");
		}

		cutOffLength = 32;
		Vector<Transition> firingSequence = new Vector<Transition>(maxLength * maxFactor);
		TShortList antiAlignment = new TShortArrayList(maxLength * maxFactor);
		//		solveByDrillingDown(maxLength * maxFactor, initialMarking, finalMarking, firingSequence, antiAlignment, -1, 0);
		solveSequential(maxLength * maxFactor, initialMarking, finalMarking, firingSequence, antiAlignment, -1);
		// FiringSequence and anti-alignment are known, but the distances to the log need to be computed.
		short[] aa = antiAlignment.toArray();
		int hd = getMinimalHammingDistanceToLog(aa, log, -1);
		cutOffLength = 5;

		//		LPMatrix lpMatrix = setupLpForFullSequence(maxLength * maxFactor, true, initialMarking, finalMarking, 0);

		//		System.out.println("IP with " + lpMatrix.getNcolumns() + " columns and " + lpMatrix.getNrows() + " rows.");

		final AntiAlignments antiAlignments = new AntiAlignments(log.length);

		if (VERBOSE) {
			System.out.println("-----------------------------------------------");
			System.out.println("Whole log. ");
			System.out.println("Maxlength: " + maxFactor * maxLength);
			System.out.flush();
		}
		//		lp.printLp();
		//		double[] result = new double[lpMatrix.getNcolumns()];
		//		solveForFullSequence(lpMatrix, -1, maxFactor * maxLength, null, result);// "D:/temp/antialignment/ilp_instance_log.mps");

		antiAlignments.getMaxMinDistances()[log.length] = hd;
		antiAlignments.getAntiAlignments()[log.length] = aa;
		antiAlignments.getTraces()[log.length] = firingSequence;
		antiAlignments.getMaxDistances()[log.length] = Math.max(getMaxTraceLength(-1), aa.length);
		if (VERBOSE) {
			System.out.println("Anti alignment: " + toString(aa));
			System.out.println("Firing Sequence: " + firingSequence);
			System.out.println("Dlog: " + hd);
			System.out.flush();
		}
		//		int[] basis = new int[lp.getNrows() + lp.getNcolumns() + 1];
		//		lp.getBasis(basis, true);

		for (int t = 0; t < log.length; t++) {
			if (VERBOSE) {
				System.out.println("-----------------------------------------------" + (t + 1) + " / " + log.length);
				System.out.println("Removed Trace: " + toString(log[t]));
				System.out.println("Maxlength: " + maxFactor * log[t].length);
				System.out.flush();
			} //			lp.setBasis(basis, true);
				//			solveForFullSequence(lpMatrix, t, maxFactor * log[t].length, null, result);//, "D:/temp/antialignment/ilp_instance_t" + t + ".mps");

			firingSequence = new Vector<Transition>(maxLength * maxFactor);
			antiAlignment = new TShortArrayList(maxLength * maxFactor);
			//			solveByDrillingDown(maxFactor * log[t].length, initialMarking, finalMarking, firingSequence, antiAlignment,
			//					t, 0);
			solveSequential(maxFactor * log[t].length, initialMarking, finalMarking, firingSequence, antiAlignment, t);
			// FiringSequence and anti-alignment are known, but the distances to the log need to be computed.
			aa = antiAlignment.toArray();
			hd = getMinimalHammingDistanceToLog(aa, log, t);

			antiAlignments.getMaxMinDistances()[t] = hd;
			antiAlignments.getAntiAlignments()[t] = aa;
			antiAlignments.getMaxDistances()[t] = Math.max(getMaxTraceLength(t), aa.length);
			antiAlignments.getTraces()[t] = firingSequence;
			antiAlignments.getTraceDistances()[t] = getHammingDistanceToTrace(aa, log[t]);
			if (VERBOSE) {
				System.out.println("Anti alignment: " + toString(aa));
				System.out.println("Firing Sequence: " + firingSequence);
				System.out.println("Dlog: " + hd + "  Drem: " + antiAlignments.getTraceDistances()[t]);
				System.out.flush();
			}
		}
		if (VERBOSE) {
			System.out.println("-----------------------------------------------");
		}
		//		System.out.println(Arrays.toString(vars));

		return antiAlignments;

	}

	private int getMaxTraceLength(int traceToIgnore) {
		int max = 0;
		for (int t = 0; t < log.length; t++) {
			if (log[t].length > max && t != traceToIgnore) {
				max = log[t].length;
			}
		}
		return max;
	}

	protected int getMinimalHammingDistanceToLog(short[] antiAlignment, short[][] log, int traceToIgnore) {
		int hd = Integer.MAX_VALUE;
		for (int t = 0; t < log.length; t++) {
			if (t != traceToIgnore) {
				int hdt = getHammingDistanceToTrace(antiAlignment, log[t]);
				if (hdt < hd) {
					hd = hdt;
				}
				if (hdt == 0) {
					break;
				}
			}
		}
		return hd;
	}

	protected int getHammingDistanceToTrace(short[] antiAlignment, short[] trace) {
		int hdt = Math.abs(trace.length - antiAlignment.length);
		for (int e = 0; e < trace.length && e < antiAlignment.length; e++) {
			if (trace[e] != antiAlignment[e]) {
				hdt++;
			}
		}

		return hdt;
	}

	protected int solveForFullSequence(LPMatrix<?> lpMatrix, int traceToIgnore, int maxLength, String filename,
			double[] vars) throws LPMatrixException {

		int row;

		row = lpMatrix.getNrows() - log.length + traceToIgnore - 1;

		if (traceToIgnore >= 0) {
			// h doesn't matter
			lpMatrix.setMat(row, lpMatrix.getNcolumns() - 2, 0);
			// minus g for the removed trace.
			lpMatrix.setMat(row, lpMatrix.getNcolumns() - 1, -1);
			lpMatrix.setConstrType(row, LpSolve.EQ);
		}
		lpMatrix.setRh(lpMatrix.getNrows() - 1, maxLength);

		int result;
		try {
			result = lpMatrix.solve(vars);
		} finally {
			if (traceToIgnore >= 0) {
				// minimize h
				lpMatrix.setMat(row, lpMatrix.getNcolumns() - 2, -1);
				// ignore g.
				lpMatrix.setMat(row, lpMatrix.getNcolumns() - 1, 0);
				lpMatrix.setConstrType(row, LpSolve.GE);
			}
		}
		return result;

	}

	protected void setupDataStructures() {
		// replay log on model (or obtain existing replay result)
		transitions = 0;
		places = 0;

		short2trans = new Transition[net.getTransitions().size()];
		trans2int = new TObjectShortHashMap<>(net.getTransitions().size() / 2 * 3, 0.7f, (short) 0);
		for (Transition t : net.getTransitions()) {
			if (t.isInvisible()) {
				trans2int.put(t, transitions);
				short2trans[transitions] = t;
				transitions++;
			}
		}
		invisibleTransitions = transitions;
		for (Transition t : net.getTransitions()) {
			if (!t.isInvisible()) {
				trans2int.put(t, transitions);
				short2trans[transitions] = t;
				transitions++;
			}
		}
		short2place = new Place[net.getPlaces().size()];
		place2int = new TObjectShortHashMap<>(net.getPlaces().size() / 2 * 3, 0.7f, (short) 0);
		for (Place p : net.getPlaces()) {
			place2int.put(p, places);
			short2place[places] = p;
			places++;
		}

		trans2label = new short[transitions];

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short t;
			Transition trans;
			if (e instanceof Arc) {
				if (e.getSource() instanceof Place) {
					trans = (Transition) e.getTarget();
					t = trans2int.get(trans);
				} else {
					trans = (Transition) e.getSource();
					t = trans2int.get(trans);
				}
			} else if (e instanceof InhibitorArc) {
				trans = (Transition) e.getTarget();
				t = trans2int.get(trans);
			} else {
				continue;
			}
			if (trans.isInvisible()) {
				trans2label[t] = -1;
			} else {
				trans2label[t] = label2short.get(trans.getLabel());
			}
		}
	}

	protected LPMatrix<?> setupLpForFullSequence(int maxLength, boolean integerVariables, Marking initialMarking,
			Marking finalMarking, int traceToIgnore, int startTracesAt) {

		//		LpSolve lp = LpSolve.makeLp((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 2);
		LPMatrix<?> lp = setupMatrix((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength
				+ invisibleTransitions + 2);

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

			for (int block = 0; block <= maxLength; block++) {
				// First the whole A matrix.
				int r = block * places + p;
				lp.setConstrType(r, type);
				lp.setRowName(r, "A" + block + "_" + p);
				for (int c = t; c < block * transitions; c += transitions) {
					// update the A matrix
					lp.setMat(r, c, lp.getMat(r, c) + dir);
				}

				// Then, the  A- matrix.
				int c = block * transitions + t;
				if ((dir < 0 || trans.isInvisible()) && block < maxLength) {
					lp.setColName(c, trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
					lp.setInt(c, integerVariables);
					lp.setUpbo(c, 1.0);
					// update the A matrix only for consumption and for invisible transitions
					// or in the last block
					lp.setMat(r, c, lp.getMat(r, c) + dir);
				}

				// Then, in the last block
				if (block == maxLength) {
					lp.setConstrType(r, LpSolve.EQ);
					if (trans.isInvisible()) {
						// Also add invisible columns to the end.
						lp.setMat(r, c, lp.getMat(r, c) + dir);
						lp.setColName(c, trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
						lp.setInt(c, integerVariables);
						lp.setUpbo(c, 1.0);
					}
				}
			}
		}
		int row = places * (maxLength + 1);
		// set up the comparisons and sums
		for (int i = 0; i < maxLength - 1; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			for (int t = (i + 1) * transitions; t < (i + 2) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : -1);
			}
			lp.setRowName(row + i, "X" + i + "X" + (i + 1));
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// sos constraints, sum max 1
		row = places * (maxLength + 1) + maxLength - 1;
		for (int i = 0; i < maxLength; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			lp.setRowName(row + i, "X" + i + ".1");
			lp.setConstrType(row + i, LpSolve.LE);
		}

		row = places * (maxLength + 1) + 2 * maxLength - 1;
		for (int t = 0; t < log.length; t++) {
			// what's the hamming distance of the x's to trace t?
			for (int e = startTracesAt; e < startTracesAt + maxLength; e++) {
				for (short tr = 0; tr < trans2label.length; tr++) {
					if (e >= log[t].length && trans2label[tr] > 0) {
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 1);
					} //					if ((e >= log[t].length || trans2label[tr] != log[t][e]) && trans2label[tr] > 0) {
					else if (e < log[t].length && !equalLabel(tr, log[t][e]) && trans2label[tr] >= 0) {
						//						lp.setMat(row + t, e * transitions + tr + 1, 1);
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 0);
					} else if (trans2label[tr] < 0) {
						// tau steps do not count
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 0);
					} else {
						// based on the index and the fact that this transition matches that label,
						// set the value of this trace's row to 0;
						//						lp.setMat(row + t, e * transitions + tr + 1, 0);
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, -1);
					}
				}
			}
			lp.setRowName(row + t, "trace_" + t);
			if (t == traceToIgnore) {
				// minus g for the removed trace.
				lp.setMat(row + t, transitions * maxLength + invisibleTransitions + 1, -1);
				lp.setConstrType(row + t, LpSolve.EQ);
			} else {
				// minus h
				lp.setMat(row + t, transitions * maxLength + invisibleTransitions, -1);
				lp.setConstrType(row + t, LpSolve.GE);
			}
		}

		// Setup maxLength constraint
		row = (maxLength + 1) * places + 2 * maxLength + log.length - 1;
		for (int t = 0; t < transitions * maxLength + invisibleTransitions; t++) {
			lp.setMat(row, t, trans2label[t % transitions] < 0 ? 0 : 1);

			// Maximize the number of VISIBLE transition firings in the anti-alignment
			// Minimize the number of INVISIBLE transition firings in the anti-alignment
			lp.setObjective(t, trans2label[t % transitions] < 0 ? -.1 : .2);
		}
		lp.setRowName(row, "SumX");
		lp.setConstrType(row, LpSolve.LE);

		// maximize h = minimal distance to the log
		lp.setInt(transitions * maxLength + invisibleTransitions, integerVariables);
		lp.setObjective(transitions * maxLength + invisibleTransitions, maxLength * maxLength);
		// minimuze g = distance to the ignored trace
		lp.setInt(transitions * maxLength + invisibleTransitions + 1, integerVariables);
		lp.setObjective(transitions * maxLength + invisibleTransitions + 1, -maxLength);

		lp.setMaxim();

		double[] rhs = new double[(maxLength + 1) * places + 2 * maxLength + log.length];

		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row < places * (maxLength + 1); row += places) {
				rhs[row] = -initialMarking.occurrences(p);
			}
		}
		for (Place p : finalMarking.baseSet()) {
			row = maxLength * places + place2int.get(p);
			rhs[row] += finalMarking.occurrences(p);
		}
		for (row = (maxLength + 1) * places + maxLength; row < (maxLength + 1) * places + 2 * maxLength; row++) {
			rhs[row - 1] = 1;
		}

		row = lp.getNrows() - log.length - 1;
		for (int t = 0; t < log.length; t++) {
			if (startTracesAt > log[t].length) {
				rhs[row + t] = 0;
			} else {
				rhs[row + t] = -(log[t].length - startTracesAt);
			}
		}

		row = (maxLength + 1) * places + 2 * maxLength + log.length - 1;
		rhs[row] = maxLength;

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
	}

	private boolean equalLabel(short transition, short event) {
		return trans2label[transition] == event;
	}

	private String toString(short[] a) {
		if (a == null)
			return "null";
		int iMax = a.length - 1;
		if (iMax == -1)
			return "[]";

		StringBuilder b = new StringBuilder();
		b.append('[');
		for (int i = 0;; i++) {
			b.append(short2label.get(a[i]));
			if (i == iMax)
				return b.append(']').toString();
			b.append(", ");
		}
	}

	protected LPMatrix<?> setupLpForSplit(int maxLengthX, int maxLengthY, boolean integerVariables,
			Marking initialMarking, Marking finalMarking, int traceToIgnore, int startTracesAt) {

		//		LpSolve lp = LpSolve.makeLp((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 2);
		LPMatrix<?> lp = setupMatrix(2 * places + 2 * log.length + 3, transitions * 2 + 3);

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short p, t;
			int dir;

			//TODO: haNDLE INHIBITORS
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
			// set up all the two A matrixes.
			lp.setMat(p, t, lp.getMat(p, t) + dir);
			if (maxLengthX > 1 || trans.isInvisible() || dir < 0) {
				lp.setMat(p + places, t, lp.getMat(p + places, t) + dir);
			}
			lp.setMat(p, t + transitions, lp.getMat(p, t + transitions) + dir);

			lp.setInt(t, integerVariables);
			lp.setInt(t + transitions, integerVariables);

			lp.setColName(t, trans.getLabel().replace("\\n$invisible$", "") + "_X");
			lp.setColName(t + transitions, trans.getLabel().replace("\\n$invisible$", "") + "_Y");
			// maximize the Y vector length
			lp.setObjective(t + transitions, 1);

			lp.setRowName(p, (maxLengthX == 1 ? "B1" : "A1") + p);
			lp.setRowName(p + places, "A2" + p);

			lp.setConstrType(p, LpSolve.EQ);
			lp.setConstrType(p + places, LpSolve.GE);

			if (!trans.isInvisible()) {
				lp.setMat(2 * places, t, 1);
				lp.setMat(2 * places + 1, t + transitions, 1);

			} else {
				// where possible, transfer tau steps to Y instead of X
				lp.setObjective(transitions + t, 1);
			}

		}
		// The length of the X vector is equal to half the maxLength
		lp.setConstrType(2 * places, LpSolve.EQ);
		// The remainder is the Y vector
		lp.setConstrType(2 * places + 1, LpSolve.LE);

		int row = 2 * places + 2;

		// we consider only part of the log that is between startTracesAt and startTracesAt + maxLengthX for X
		for (int t = 0; t < log.length; t++) {
			// Set up the constraints to estimate Hamming Distance, first for X vectors
			for (int e = startTracesAt; e < startTracesAt + maxLengthX; e++) {
				if (e < log[t].length) {
					// There is an event at position e in trace t
					// Now find similarly labeled transitions
					for (short tr = 0; tr < trans2label.length; tr++) {
						if (equalLabel(tr, log[t][e])) {
							// add one to transition tr, labeled with log[t][e]
							lp.setMat(row + 2 * t, tr, lp.getMat(row + 2 * t, tr) + 1);
						}
					}
				}
			}
			lp.setRowName(row + 2 * t, "t_" + t + "_X");
			// minus h
			if (t != traceToIgnore) {
				lp.setMat(row + 2 * t, 2 * transitions, -1);
				lp.setConstrType(row + 2 * t, LPMatrix.LE);
			} else {
				lp.setConstrType(row + 2 * t, LPMatrix.GE);
			}
		}

		double maxOccInY = 0;
		// we consider only part of the log that is between maxLengthX and maxLengthY + maxLengthX for Y
		for (int t = 0; t < log.length; t++) {
			// Set up the constraints to estimate Hamming Distance, second for Y vectors
			for (int e = startTracesAt + maxLengthX; e < startTracesAt + maxLengthY + maxLengthX; e++) {
				if (e < log[t].length) {
					// There is an event at position e in trace t
					// Now find similarly labeled transitions
					for (short tr = 0; tr < trans2label.length; tr++) {
						if (equalLabel(tr, log[t][e])) {
							// add one to transition tr, labeled with log[t][e]
							double val = lp.getMat(row + 2 * t + 1, tr + transitions) + 1;
							lp.setMat(row + 2 * t + 1, tr + transitions, val);
							if (t != traceToIgnore && val > maxOccInY) {
								maxOccInY = val;
							}
						}
					}
				}
			}

			lp.setRowName(row + 2 * t + 1, "t_" + t + "_Y");
			// minus g
			if (t != traceToIgnore) {
				lp.setMat(row + 2 * t + 1, 2 * transitions + 1, -1);
				lp.setConstrType(row + 2 * t + 1, LPMatrix.LE);
			} else {
				lp.setConstrType(row + 2 * t + 1, LPMatrix.GE);
			}
		}

		// minimize distance to removed trace
		if (traceToIgnore >= 0 && traceToIgnore < log.length) {
			for (int t = 0; t < transitions; t++) {
				lp.setMat(2 * places + 2 * log.length + 2, t, lp.getMat(2 * places + 2 * traceToIgnore, t));
				lp.setMat(2 * places + 2 * log.length + 2, t + transitions,
						lp.getMat(2 * places + 2 * traceToIgnore + 1, t + transitions));
			}
			lp.setMat(2 * places + 2 * log.length + 2, 2 * transitions + 2, -1);
		} else {
			lp.setMat(2 * places + 2 * log.length + 2, 2 * transitions + 2, 1);
		}
		lp.setConstrType(2 * places + 2 * log.length + 2, LPMatrix.EQ);

		// minimize the distance to the removed trace by maximizing the overlap in labels
		lp.setObjective(2 * transitions + 2, +maxLengthY);

		// minimize h (HX is an underestimate of the minimal hamming distance from X to any trace)
		lp.setColName(2 * transitions, "HX");
		lp.setObjective(2 * transitions, -maxLengthX * maxLengthY);
		// minimize g (GX is an underestimate of the minimal hamming distance from Y to any trace)
		for (int t = transitions; t < 2 * transitions; t++) {
			if (trans2label[t - transitions] >= 0) {
				lp.setObjective(t, maxOccInY * maxOccInY);
			}
		}
		lp.setColName(2 * transitions + 1, "HY");
		lp.setObjective(2 * transitions + 1, -1);

		lp.setMaxim();

		double[] rhs = new double[lp.getNrows()];

		for (Place p : initialMarking.baseSet()) {
			rhs[place2int.get(p)] = -initialMarking.occurrences(p);
			rhs[place2int.get(p) + places] = -initialMarking.occurrences(p);
		}
		for (Place p : finalMarking.baseSet()) {
			rhs[place2int.get(p)] += finalMarking.occurrences(p);
		}
		rhs[2 * places] = maxLengthX;
		rhs[2 * places + 1] = maxLengthY;

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
	}

	private LPMatrix<?> setupMatrix(int rows, int columns) {
		if (mode == MODE_GUROBI) {
			return new LPMatrix.GUROBI(gbEnv, rows, columns);
		} else {
			return new LPMatrix.LPSOLVE(rows, columns);
		}
	}

	protected LPMatrix<?> setupLpForFinalInvisibleSteps(boolean integerVariables, Marking initialMarking,
			Marking finalMarking) {

		//		LpSolve lp = LpSolve.makeLp((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 2);
		LPMatrix<?> lp = setupMatrix(places, invisibleTransitions);

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short p, t;
			int dir;

			//TODO: Handle Inhibitors
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
			if (!trans.isInvisible()) {
				continue;
			}
			// set up the A matrixes.
			lp.setMat(p, t, lp.getMat(p, t) + dir);

			lp.setInt(t, integerVariables);

			lp.setColName(t, trans.getLabel().replace("\\n$invisible$", "") + "_X");

			// minimize the X vector length
			lp.setObjective(t, 1);

			lp.setRowName(p, "A" + p);

			lp.setConstrType(p, LpSolve.EQ);

		}

		lp.setMinim();

		double[] rhs = new double[places];

		for (Place p : initialMarking.baseSet()) {
			rhs[place2int.get(p)] = -initialMarking.occurrences(p);
		}
		for (Place p : finalMarking.baseSet()) {
			rhs[place2int.get(p)] += finalMarking.occurrences(p);
		}

		lp.setRhVec(rhs);

		return lp;
	}

	protected LPMatrix<?> setupLpForHybrid(int maxLengthX, int maxLengthY, boolean integerVariables,
			Marking initialMarking, Marking finalMarking, int traceToIgnore, int startTracesAt) {

		int useHY = 1; // set to 1 if not using HY

		//-                            B    >= m0
		//-                            AB   >= m0
		//-                            AAB  >= m0
		//-                            AAAA == mf             Xi - Xi+1<=0   Xi <= 1   Y.1 <= maxLength-lengthX 
		LPMatrix<?> lp = setupMatrix(places * (maxLengthX + 1) + (maxLengthX - 1) + maxLengthX + 1/*                     */
				//- y.1 - maxLength X_(lengthX-1).1 <=0 traceXHD     traceYHD
				+ 1 /*                             */+ useHY * log.length, transitions * maxLengthX + transitions + 3);

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

			for (int block = 0; block <= maxLengthX; block++) {
				// First the whole A matrices.
				int r = block * places + p;
				lp.setConstrType(r, type);
				lp.setRowName(r, "A" + block + "_" + p);
				for (int c = t; c < block * transitions; c += transitions) {
					// update the A matrix
					lp.setMat(r, c, lp.getMat(r, c) + dir);
				}

				// Then, the  A- matrix.
				int c = block * transitions + t;
				if ((dir < 0 || trans.isInvisible()) && block < maxLengthX) {
					lp.setColName(c, trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
					lp.setInt(c, integerVariables);
					lp.setUpbo(c, 1.0);
					// update the A matrix only for consumption and for invisible transitions
					// or in the last block
					lp.setMat(r, c, lp.getMat(r, c) + dir);
				}

				// Then, in the last block
				if (block == maxLengthX) {
					lp.setConstrType(r, LpSolve.EQ);
					lp.setMat(r, c, lp.getMat(r, c) + dir);
					lp.setColName(c, trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
					lp.setInt(c, integerVariables);
				}
			}
		}
		int row = places * (maxLengthX + 1);
		// set up the comparisons and sums
		for (int i = 0; i < maxLengthX - 1; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			for (int t = (i + 1) * transitions; t < (i + 2) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : -1);
			}
			lp.setRowName(row + i, "X" + i + "~X" + (i + 1));
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// sos constraints, sum max 1
		row += maxLengthX - 1;
		for (int i = 0; i < maxLengthX; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			lp.setRowName(row + i, "X" + i + ".1");
			lp.setConstrType(row + i, LpSolve.LE);
		}

		row += maxLengthX;
		// setup Y length
		for (int t = 0; t < transitions; t++) {
			lp.setMat(row, maxLengthX * transitions + t, trans2label[t] < 0 ? 0 : 1);
		}
		lp.setRowName(row, "Y.1");
		lp.setConstrType(row, LpSolve.LE);

		row++;
		// setup Y only if X_{lengthX}.1==1;
		for (int t = 0; t < transitions; t++) {
			lp.setMat(row, (maxLengthX - 1) * transitions + t, trans2label[t] < 0 ? 0 : -maxLengthY);
			lp.setMat(row, maxLengthX * transitions + t, trans2label[t] < 0 ? 0 : 1);
		}
		lp.setRowName(row, "Y~X_l");
		lp.setConstrType(row, LpSolve.LE);

		double maxOccInY = 0;
		//Hamming distances for the X and Y vectors
		row++;
		for (short tr = 0; tr < trans2label.length; tr++) {
			for (int t = 0; t < log.length; t++) {
				// what's the hamming distance of the x's to trace t?
				for (int e = startTracesAt; e < startTracesAt + maxLengthX; e++) {
					if (e >= log[t].length && trans2label[tr] > 0) {
						// e is passed the length of the trace
						lp.setMat(row + useHY * t, (e - startTracesAt) * transitions + tr, 1);
					} else if (e < log[t].length && !equalLabel(tr, log[t][e])) {
						// e is in the trace window, but label does not match
						lp.setMat(row + useHY * t, (e - startTracesAt) * transitions + tr, 0);
					} else {
						//e is in the trace window and label matches.
						lp.setMat(row + useHY * t, (e - startTracesAt) * transitions + tr, -1);
					}
				}

				if (useHY == 2) {
					// what's a guess for the hamming distance of the y's to trace t?
					for (int e = startTracesAt + maxLengthX; e < startTracesAt + maxLengthX + maxLengthY
							&& e < log[t].length; e++) {
						if (equalLabel(tr, log[t][e])) {
							// e is in the trace window and label matches
							double val = lp.getMat(row + useHY * t + 1, maxLengthX * transitions + tr) + 1;
							lp.setMat(row + useHY * t + 1, maxLengthX * transitions + tr, val);
							if (t != traceToIgnore && val > maxOccInY) {
								maxOccInY = val;
							}
						}
					}
				}

				if (tr == 0) {
					// Setup the rows for this trace;
					lp.setRowName(row + useHY * t, "Hx-t" + t);

					if (t == traceToIgnore) {
						// minus g for the removed trace.
						lp.setMat(row + useHY * t, transitions * (maxLengthX + 1) + 2, -1);
						lp.setConstrType(row + useHY * t, LpSolve.EQ);
					} else {
						// minus hx
						lp.setMat(row + useHY * t, transitions * (maxLengthX + 1), -1);
						// minus hy
						lp.setConstrType(row + useHY * t, LpSolve.GE);
					}
					if (useHY == 2) {
						lp.setRowName(row + useHY * t + 1, "Hy-t" + t);
						if (t == traceToIgnore) {
							// minus g for the removed trace.
							lp.setConstrType(row + useHY * t + 1, LpSolve.GE);
						} else {
							// minus hy
							lp.setMat(row + useHY * t + 1, transitions * (maxLengthX + 1) + 1, -1);
							lp.setConstrType(row + useHY * t + 1, LpSolve.GE);
						}
					}

				}

			}
		}

		// Setup objective function

		// Minimize all invisible transitions
		for (int t = 0; t < transitions * maxLengthX; t++) {
			// Minimize the number of INVISIBLE transition firings in the anti-alignment
			lp.setObjective(t, trans2label[t % transitions] < 0 ? -1 : 0);
		}
		for (int t = transitions * maxLengthX; t < transitions * (maxLengthX + 1); t++) {
			if (trans2label[t % transitions] >= 0) {
				// Maximize the number of Y transitions.
				lp.setObjective(t, maxOccInY * maxLengthY);
			} else {
				// Minimize the number of INVISIBLE transition firings in Y part of the anti-alignment
				lp.setObjective(t, -1.0 / ((maxLengthX + maxLengthY) * (maxLengthX + maxLengthY)));
			}
		}
		if (useHY == 1) {

			double max = 0.0;
			//		 Minimize overlap in Y vector
			for (short tr = 0; tr < trans2label.length; tr++) {
				for (int t = 0; t < log.length; t++) {
					for (int e = startTracesAt + maxLengthX; e < startTracesAt + maxLengthX + maxLengthY
							&& e < log[t].length; e++) {
						if (equalLabel(tr, log[t][e])) {
							// label at index e in trace t matched transition tr
							// increase the objective by 1 for transition tr in vector Y
							lp.setObjective(maxLengthX * transitions + tr,
									lp.getObjective(maxLengthX * transitions + tr) + 1);
						}
					}
				}
				if (lp.getObjective(maxLengthX * transitions + tr) > max) {
					max = lp.getObjective(maxLengthX * transitions + tr);
				}
			}
			max += 1.0;
			for (short tr = 0; tr < trans2label.length; tr++) {
				lp.setObjective(maxLengthX * transitions + tr, max - lp.getObjective(maxLengthX * transitions + tr));
			}
		} else {
			// minimize hy = similarity to y
			lp.setInt(transitions * (maxLengthX + 1) + 1, integerVariables);
			lp.setObjective(transitions * (maxLengthX + 1) + 1, -1);
		}

		// maximize hx = minimal distance to the log
		lp.setInt(transitions * (maxLengthX + 1), integerVariables);
		lp.setObjective(transitions * (maxLengthX + 1), (maxLengthX + maxLengthY) * (maxLengthX + maxLengthY));

		// minimize gx = distance to the ignored trace
		lp.setInt(transitions * (maxLengthX + 1) + 2, integerVariables);
		lp.setObjective(transitions * (maxLengthX + 1) + 2, -(maxLengthX + maxLengthY));

		lp.setMaxim();

		double[] rhs = new double[lp.getNrows()];

		// first -initial Marking
		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row < places * (maxLengthX + 1); row += places) {
				rhs[row] = -initialMarking.occurrences(p);
			}
		}
		// then add final marking
		for (Place p : finalMarking.baseSet()) {
			row = maxLengthX * places + place2int.get(p);
			rhs[row] += finalMarking.occurrences(p);
		}
		// then X.i - X.i+1 >= 0

		// then X.1 <= 1
		for (row = (maxLengthX + 1) * places + maxLengthX - 1; row < (maxLengthX + 1) * places + 2 * maxLengthX - 1; row++) {
			rhs[row] = 1;
		}
		// Y.1 <= maxLengthY
		rhs[row] = maxLengthY;
		row++;
		// Y.1 - l*X_l <= 0
		row++;

		for (int t = 0; t < log.length; t++) {
			if (startTracesAt > log[t].length) {
				// trace is outside window, length = 0;
				rhs[row + useHY * t] = 0;
			} else if (log[t].length > startTracesAt + maxLengthX) {
				// trace spans full window
				rhs[row + useHY * t] = -maxLengthX;
			} else {
				rhs[row + useHY * t] = -(log[t].length - startTracesAt);
			}
		}

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
	}

	protected HybridEquationResult determineSplitMarkingForHybrid(LPMatrix<?> matrix, int maxLengthX,
			boolean includeTrailingTaus) throws LPMatrixException {

		double[] vars = new double[matrix.getNcolumns()];
		int result = matrix.solve(vars);

		if (result == LPMatrix.OPTIMAL) {
			// compute the expected length of Y
			int lengthY = 0;
			for (int t = maxLengthX * transitions + invisibleTransitions; t < matrix.getNcolumns() - 3; t++) {
				if (vars[t] > 0.5) {
					lengthY += (int) (vars[t] + 0.5);
				}
			}

			TShortList antiAlignment = new TShortArrayList(maxLengthX);
			Stack<Transition> firingSequence = new Stack<Transition>();

			int lengthX = 0;
			// get the intermediate marking:
			Marking intermediate = new Marking();
			// Iterate over the X vector
			for (int p = 0; p < places; p++) {
				double effect = -matrix.getRh(p);
				for (int t = 0; t < maxLengthX * transitions + (includeTrailingTaus ? invisibleTransitions : 0); t++) {
					effect += vars[t] * matrix.getMat(p + places, t % transitions);
					if (p == 0 && vars[t] > 0.5) {
						if (trans2label[t % transitions] < 0) {
							double v = vars[t];
							do {
								firingSequence.add(short2trans[t % transitions]);
								v = v - 1;
							} while (v > 0.5);
						} else {
							lengthX++;
							// single visible transition
							assert ((int) (vars[t] + 0.5)) == 1;
							firingSequence.add(short2trans[t % transitions]);
							antiAlignment.add(label2short.get(short2trans[t % transitions].getLabel()));
						}
					}
				}
				assert (effect + 1E-10) >= 0;
				if (effect > 0) {
					intermediate.add(short2place[p], (int) (effect + 0.5));
				}
			}

			return new HybridEquationResult(intermediate, lengthX, lengthY, antiAlignment, firingSequence);
		} else {
			return null;
		}
	}
}
