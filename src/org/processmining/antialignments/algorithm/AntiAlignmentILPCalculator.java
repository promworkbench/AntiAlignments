package org.processmining.antialignments.algorithm;

import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.Vector;

import lpsolve.LpSolve;
import nl.tue.astar.util.LPMatrix;
import nl.tue.astar.util.LPMatrix.LPMatrixException;

import org.deckfour.xes.classification.XEventClass;
import org.processmining.antialignments.algorithm.ilp.HybridEquationResult;
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
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;

public class AntiAlignmentILPCalculator extends AbstractILPCalculator {

	private final int maxFactor;

	private final int maxLength;

	public AntiAlignmentILPCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<XEventClass> label2short, TShortObjectMap<XEventClass> short2label,
			TransEvClassMapping mapping, short[][] log, int maxLength, int maxFactor) {

		super(net, initialMarking, finalMarking, label2short, short2label, mapping, log);

		this.maxLength = maxLength;
		this.maxFactor = maxFactor;

	}

	protected void solveSequential(int maxLength, Marking initialMarking, Marking finalMarking,
			final Stack<Transition> firingSequence, final TShortList antiAlignment, final int traceToIgnore)
			throws LPMatrixException {

		int[] hammingDistances = new int[log.length];

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		LPMatrix<?> matrix;
		HybridEquationResult intermediate = new HybridEquationResult(initialMarking, 0, 0);
		int startTracesAt = 0;
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

			if (startTracesAt < 0) {
				// So we backtracked a bit to far :)
				System.err.println("HMM");
			}

			matrix = setupLpForHybrid(maxLengthX, maxLengthY, true, marking, finalMarking, traceToIgnore, startTracesAt);

			//			FileWriter writer;
			//			try {
			//				writer = new FileWriter("D:/temp/antialignment/debugLP.csv");
			//				matrix.printLp(writer, ";");
			//				writer.close();
			//				((LpSolve) matrix.toSolver()).writeLp("D:/temp/antialignment/debugLP.lp");
			//
			//			} catch (Exception e1) {
			//				return;
			//			}

			int[] hd = new int[log.length];
			// Determine an intermediate marking at lengthX visible steps
			HybridEquationResult nextResult = determineSplitMarkingForHybrid(matrix, maxLengthX, startTracesAt,
					maxLengthY == 0, antiAlignment, firingSequence, hd);

			if (VERBOSE && nextResult != null) {
				System.out.println("Marking reached: " + nextResult.getMarking() + " in " + nextResult.getLengthX()
						+ " steps and estimating " + nextResult.getLengthY() + " remaining steps.");
			}
			if (nextResult == null
					|| nextResult.getLengthX() + nextResult.getLengthY() < intermediate.getLengthY()
							/ backtrackThreshold) {
				// nextResult did not live up to its expectations. 
				// go back one step.
				if (nextResult == null) {
					if (VERBOSE) {
						System.out.println("Backtracking because of infeasibility.");
					}
				} else {
					if (VERBOSE) {
						System.out.println("Backtracking since " + (nextResult.getLengthX() + nextResult.getLengthY())
								+ " < " + intermediate.getLengthY());
					}
					nextResult.undo(antiAlignment, firingSequence);
				}

				if (antiAlignment.isEmpty()) {
					System.err.println("Infeasible model: " + matrix.getLastSolverResult());
					FileWriter writer;
					try {
						writer = new FileWriter("D:/temp/antialignment/debugLP-AntiAlignment.csv");
						matrix.printLp(writer, ";");
						writer.close();
						((LpSolve) matrix.toSolver()).writeLp("D:/temp/antialignment/debugLP-AntiAlignment.lp");

					} catch (Exception e1) {
						e1.printStackTrace();
					}

				}

				// remove the last transition from the stack.
				Transition last;
				short label = antiAlignment.get(antiAlignment.size() - 1);
				antiAlignment.removeAt(antiAlignment.size() - 1);
				// correct also the hamming distances...
				int evt = maxLength - maxLengthX - maxLengthY - 1;
				for (int t = 0; t < log.length; t++) {
					if (evt >= log[t].length || log[t][evt] != label) {
						hammingDistances[t]--;
					}
				}

				do {
					last = firingSequence.pop();

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
				} while (!firingSequence.isEmpty() && firingSequence.peek().isInvisible());

				// We've removed the last transition in the row, 
				// to avoid loops in backtracking, require one more step in X
				maxLengthX++;
				// and include one additional event.
				startTracesAt--;
			} else {
				for (int t = 0; t < log.length; t++) {
					hammingDistances[t] += hd[t];
				}
				int[] hdComputed = new int[log.length];
				for (int t = 0; t < log.length; t++) {
					hdComputed[t] = getHammingDistanceToTrace(log[t], antiAlignment.toArray())
							- (log[t].length > antiAlignment.size() ? log[t].length - antiAlignment.size() : 0);
				}
				assert Arrays.equals(hdComputed, hammingDistances);

				// Setup for next iteration
				intermediate = nextResult;
				splits.add(intermediate);
				maxLengthX = cutOffLength;
				marking = intermediate.getMarking();
				maxLengthY -= maxLengthX;
				startTracesAt += maxLengthX;
			}
		} while (maxLengthY >= 0 && !intermediate.getMarking().equals(finalMarking));

		if (!intermediate.getMarking().equals(finalMarking)) {
			if (VERBOSE) {
				System.out.println("Trying to get from " + marking + " to " + finalMarking + " in "
						+ (maxLengthX + maxLengthY) + " exact steps.");
			}
			matrix = setupLpForHybrid(maxLengthX + maxLengthY, 0, true, marking, finalMarking, traceToIgnore,
					startTracesAt);
			intermediate = determineSplitMarkingForHybrid(matrix, maxLengthX + maxLengthY, startTracesAt, true,
					antiAlignment, firingSequence, hammingDistances);
			splits.add(intermediate);
			if (VERBOSE) {
				System.out.println("Final marking reached.");
			}
		}

		assert intermediate.getMarking().equals(finalMarking);
		for (int t = 0; t < log.length; t++) {
			hammingDistances[t] += log[t].length > antiAlignment.size() ? log[t].length - antiAlignment.size() : 0;
		}

		if (VERBOSE) {
			int[] hdComputed = new int[log.length];
			for (int t = 0; t < log.length; t++) {
				hdComputed[t] = getHammingDistanceToTrace(log[t], antiAlignment.toArray());
			}
			assert Arrays.equals(hdComputed, hammingDistances);

			System.out.println("Anti alignment: " + toString(firingSequence));
			System.out.println("Hamming Distances solved:   " + Arrays.toString(hammingDistances));
			System.out.println("Hamming Distances computed: " + Arrays.toString(hdComputed));
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
		VERBOSE = false;
		mode = MODE_LPSOLVE;
		mode = MODE_GUROBI;

		String sep = "\t";
		double exp = 1;

		cutOffLength = 1;

		backtrackThreshold = 2.0;

		int maxCutOffLength = 40;// maxLength * maxFactor;

		System.out.println("mode" + sep + "cutOffLength" + sep + "time(ms)" + sep + "Dlog");
		while (cutOffLength < maxCutOffLength) {
			cutOffLength += 1;
			long sum = 0;
			long dist = 0;

			Stack<Transition> firingSequence;
			TShortList antiAlignment = null;

			for (int i = 0; i < exp; i++) {
				firingSequence = new Stack<Transition>();
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

		//		if (doAntiAlignmentTest(initialMarking, finalMarking)) {
		//			return null;
		//		}

		if (VERBOSE) {
			System.out.println("Solving using sequential hybrid approach.");
		}

		cutOffLength = 15;
		Stack<Transition> firingSequence = new Stack<Transition>();
		TShortList antiAlignment = new TShortArrayList(maxLength * maxFactor);
		//		solveByDrillingDown(maxLength * maxFactor, initialMarking, finalMarking, firingSequence, antiAlignment, -1, 0);
		solveSequential(maxLength * maxFactor, initialMarking, finalMarking, firingSequence, antiAlignment, -1);
		// FiringSequence and anti-alignment are known, but the distances to the log need to be computed.
		short[] aa = antiAlignment.toArray();
		int hd = getMinimalHammingDistanceToLog(aa, log, -1);
		cutOffLength = 15;

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
			System.out.println("Anti alignment: " + toString(firingSequence));
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

			firingSequence = new Stack<Transition>();
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
		return hd == Integer.MAX_VALUE ? antiAlignment.length : hd;
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
					t = trans2short.get(trans);
					dir = -((Arc) e).getWeight();
				} else {
					trans = (Transition) e.getSource();
					t = trans2short.get(trans);
					p = place2int.get(e.getTarget());
					dir = ((Arc) e).getWeight();
				}
			} else if (e instanceof InhibitorArc) {
				p = place2int.get(e.getSource());
				trans = (Transition) e.getTarget();
				t = trans2short.get(trans);
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
					if (trans2label[tr] < 0) {
						//tr is a tau step. They don't add to the hamming distance
						// this adds 1 to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 0);
					} else if (e >= log[t].length) {
						// e is passed the length of the trace
						// this adds 1 to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 1); // -1
					} else if (e < log[t].length && !equalLabel(tr, log[t][e])) {
						// e is in the trace window, but label does not match
						// this adds 1 to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 1); // 0
					} else {
						//e is in the trace window and label matches.
						// this does not add to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 0); //0
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

		row = (maxLength + 1) * places + 2 * maxLength + log.length - 1;
		rhs[row] = maxLength;

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
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

	private String toString(Vector<Transition> firingSequence) {
		if (firingSequence == null)
			return "null";
		int iMax = firingSequence.size() - 1;
		if (iMax == -1)
			return "[]";

		StringBuilder b = new StringBuilder();
		b.append('[');
		for (int i = 0; i < firingSequence.size(); i++) {
			if (!firingSequence.get(i).isInvisible()) {
				XEventClass clazz = mapping.get(firingSequence.get(i));
				if (clazz != null && !clazz.equals(mapping.getDummyEventClass())) {
					b.append(clazz);
					if (i == iMax)
						return b.append(']').toString();
					b.append(", ");
				}
			}
		}
		return b.toString();
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
					t = trans2short.get(trans);
					dir = -((Arc) e).getWeight();
				} else {
					trans = (Transition) e.getSource();
					t = trans2short.get(trans);
					p = place2int.get(e.getTarget());
					dir = ((Arc) e).getWeight();
				}
			} else if (e instanceof InhibitorArc) {
				p = place2int.get(e.getSource());
				trans = (Transition) e.getTarget();
				t = trans2short.get(trans);
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
					t = trans2short.get(trans);
					dir = -((Arc) e).getWeight();
				} else {
					trans = (Transition) e.getSource();
					t = trans2short.get(trans);
					p = place2int.get(e.getTarget());
					dir = ((Arc) e).getWeight();
				}
			} else if (e instanceof InhibitorArc) {
				p = place2int.get(e.getSource());
				trans = (Transition) e.getTarget();
				t = trans2short.get(trans);
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

		//-                            B    >= m0
		//-                            AB   >= m0
		//-                            AAB  >= m0
		//-                            AAAA == mf             Xi - Xi+1<=0   Xi <= 1   Y.1 <= maxLength-lengthX 
		LPMatrix<?> lp = setupMatrix(places * (maxLengthX + 1) + (maxLengthX - 1) + maxLengthX + 1/*                     */
				//- y.1 - maxLength X_(lengthX-1).1 <=0 traceXHD  rem.traceHD   
				+ 1 /*                             */+ log.length, transitions * maxLengthX + transitions + 2);

		for (int block = 0; block < maxLengthX + 1; block++) {
			int row = block * places;
			for (int j = 0; j < block; j++) {
				// Setup the A matrices
				matrixA.copyIntoMatrix(lp, row, j * transitions);
			}

			// Setup the A- matrices
			if (block < maxLengthX) {
				matrixAMin.copyIntoMatrix(lp, row, block * transitions);
			} else {
				matrixA.copyIntoMatrix(lp, row, block * transitions);
			}
			for (int p = 0; p < places; p++) {
				lp.setConstrType(block * places + p, LPMatrix.GE);
				if (NAMES) {
					lp.setRowName(block * places + p, "r" + block + "P" + p);
				}
			}

			for (int t = 0; t < transitions; t++) {
				if (NAMES) {
					lp.setColName(block * transitions + t, "c" + block + "T" + t);
				}
				if (block < maxLengthX) {
					lp.setBinary(block * transitions + t, integerVariables);
					lp.setUpbo(block * transitions + t, 1);
				} else {
					lp.setInt(block * transitions + t, integerVariables);
				}

				// OBJECTIVE IS SET LATER
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
			if (NAMES) {
				lp.setRowName(row + i, "X" + i + "~X" + (i + 1));
			}
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// sos constraints, sum max 1
		row += maxLengthX - 1;
		for (int i = 0; i < maxLengthX; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			if (NAMES) {
				lp.setRowName(row + i, "X" + i + ".1");
			}
			lp.setConstrType(row + i, LpSolve.LE);
		}

		row += maxLengthX;
		// setup Y length
		for (int t = 0; t < transitions; t++) {
			lp.setMat(row, maxLengthX * transitions + t, trans2label[t] < 0 ? 0 : 1);
		}
		if (NAMES) {
			lp.setRowName(row, "Y.1");
		}
		lp.setConstrType(row, LpSolve.LE);

		row++;
		// setup Y only if X_{lengthX}.1==1;
		for (int t = 0; t < transitions; t++) {
			lp.setMat(row, (maxLengthX - 1) * transitions + t, trans2label[t] < 0 ? 0 : -maxLengthY);
			lp.setMat(row, maxLengthX * transitions + t, trans2label[t] < 0 ? 0 : 1);
		}
		if (NAMES) {
			lp.setRowName(row, "Y~X_l");
		}
		lp.setConstrType(row, LpSolve.LE);

		double maxOccInY = 0;
		//Hamming distances for the X and Y vectors
		row++;
		for (short tr = 0; tr < trans2label.length; tr++) {
			for (int t = 0; t < log.length; t++) {
				// what's the hamming distance of the x's to trace t?
				for (int e = startTracesAt; e < startTracesAt + maxLengthX; e++) {
					if (trans2label[tr] < 0) {
						//tr is a tau step. They don't add to the hamming distance
						// this adds 0 to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 0);
					} else if (e >= log[t].length) {
						// e is passed the length of the trace
						// this adds 1 to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 1);
					} else if (e < log[t].length && !equalLabel(tr, log[t][e])) {
						// e is in the trace window, but label does not match
						// this adds 1 to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 1);
					} else {
						//e is in the trace window and label matches.
						// this does not add to the hamming distance
						lp.setMat(row + t, (e - startTracesAt) * transitions + tr, 0);
					}
				}

				if (tr == 0) {
					// Setup the rows for this trace;
					if (NAMES) {
						lp.setRowName(row + t, "Hx-t" + t);
					}

					if (t == traceToIgnore) {
						// minus HD for the removes trace
						lp.setMat(row + t, transitions * (maxLengthX + 1) + 1, -1);
						lp.setConstrType(row + t, LpSolve.EQ);
					} else {
						// minus HX
						lp.setMat(row + t, transitions * (maxLengthX + 1), -1);
						lp.setConstrType(row + t, LpSolve.GE);
					}

				}

			}
		}

		// Setup objective function

		// Minimize all invisible transitions
		for (int t = 0; t < transitions * maxLengthX; t++) {
			// Minimize the number of INVISIBLE transition firings in the anti-alignment
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
		double max = 0.0;
		//		 Minimize overlap in Y vector
		for (short tr = 0; tr < trans2label.length; tr++) {
			for (int t = 0; t < log.length; t++) {
				if (t == traceToIgnore) {
					continue;
				}
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

		if (NAMES) {
			lp.setColName(transitions * (maxLengthX + 1), "Dlog");
			lp.setColName(transitions * (maxLengthX + 1) + 1, "Drem");
		}
		// maximize Dlog
		if (log.length == 0 || (log.length == 1 && traceToIgnore >= 0)) {
			// a DLog is a 0-column. A positive weight would make the model unbounded.
			lp.setObjective(transitions * (maxLengthX + 1), -1);
			// maximize Drem in this case.
			lp.setObjective(transitions * (maxLengthX + 1) + 1, (maxLengthX + maxLengthY));
		} else {
			// maximize Dlog
			lp.setObjective(transitions * (maxLengthX + 1), (maxLengthX + maxLengthY) * (maxLengthX + maxLengthY));
			// minimize Drem
			lp.setObjective(transitions * (maxLengthX + 1) + 1, -(maxLengthX + maxLengthY));
		}

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

		// Y.1 - l*X_l <= 0

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
	}

	protected HybridEquationResult determineSplitMarkingForHybrid(LPMatrix<?> matrix, int maxLengthX,
			int startTracesAt, boolean includeTrailingTaus, TShortList antiAlignment, Stack<Transition> firingSequence,
			int[] hammingDistances) throws LPMatrixException {

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
								firingSequence.push(short2trans[t % transitions]);
								v = v - 1;
							} while (v > 0.5);
						} else {
							lengthX++;
							// single visible transition
							assert ((int) (vars[t] + 0.5)) == 1;
							firingSequence.push(short2trans[t % transitions]);
							antiAlignment.add(trans2label[t % transitions]);
						}
					}
				}
				assert (effect + 1E-10) >= 0;
				if (effect > 0) {
					intermediate.add(short2place[p], (int) (effect + 0.5));
				}
			}

			// compute the hamming distances for the prefixes
			int row = (maxLengthX + 1) * places + 2 * maxLengthX + 1;
			for (int t = 0; t < log.length; t++) {
				double hd_t = 0;// -matrix.getRh(row);
				for (int v = 0; v < maxLengthX * transitions; v++) {
					hd_t += vars[v] * matrix.getMat(row, v);
				}

				// round hamming distance
				hammingDistances[t] += (int) (hd_t + 0.5);
				row++;
			}

			return new HybridEquationResult(intermediate, lengthX, lengthY);
		} else {
			return null;
		}
	}
}
