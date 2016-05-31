package org.processmining.antialignments.pathfinder;

import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Vector;

import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class AntiAlignmentFinder {

	private final Petrinet net;
	private final PetrinetSemantics semantics;
	private final Marking initialMarking;
	private final TObjectShortMap<String> label2short;
	private final Marking finalMarking;

	public static int SCALING = 1000;
	private static final boolean DOT = false;

	public AntiAlignmentFinder(Petrinet net, Marking initialMarking,
			Marking finalMarking,
			// Map<Marking, TShortSet> statesVisitedPerTrace,
			TObjectShortMap<String> label2short) {
		this.net = net;
		this.initialMarking = initialMarking;
		this.finalMarking = finalMarking;
		this.label2short = label2short;

		this.semantics = PetrinetSemanticsFactory
				.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);
	}

	protected Map<State, State> toDoSet = new HashMap<State, State>();
	protected TObjectIntMap<State> stateMap = new TObjectIntHashMap<State>();

	protected PriorityQueue<State> toDo = new PriorityQueue<>(
			new Comparator<State>() {

				@Override
				public int compare(State o1, State o2) {
					// States should be compared on their delta-distance
					double d1 = o2.getPathLength();
					double d2 = o1.getPathLength();

					// int d1 = o1.getMinDistance();
					// int d2 = o2.getMinDistance();

					// first order sorting based on distance
					//
					if (d1 == d2) {
						d1 = o1.getAntiAlignmentLength();
						d2 = o2.getAntiAlignmentLength();
					}

					return d2 < d1 ? -1 : d2 > d1 ? 1 : 0;
				}

			});

	public AntiAlignments getAntiAlignment(final short[][] log, int maxLength,
			final double maxFactor) {

		AntiAlignments aa = new AntiAlignments(log.length);

		int[][] current = new int[log.length][];

		int j = 0;
		int min = Integer.MAX_VALUE;
		for (short[] trace : log) {
			int i = 0;
			current[j] = new int[trace.length + 1];
			for (short a : trace) {
				current[j][++i] = i;
			}
			if (trace.length < min) {
				min = trace.length;
			}
			j++;
		}

		maxLength = (int) (maxLength * maxFactor + 0.5);

		final State initialState = new State(initialMarking, null, null,
				current, 0, 0);// Math.abs(min - maxLength));
		// System.out.println("InitialState: " + initialState.getMinDistance());

		for (int t = log.length; t >= 0; t--) {
			if (DOT) {
				System.out.println("Digraph D {");
			}
			toDo.clear();
			toDoSet.clear();

			toDo.add(initialState);
			toDoSet.put(initialState, initialState);

			while (!toDo.isEmpty()) {
				State s = toDo.poll();
				toDoSet.remove(s);

				assert toDo.isEmpty()
						|| toDo.peek().getPathLength() >= s.getPathLength();

				if (s.getAntiAlignmentLength() > maxLength) {
					// No final marking reached with given length...
					// stop search from this state and continue the search;
					continue;
				}

				if (expand(log, s, t, maxLength)) {
					// Done!
					// s is the final state and we found an anti-alignment.
					// if (s.getMinDistance() > aa.getMaxDistances()[t]) {

					// store the minimal distance of the anti-alignment found to
					// the Log
					// aa.getMaxDistances()[t] = 1 - s.getPathLength();

					// TODO!
					aa.getMaxMinDistances()[t] = s.getMinMinDistance(t);

					// store the distance of the anti-alignment found to the
					// removed trace
					if (t < log.length) {
						aa.getTraceDistances()[t] = s.getDistance(t);
					}

					Vector<Transition> firingSequence = new Vector<Transition>(
							s.getAntiAlignmentLength());
					aa.getTraces()[t] = firingSequence;
					aa.getAntiAlignments()[t] = new short[s
							.getAntiAlignmentLength()];
					// int[][] matrix = new int[s.getAntiAlignmentLength() +
					// 1][];
					int i = s.getAntiAlignmentLength();
					do {
						Transition transition = s.getFiredTransition();
						firingSequence.add(0, transition);
						if (!transition.isInvisible()) {
							// matrix[i] = s.getCurrent()[t];
							aa.getAntiAlignments()[t][--i] = label2short
									.get(transition.getLabel());
						}
						s = s.getPredecessor();
					} while (s.getFiredTransition() != null);
					// matrix[0] = s.getCurrent()[t];
					//
					// System.out.println("    " + Arrays.toString(log[t]));
					// for (i = 0; i < matrix.length; i++) {
					// if (i > 0) {
					// System.out.print(aa.getAntiAlignments()[t][i - 1]);
					// } else {
					// System.out.print(" ");
					// }
					// System.out.println(Arrays.toString(matrix[i]));
					// }
					// }
					toDo.clear();
					toDoSet.clear();
					// System.exit(0);
				}
			}
		}

		return aa;
	}

	// private short[][] forced = new short[][] { { 8, 6, 4, 2, 0, 1, 5, 3, 7 },
	// { 8, 6, 4, 2, 1, 7, 0, 3, 5 }, { 8, 6, 4, 2, 1, 7, 0, 3, 5 },
	// { 8, 6, 4, 2, 1, 7, 0, 3, 5 }, { 8, 6, 4, 2, 1, 7, 0, 3, 5 },
	// { 8, 6, 4, 2, 1, 7, 0, 3, 5 } };

	// RETURN TRUE if state s is the final state!
	protected boolean expand(final short[][] log, State s, int traceToIgnore,
			int maxLength) {

		Marking currentMarking = s.getMarking();
		semantics.setCurrentState(currentMarking);
		Collection<Transition> enabled = semantics.getExecutableTransitions();

		if (currentMarking.equals(finalMarking)
				&& (enabled.isEmpty() || s.getAntiAlignmentLength() == maxLength)) {
			// final marking reached, nothing to do in the model, but there may
			// be better anti-alignments in the toDo list that can still expand
			// to a greater distance.
			//
			// Reschedule this state with definitive potential.
			// if (s.getPathLength() == s.getMinDistance()) {
			return true;
			// } else {
			// s.setPotential(s.getMinDistance());
			// toDo.add(s);
			// toDoSet.put(s, s.getPotential());
			// return false;
			// }
		}

		for (Transition t : enabled) {
			try {
				semantics.setCurrentState(currentMarking);
				semantics.executeExecutableTransition(t);
			} catch (IllegalTransitionException e) {
				e.printStackTrace();
			}

			double len = 0;
			State newState;
			// boolean explain = false;
			if (t.isInvisible()) {
				// ignore invisible transitions in the distance metric
				// invisible transitions are considered non-edits.
				newState = new State(semantics.getCurrentState(), s, t,
						s.getCurrent(), s.getFiringSequenceLength() + 1,
						s.getPathLength());
			} else {
				short label = label2short.get(t.getLabel());

				// if (s.getAntiAlignmentLength() < forced[traceToIgnore].length
				// && label == forced[traceToIgnore][s
				// .getAntiAlignmentLength()]) {
				// State temp = s;
				// int i;
				// explain = true;
				// if (s.getAntiAlignmentLength() > 0) {
				// do {
				// i = temp.getAntiAlignmentLength();
				// explain &= forced[traceToIgnore][i - 1] == label2short
				// .get(temp.getFiredTransition().getLabel());
				// temp = temp.getPredecessor();
				// } while (explain && temp.getAntiAlignmentLength() > 0);
				// }
				// }

				int[][] previous = s.getCurrent();
				int[][] current = new int[log.length][];

				double traceMin = Double.MAX_VALUE;
				int partMin = Integer.MAX_VALUE;

				for (int tr = 0; tr < log.length; tr++) {
					short[] trace = log[tr];

					current[tr] = new int[trace.length + 1];
					current[tr][0] = previous[tr][0] + 1;
					if (tr != traceToIgnore && current[tr][0] < partMin) {
						partMin = current[tr][0];
					}

					for (int i = 1; i < trace.length + 1; i++) {
						if (trace[i - 1] == label) {
							current[tr][i] = previous[tr][i - 1];
						} else {
							int min = Integer.MAX_VALUE;
							if (previous[tr][i - 1] < min) {
								// replace
								min = previous[tr][i - 1];
							}
							if (previous[tr][i] < min) {
								// insert
								min = previous[tr][i];
							}
							if (current[tr][i - 1] < min) {
								// delete
								min = current[tr][i - 1];
							}
							current[tr][i] = min + 1;
						}
						if (tr != traceToIgnore && current[tr][i] < partMin) {
							partMin = current[tr][i];
						}

					}
					if (tr != traceToIgnore) {

						if (current[tr][trace.length]
								/ (double) Math.max(trace.length,
										current[tr][0]) < traceMin) {
							traceMin = current[tr][trace.length]
									/ (double) Math.max(trace.length,
											current[tr][0]);
						}

						// if (current[tr][0] < current[tr].length) {
						// if (current[tr][current[tr][0]] < partMin) {
						// partMin = current[tr][current[tr][0]];
						// }
						// } else {
						// if (traceMin < partMin) {
						// partMin = traceMin;
						// }
						// }
					}

				}

				// The distance of this edge is defined by the loss of
				// potential.

				newState = new State(semantics.getCurrentState(), s, t,
				// current, s.getFiringSequenceLength() + 1, 1 - traceMin);
						current, s.getFiringSequenceLength() + 1, -partMin);

			}
			State existing = toDoSet.get(newState);

			if (existing == null
					|| existing.getPathLength() > +newState.getPathLength()) {

				// state is either new, or is in the queue with greater delta
				// indicating that it was reached before trough a longer path.
				toDo.remove(newState);
				toDoSet.remove(newState);
				toDo.add(newState);
				toDoSet.put(newState, newState);

				if (DOT && existing == null) {
					System.out.println("S" + stateMap.get(s) + " -> S"
							+ stateMap.get(newState)
							+ " [color=\"green\",label=\"" + t + "," + len
							+ "\"];");
				} else if (DOT) {
					System.out.println("S" + stateMap.get(s) + " -> S"
							+ stateMap.get(newState) + " [label=\"" + t + ","
							+ len + "\"];");
				}
			}
			// TODO: IF A STATE IS REACHED AGAIN WITH NEW, EQUAL PATH LENGHT,
			// HOW TO HANDLE THE CASE WHERE THE NEW PATH IS THE ONE NEEDED TO
			// GET
			// BETTER EDIT DISTANCE IN THE FUTURE?????

			// else if (explain && existing != null
			// && existing.getPathLength() == newState.getPathLength()) {
			// System.out.println("State " + newState
			// + " was not scheduled when ignoring trace "
			// + traceToIgnore);
			// System.out.println("Existing  "
			// + Arrays.deepToString(existing.getCurrent())
			// + " (min: " + existing.getMinDistance());
			// System.out.println("New state "
			// + Arrays.deepToString(newState.getCurrent())
			// + " (min: " + existing.getMinDistance());
			// }
		}

		return false;

	}
}
