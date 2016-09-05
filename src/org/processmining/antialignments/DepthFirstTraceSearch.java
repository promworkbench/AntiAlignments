package org.processmining.antialignments;

import gnu.trove.map.TObjectShortMap;
import gnu.trove.procedure.TObjectProcedure;
import gnu.trove.set.hash.TCustomHashSet;
import gnu.trove.strategy.HashingStrategy;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Stack;
import java.util.Vector;

import org.processmining.antialignments.TestAntiAlignment.MarkedNet;
import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class DepthFirstTraceSearch {

	private static interface TraceHandler {
		/**
		 * Handle the observed trace and return true if the search should
		 * continue, false if the search can be aborted.
		 * 
		 * @param firingSequence
		 * @param word
		 * @return
		 */
		public boolean handleTrace(Vector<Transition> firingSequence, short[] word);
	}

	private final PetrinetSemantics semantics;
	private final Petrinet net;
	private final Marking initialMarking;

	private final Stack<Marking> currentStates;
	private final Stack<Transition> currentTransitions;
	private short[] currentLabels;
	private final TObjectShortMap<String> label2short;
	private final Marking finalMarking;

	public DepthFirstTraceSearch(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short) {
		this.net = net;
		this.initialMarking = initialMarking;
		this.finalMarking = finalMarking;
		this.label2short = label2short;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

		this.currentStates = new Stack<>();
		currentStates.push(initialMarking);

		this.currentTransitions = new Stack<>();
	}

	public DepthFirstTraceSearch(MarkedNet markedNet, TObjectShortMap<String> label2short) {
		this(markedNet.net, markedNet.initialMarking, markedNet.finalMarking, label2short);
	}

	/**
	 * Produces the language of the model in terms of sequences of observable
	 * transitions. Note that this procedure may return a very large number of
	 * sequences depending on the model's language. Use with care.
	 * 
	 * Be careful. This method does not terminate if there are loops of
	 * invisible transitions in the model!
	 * 
	 * @param maxLength
	 * @return
	 */
	public short[][] getLanguage(int maxLength) {
		// build all traces of length "max" for the model.

		final TCustomHashSet<short[]> language = new TCustomHashSet<short[]>(new HashingStrategy<short[]>() {

			private static final long serialVersionUID = -5834612520934860166L;

			@Override
			public int computeHashCode(short[] arg0) {
				return Arrays.hashCode(arg0);
			}

			@Override
			public boolean equals(short[] arg0, short[] arg1) {
				return Arrays.equals(arg0, arg1);
			}
		}) {

			public String toString() {
				final StringBuilder buf = new StringBuilder("{");
				forEach(new TObjectProcedure<short[]>() {
					private boolean first = true;

					public boolean execute(short[] value) {
						if (first) {
							first = false;
						} else {
							buf.append(", ");
						}

						buf.append(Arrays.toString(value));
						return true;
					}
				});
				buf.append("}");
				return buf.toString();
			}

		};

		depthFirstSearch(maxLength, new TraceHandler() {

			@Override
			public boolean handleTrace(Vector<Transition> firingSequence, short[] word) {
				language.add(word);
				return true;
			}
		});

		return language.toArray(new short[language.size()][]);
	}

	/**
	 * Search through the trace-space of the model and find a trace for which
	 * the minimal distance to any of the traces in the log is as high as
	 * possible. This procedure considers the whole log as well as each sublog
	 * where one trace was removed.
	 * 
	 * The search continues exhaustively until all traces of maximal length
	 * maxLength resulting in a deadlock have been processed, or until, for each
	 * sublog/log, and anti alignment is found with maximal minimal distance,
	 * i.e. where minimal distance is equal to the length of the shortest trace
	 * in the sublog.
	 * 
	 * Be careful. This method does not terminate if there are loops of
	 * invisible transitions in the model!
	 * 
	 * @param log
	 * @param maxLength
	 * @param metric
	 */
	public AntiAlignments getAntiAlignments(final short[][] log, final int maxLength, final double maxFactor,
			final DistanceMetric metric) {

		final AntiAlignments antiAlignments = new AntiAlignments(log.length, maxLength, maxFactor);

		Arrays.fill(antiAlignments.getMaxDistances(), Integer.MAX_VALUE);
		for (int t = 0; t < log.length; t++) {
			for (int u = log.length + 1; u-- > 0;) {
				if (u != t && antiAlignments.getMaxDistances()[u] > (int) (log[t].length * maxFactor)) {
					antiAlignments.getMaxDistances()[u] = (int) (log[t].length * maxFactor);
				}
			}
		}
		if (metric instanceof DistanceMetric.Hamming) {
			// For hamming distances, the distance is always bound by the
			// minimum of the shortest trace in the log and the removed trace
			for (int t = 0; t < log.length; t++) {
				if (antiAlignments.getMaxDistances()[t] > (int) (log[t].length * maxFactor)) {
					antiAlignments.getMaxDistances()[t] = (int) (log[t].length * maxFactor);
				}
			}
		} else if (metric instanceof DistanceMetric.Edit) {
			// For edit distances, the distance is always bound to maximum
			// of the removed trace and the minimum trace length in the log
			for (int t = 0; t < log.length; t++) {
				if (antiAlignments.getMaxDistances()[t] < (int) (log[t].length * maxFactor)) {
					antiAlignments.getMaxDistances()[t] = (int) (log[t].length * maxFactor);
				}
				if (antiAlignments.getMaxDistances()[log.length] < (int) (log[t].length * maxFactor)) {
					antiAlignments.getMaxDistances()[log.length] = (int) (log[t].length * maxFactor);
				}
			}
		}

		depthFirstSearch((int) (maxLength * maxFactor), new TraceHandler() {
			int calls = 0;
			final int length = log.length;
			final int[] currentDistances = new int[log.length];
			final int[] currentMin = new int[log.length + 1];

			@Override
			public boolean handleTrace(Vector<Transition> firingSequence, short[] word) {
				calls++;
				// determine if minimum higher than current minimum for:
				// 1) whole log,
				// 2) log with each trace removed
				// Store firing sequence and word if indeed better
				Arrays.fill(currentMin, Integer.MAX_VALUE);
				for (int t = length; t-- > 0;) {
					currentDistances[t] = metric.getDistance(log[t], word);
					for (int u = length + 1; u-- > 0;) {
						if (u != t && currentDistances[t] < currentMin[u]) {
							currentMin[u] = currentDistances[t];
						}
					}
				}
				// the array min contains the minimum distances for the log
				// with one trace removed and the minimal distance to the entire
				// log in the last position;
				boolean stop = true;
				for (int t = length + 1; t-- > 0;) {
					// ignore if word is longer than the max for the removed
					// trace!
					if ((t == length || word.length <= (int) (log[t].length * maxFactor))
							&& currentMin[t] > antiAlignments.getMaxMinDistances()[t]) {
						// new anti-alignment found!
						antiAlignments.getMaxMinDistances()[t] = currentMin[t];
						antiAlignments.getAntiAlignments()[t] = word;
						antiAlignments.getTraces()[t] = firingSequence;
						if (t < length) {
							antiAlignments.getTraceDistances()[t] = metric.getDistance(log[t], word);
						}
					}
					stop &= antiAlignments.getMaxMinDistances()[t] == antiAlignments.getMaxDistances()[t];
				}
				// if (stop) {
				// System.out.println("Processed " + calls + " traces.");
				// }
				return !stop;
			}
		});

		return antiAlignments;
	}

	protected void depthFirstSearch(final int maxLength, final TraceHandler handler) {
		this.currentLabels = new short[maxLength];
		depthFirstSearchInternal(maxLength, 0, handler);

	}

	protected boolean depthFirstSearchInternal(final int maxLength, int length, final TraceHandler handler) {

		semantics.setCurrentState(currentStates.peek());
		List<Transition> transitions = new ArrayList<Transition>(semantics.getExecutableTransitions());
		Collections.shuffle(transitions, new Random(1134246));
		if (currentStates.peek().equals(finalMarking)) {
			// complete trace found. Add copy of currentlabels to language.
			if (!handler.handleTrace((Vector<Transition>) currentTransitions.clone(),
					Arrays.copyOf(currentLabels, length))) {
				// abort the search!
				return false;
			}
		}

		for (Transition t : transitions) {
			semantics.setCurrentState(currentStates.peek());

			try {
				semantics.executeExecutableTransition(t);
			} catch (IllegalTransitionException e) {
				e.printStackTrace();
			}
			currentStates.push(semantics.getCurrentState());
			currentTransitions.push(t);
			boolean continueNext = true;
			if (t.isInvisible()) {
				// do not add the label to the word and continue with
				// same length
				continueNext = depthFirstSearchInternal(maxLength, length, handler);
			} else if (length < maxLength) {
				// add the label to the word and continue with
				// length +1
				currentLabels[length] = label2short.get(t.getLabel());
				continueNext = depthFirstSearchInternal(maxLength, length + 1, handler);
			}
			currentTransitions.pop();
			currentStates.pop();
			if (!continueNext) {
				return false;
			}
		}
		return true;
	}
}
