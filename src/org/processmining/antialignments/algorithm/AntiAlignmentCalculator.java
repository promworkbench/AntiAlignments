package org.processmining.antialignments.algorithm;

import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.hash.TObjectShortHashMap;

import java.util.ArrayList;
import java.util.List;

import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class AntiAlignmentCalculator {

	protected final static boolean VERBOSE = false;

	private final Petrinet net;
	private final Marking initialMarking;
	private final PetrinetSemantics semantics;
	private final Marking finalMarking;
	private final TObjectShortMap<String> label2short;

	public AntiAlignmentCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short) {
		this.net = net;
		this.initialMarking = initialMarking;
		this.finalMarking = finalMarking;
		this.label2short = label2short;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

	}

	public AntiAlignments getAntiAlignments(final short[][] log, int maxLength, int maxFactor) {
		int max = 0;
		for (int t = 0; t < log.length; t++) {
			if (log[t].length > max) {
				max = log[t].length;
			}
		}
		//		State initialState = new HammingState(initialMarking, log.length);
		State initialState = new EditDistanceState(initialMarking, log);

		short i = 0;
		TObjectShortMap<Transition> trans2short = new TObjectShortHashMap<>(net.getTransitions().size() * 2, 0.7f,
				(short) -1);
		for (Transition t : net.getTransitions()) {
			trans2short.put(t, i++);
		}
		TObjectShortMap<Place> place2short = new TObjectShortHashMap<>(net.getPlaces().size() * 2, 0.7f, (short) -1);
		i = 0;
		for (Place p : net.getPlaces()) {
			place2short.put(p, i++);
		}

		final AntiAlignments antiAlignments = new AntiAlignments(log.length);

		for (int traceToIgnore = log.length; traceToIgnore >= 0; traceToIgnore--) {
			int length = traceToIgnore < log.length ? log[traceToIgnore].length * maxFactor : maxLength * maxFactor;

			// Initialize the initial generation (note that initialMarking may be the final Marking
			SearchQueue<State> current = new SearchQueue<State>();
			current.add(initialState);

			if (VERBOSE) {
				System.out.println("Looking for anti alignment of maximum length: " + length);
			}
			State finalState = null;
			do {
				finalState = update(current, log, traceToIgnore, length);
				//System.out.println(current.toString());
			} while (finalState == null && !current.isEmpty());

			antiAlignments.getMaxMinDistances()[traceToIgnore] = finalState.getMinimalDistance(log, traceToIgnore);
			antiAlignments.getAntiAlignments()[traceToIgnore] = finalState.getAntiAlignment();
			antiAlignments.getTraces()[traceToIgnore] = finalState.getFiringSequence();

			if (VERBOSE) {
				System.out.print("AA:   " + finalState.getAntiAlignmentString());
				System.out.print("  Dlog: " + antiAlignments.getMaxMinDistances()[traceToIgnore]);
			}
			if (traceToIgnore < log.length) {
				antiAlignments.getTraceDistances()[traceToIgnore] = finalState.getDistance(log, traceToIgnore);
				if (VERBOSE) {
					System.out.print("  Drem: " + finalState.getDistance(log, traceToIgnore));
				}
			}
			if (VERBOSE) {
				System.out.println();
				finalState.printMatrix(System.out);
			}

		}

		return antiAlignments;
	}

	protected <S extends State> S update(SearchQueue<S> current, final short[][] log, final int traceToIgnore,
			int maxLength) {

		S s = current.pull();

		// what if s is a final state?
		if (s.getMarking().equals(finalMarking)) {
			return s;
		}

		// Iterate over the states
		semantics.setCurrentState(s.getMarking());
		List<Transition> transitions = new ArrayList<Transition>(semantics.getExecutableTransitions());
		//Collections.shuffle(transitions, new Random(1134246));

		// for each executable transition
		for (Transition t : transitions) {

			// Execute the transition
			semantics.setCurrentState(s.getMarking());
			try {
				semantics.executeExecutableTransition(t);
			} catch (IllegalTransitionException e) {
				e.printStackTrace();
			}

			// get the new Marking
			Marking newMarking = semantics.getCurrentState();

			short label = State.NOLABEL;
			if (!t.isInvisible()) {
				label = label2short.get(t.getLabel());
			}

			// determine new state
			S newState = s.getNextState(newMarking, log, traceToIgnore, label, t);

			if (newState.getLength() > maxLength) {
				// ignore too long sequences.
				continue;
			}
			if (newMarking.equals(finalMarking)) {
				// update the state
				newState.setFinalMarkingReached(log, traceToIgnore);
			}

			// add to next generation.
			if (VERBOSE) {
				System.out.println(s + " --> " + newState);
			}
			current.add(newState);

		}
		return null;
	}
}
