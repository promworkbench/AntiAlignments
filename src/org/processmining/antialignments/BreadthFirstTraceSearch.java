package org.processmining.antialignments;

import gnu.trove.set.TShortSet;

import java.util.HashSet;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;

import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class BreadthFirstTraceSearch {

	private final Petrinet net;
	private final Map<Marking, TShortSet> statesVisitedPerTrace;
	private final PetrinetSemantics semantics;

	public BreadthFirstTraceSearch(Petrinet net, Marking initialMarking,
			Map<Marking, TShortSet> statesVisitedPerTrace) {
		this.net = net;
		this.statesVisitedPerTrace = statesVisitedPerTrace;

		this.semantics = PetrinetSemanticsFactory
				.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);
	}

	private class MarkingVisit {
		public final Marking marking;
		public final int distance;

		public MarkingVisit(Marking marking, int distance) {
			this.marking = marking;
			this.distance = distance;
		}
	}

	public int getShortestDistance(Marking currentMarking, int maxDistance,
			short traceToIgnore) {
		Set<Marking> done = new HashSet<Marking>();
		Queue<MarkingVisit> todo = new LinkedBlockingQueue<MarkingVisit>();

		todo.add(new MarkingVisit(currentMarking, 0));
		while (!todo.isEmpty()) {
			MarkingVisit m = todo.poll();
			if (done.contains(m.marking)) {
				continue;
			}

			TShortSet traces = statesVisitedPerTrace.get(m.marking);
			if (traces != null
					&& (traces.size() > 1 || !traces.contains(traceToIgnore))) {
				return m.distance;
			}

			done.add(m.marking);
			semantics.setCurrentState(m.marking);
			for (Transition t : semantics.getExecutableTransitions()) {
				try {
					semantics.executeExecutableTransition(t);
				} catch (IllegalTransitionException e) {
					e.printStackTrace();
				}
				todo.add(new MarkingVisit(semantics.getCurrentState(),
						m.distance + 1));
				semantics.setCurrentState(m.marking);
			}
		}
		return Integer.MAX_VALUE;
	}
}
