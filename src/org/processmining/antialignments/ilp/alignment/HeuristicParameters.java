package org.processmining.antialignments.ilp.alignment;

import java.util.Collection;
import java.util.Map;

import org.deckfour.xes.classification.XEventClass;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.plugins.petrinet.replayer.algorithms.costbasedcomplete.CostBasedCompleteParam;

public class HeuristicParameters extends CostBasedCompleteParam {

	private int expectedModelMoves;

	public HeuristicParameters(Collection<XEventClass> evClassCol, XEventClass dummyEvClass,
			Collection<Transition> transCol) {
		super(evClassCol, dummyEvClass, transCol);
	}

	public HeuristicParameters(Collection<XEventClass> evClassCol, XEventClass dummyEvClass,
			Collection<Transition> transCol, int defMoveOnModelCost, int defMoveOnLogCost) {
		super(evClassCol, dummyEvClass, transCol, defMoveOnModelCost, defMoveOnLogCost);
	}

	public HeuristicParameters(Map<XEventClass, Integer> mapEvClass2Cost, Map<Transition, Integer> mapTrans2Cost) {
		super(mapEvClass2Cost, mapTrans2Cost);
	}

	public HeuristicParameters(Map<XEventClass, Integer> mapEvClass2Cost, Map<Transition, Integer> mapTrans2Cost,
			Map<Transition, Integer> mapSync2Cost) {
		super(mapEvClass2Cost, mapTrans2Cost, mapSync2Cost);
	}

	public void setExpectedModelMoves(int expectedModelMoves) {
		this.expectedModelMoves = expectedModelMoves;

	}

	public int getExpecteModelMoves() {
		return this.expectedModelMoves;
	}
}
