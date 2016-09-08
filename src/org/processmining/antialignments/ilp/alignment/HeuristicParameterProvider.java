package org.processmining.antialignments.ilp.alignment;

import javax.swing.JComponent;

import org.deckfour.xes.model.XLog;
import org.processmining.framework.plugin.PluginContext;
import org.processmining.models.graphbased.directed.petrinet.PetrinetGraph;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;
import org.processmining.plugins.petrinet.replayer.algorithms.IPNReplayParameter;
import org.processmining.plugins.petrinet.replayer.algorithms.costbasedcomplete.CostBasedCompleteParamProvider;

public class HeuristicParameterProvider extends CostBasedCompleteParamProvider {

	public HeuristicParameterProvider(PluginContext context, PetrinetGraph net, XLog log, TransEvClassMapping mapping) {
		super(context, net, log, mapping);

	}

	public IPNReplayParameter constructReplayParameter(JComponent ui) {
		if (ui instanceof HeuristicParameterUI) {
			HeuristicParameterUI hui = (HeuristicParameterUI) ui;

			HeuristicParameters paramObj = new HeuristicParameters(hui.getMapEvClassToCost(), hui.getTransitionWeight());
			paramObj.setMapSync2Cost(hui.getSyncCost());
			paramObj.setMaxNumOfStates(hui.getMaxNumOfStates());
			paramObj.setInitialMarking(initMarking);
			paramObj.setFinalMarkings(finalMarkings);
			//			paramObj.setUsePartialOrderedEvents(hui.isUsePartialOrderedEvents());

			paramObj.setExpectedModelMoves(hui.getExpectedModelMoves());
			paramObj.setBacktrackLimit(hui.getBacktrackLimit());
			paramObj.setBacktrackThreshold(hui.getBacktrackThreshold());

			return paramObj;
		} else {
			return null;
		}
	}

	public JComponent constructUI() {
		return new HeuristicParameterUI(transCol, evClassCol, transCol.size());
	}
}
