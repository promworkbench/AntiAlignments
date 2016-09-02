package org.processmining.antialignments.alignments;

import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gnu.trove.map.hash.TShortObjectHashMap;

import java.util.Map;
import java.util.Set;

import nl.tue.astar.AStarException;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.classification.XEventClasses;
import org.deckfour.xes.classification.XEventClassifier;
import org.deckfour.xes.info.XLogInfo;
import org.deckfour.xes.info.XLogInfoFactory;
import org.deckfour.xes.model.XEvent;
import org.deckfour.xes.model.XLog;
import org.deckfour.xes.model.XTrace;
import org.processmining.antialignments.algorithm.ilp.LPMatrix.LPMatrixException;
import org.processmining.framework.plugin.PluginContext;
import org.processmining.framework.plugin.annotations.KeepInProMCache;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetGraph;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;
import org.processmining.plugins.petrinet.replayer.algorithms.IPNReplayAlgorithm;
import org.processmining.plugins.petrinet.replayer.algorithms.IPNReplayParamProvider;
import org.processmining.plugins.petrinet.replayer.algorithms.IPNReplayParameter;
import org.processmining.plugins.petrinet.replayer.algorithms.costbasedcomplete.CostBasedCompleteParam;
import org.processmining.plugins.petrinet.replayer.algorithms.costbasedcomplete.CostBasedCompleteParamProvider;
import org.processmining.plugins.petrinet.replayer.annotations.PNReplayAlgorithm;
import org.processmining.plugins.petrinet.replayresult.PNRepResult;

@KeepInProMCache
@PNReplayAlgorithm
public class HeuristicPNetReplayerAlgorithm implements IPNReplayAlgorithm {

	private Map<Transition, Integer> mapTrans2Cost;
	private Integer maxNumOfStates;
	private Map<XEventClass, Integer> mapEvClass2Cost;
	private Map<Transition, Integer> mapSync2Cost;
	private Marking initialMarking;
	private Marking finalMarking;
	private boolean usePartialOrderEvents;
	private XEventClassifier classifier;
	private TShortObjectMap<String> short2label;
	private TObjectShortMap<String> label2short;

	public PNRepResult replayLog(PluginContext context, PetrinetGraph net, XLog xLog, TransEvClassMapping mapping,
			IPNReplayParameter parameters) throws AStarException {

		importParameters((CostBasedCompleteParam) parameters);
		classifier = mapping.getEventClassifier();

		TObjectIntMap<TShortList> tempLog = new TObjectIntHashMap<>(xLog.size());
		short2label = new TShortObjectHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);
		label2short = new TObjectShortHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);

		short c = 0;
		for (XTrace trace : xLog) {
			TShortList list = new TShortArrayList(trace.size());
			for (XEvent event : trace) {
				String clazz = classifier.getClassIdentity(event);
				short id = label2short.putIfAbsent(clazz, c);
				if (id == label2short.getNoEntryValue()) {
					short2label.put(c, clazz);
					id = c;
					c++;
				}
				list.add(id);
			}
			tempLog.adjustOrPutValue(list, 1, 1);
		}
		for (Transition t : net.getTransitions()) {
			XEventClass clazz = mapping.get(t);
			short id = label2short.putIfAbsent(clazz.getId(), c);
			if (id == label2short.getNoEntryValue()) {
				short2label.put(c, clazz.getId());
				id = c;
				c++;
			}
			label2short.putIfAbsent(t.getLabel(), id);
		}

		int t = 0;
		int[] frequencies = new int[tempLog.size()];
		short[][] log = new short[tempLog.size()][];
		TObjectIntIterator<TShortList> it = tempLog.iterator();
		while (it.hasNext()) {
			it.advance();
			log[t] = it.key().toArray();
			frequencies[t] = it.value();
			t++;
		}

		AlignmentILPCalculator calculator = new AlignmentILPCalculator(net, initialMarking, finalMarking, label2short,
				short2label, log);
		try {
			for (int tr = 1; tr < log.length; tr++) {
				calculator.solveSequential(initialMarking, finalMarking, tr);
			}
		} catch (LPMatrixException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return null;
	}

	public String getHTMLInfo() {
		return "<html>This is an algorithm to calculate cost-based fitness between a log and a Petri net. <br/><br/>"
				+ "Given a trace and a Petri net, this algorithm "
				+ "return a matching between the trace and an allowed firing sequence of the net with the"
				+ "least deviation cost using a heuristic technique. The firing sequence has to reach proper "
				+ "termination of the net, specified by 1 final marking. <br/><br/>"
				+ "The algorithm does not guarantee optimal results.</html>";
	}

	public IPNReplayParamProvider constructParamProvider(PluginContext context, PetrinetGraph net, XLog log,
			TransEvClassMapping mapping) {
		return new CostBasedCompleteParamProvider(context, net, log, mapping);
	}

	/**
	 * Return true if all replay inputs are correct
	 */
	public boolean isAllReqSatisfied(PluginContext context, PetrinetGraph net, XLog log, TransEvClassMapping mapping,
			IPNReplayParameter parameter) {

		if (isReqWOParameterSatisfied(context, net, log, mapping)) {
			if (isParameterReqCorrect(net, log, mapping, parameter)) {
				Marking[] finalMarking = ((CostBasedCompleteParam) parameter).getFinalMarkings();
				if ((finalMarking != null) && (finalMarking.length == 1)) {
					return true;
				}
			}
		}

		return false;
	}

	/**
	 * Return true if input of replay without parameters are correct
	 */
	public boolean isReqWOParameterSatisfied(PluginContext context, PetrinetGraph net, XLog log,
			TransEvClassMapping mapping) {
		//		if ((net instanceof ResetInhibitorNet) || (net instanceof InhibitorNet) || (net instanceof ResetNet)
		//				|| (net instanceof Petrinet) || (net instanceof OpenNet)) {
		if ((net instanceof Petrinet)) {
			// check number of transitions, places, and event classes, should be less than Short.MAX_VALUE
			if ((net.getTransitions().size() < Short.MAX_VALUE) && (net.getPlaces().size() < Short.MAX_VALUE)) {
				// check the number of event classes, should be less than Short.MAX_VALUE
				XLogInfo summary = XLogInfoFactory.createLogInfo(log, mapping.getEventClassifier());
				return (summary.getEventClasses().getClasses().size() < Short.MAX_VALUE);
			}
		}
		return false;
	}

	/**
	 * Return true if all replay inputs are correct: parameter type is correct
	 * and non empty (no null); all transitions are mapped to cost; all event
	 * classes (including dummy event class, i.e. an event class that does not
	 * exist in log, any transitions that are NOT silent and not mapped to any
	 * event class in the log is mapped to it) are mapped to cost; all costs
	 * should be non negative; numStates is non negative
	 */
	public boolean isParameterReqCorrect(PetrinetGraph net, XLog log, TransEvClassMapping mapping,
			IPNReplayParameter parameter) {
		if (parameter instanceof CostBasedCompleteParam) {
			CostBasedCompleteParam param = (CostBasedCompleteParam) parameter;
			if ((param.getMapTrans2Cost() != null) && (param.getMaxNumOfStates() != null)
					&& (param.getMapEvClass2Cost() != null) && (param.getInitialMarking() != null)
					&& (param.getFinalMarkings() != null)) {
				// check all transitions are indeed mapped to cost
				if ((param.getMaxNumOfStates() >= 0)
						&& (param.getMapTrans2Cost().keySet().containsAll(net.getTransitions()))) {
					Set<XEventClass> evClassWithCost = param.getMapEvClass2Cost().keySet();
					// check all event classes are mapped to cost
					XEventClassifier classifier = mapping.getEventClassifier();
					XLogInfo summary = XLogInfoFactory.createLogInfo(log, classifier);
					XEventClasses eventClassesName = summary.getEventClasses();

					if (evClassWithCost.containsAll(eventClassesName.getClasses())) {

						// all cost should be non negative
						for (Integer costVal : param.getMapEvClass2Cost().values()) {
							if (costVal < 0) {
								return false;
							}
						}
						for (Integer costVal : param.getMapTrans2Cost().values()) {
							if (costVal < 0) {
								return false;
							}
						}
						return true;
					}
				}
			}
		}
		return false;
	}

	public String toString() {
		return "A Heuristic Cost-based Fitness with ILP, assuming at most " + Short.MAX_VALUE
				+ " tokens in each place.";
	}

	protected void importParameters(CostBasedCompleteParam parameters) {
		// replay parameters
		mapTrans2Cost = parameters.getMapTrans2Cost();
		maxNumOfStates = parameters.getMaxNumOfStates();
		mapEvClass2Cost = parameters.getMapEvClass2Cost();
		mapSync2Cost = parameters.getMapSync2Cost();
		initialMarking = parameters.getInitialMarking();
		finalMarking = parameters.getFinalMarkings()[0];
		usePartialOrderEvents = parameters.isPartiallyOrderedEvents();
	}

}
