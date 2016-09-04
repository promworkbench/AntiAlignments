package org.processmining.antialignments.alignments;

import gnu.trove.iterator.TObjectByteIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.TShortList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectByteHashMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gnu.trove.map.hash.TShortObjectHashMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nl.tue.astar.AStarException;
import nl.tue.astar.util.LPMatrix.LPMatrixException;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.classification.XEventClasses;
import org.deckfour.xes.classification.XEventClassifier;
import org.deckfour.xes.info.XLogInfo;
import org.deckfour.xes.info.XLogInfoFactory;
import org.deckfour.xes.model.XEvent;
import org.deckfour.xes.model.XLog;
import org.deckfour.xes.model.XTrace;
import org.processmining.framework.plugin.PluginContext;
import org.processmining.framework.plugin.annotations.KeepInProMCache;
import org.processmining.framework.util.Pair;
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
import org.processmining.plugins.petrinet.replayresult.StepTypes;
import org.processmining.plugins.replayer.replayresult.SyncReplayResult;

@KeepInProMCache
@PNReplayAlgorithm
public class HeuristicPNetReplayerAlgorithm implements IPNReplayAlgorithm {

	private Map<Transition, Integer> mapTrans2Cost;
	private Map<XEventClass, Integer> mapEvClass2Cost;
	private Map<Transition, Integer> mapSync2Cost;
	private Integer maxNumOfStates;
	private Marking initialMarking;
	private Marking finalMarking;
	private boolean usePartialOrderEvents;
	private XEventClassifier classifier;
	private TShortObjectMap<XEventClass> short2label;
	private TObjectShortMap<XEventClass> label2short;
	private TransEvClassMapping mapping;
	private XLog xLog;
	private XEventClasses classes;
	private Representative[] log2xLog;

	private static class LookupMap extends TObjectByteHashMap<Representative> {

		public LookupMap(int size) {
			super(size);
		}

		public Representative getKeyIfPresent(Representative key) {
			int i = index(key);
			if (i == -1) {
				return null;
			} else {
				return (Representative) _set[i];
			}
		}
	}

	private static class Representative {
		private final int number;
		private final TIntList represented = new TIntArrayList(10);
		private final TShortList trace;

		public Representative(TShortList trace, int number) {
			this.trace = trace;
			this.number = number;
		}

		public void addRepresentedTrace(int trace) {
			represented.add(trace);
		}

		public TIntList getRepresented() {
			return represented;
		}

		public TShortList getTrace() {
			return trace;
		}

		public boolean equals(Object o) {
			return o != null && (o instanceof Representative ? ((Representative) o).trace.equals(trace) : false);
		}

		public int hashCode() {
			return trace.hashCode();
		}

		public String toString() {
			return trace.toString() + " representing " + represented.toString();
		}

		public int getNumber() {
			return number;
		}
	}

	public PNRepResult replayLog(PluginContext context, PetrinetGraph net, XLog xLog, TransEvClassMapping mapping,
			IPNReplayParameter parameters) throws AStarException {

		context.getProgress().setMaximum(xLog.size() + 1);
		this.xLog = xLog;

		this.mapping = mapping;
		importParameters((CostBasedCompleteParam) parameters);
		classifier = mapping.getEventClassifier();
		final XLogInfo summary = XLogInfoFactory.createLogInfo(xLog, classifier);
		this.classes = summary.getEventClasses();

		LookupMap tempLog = new LookupMap(xLog.size());
		short2label = new TShortObjectHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);
		label2short = new TObjectShortHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);

		short c = 0;
		int number = 0;
		for (int t = 0; t < xLog.size(); t++) {
			XTrace trace = xLog.get(t);
			TShortList list = new TShortArrayList(trace.size());
			for (XEvent event : trace) {
				XEventClass clazz = classes.getClassOf(event);
				short id = label2short.putIfAbsent(clazz, c);
				if (id == label2short.getNoEntryValue()) {
					short2label.put(c, clazz);
					id = c;
					c++;
				}
				list.add(id);
			}
			Representative rep = new Representative(list, number);
			Representative existing = tempLog.getKeyIfPresent(rep);
			if (null == existing) {
				tempLog.put(rep, (byte) 1);
				number++;
			} else {
				rep = existing;
			}
			rep.addRepresentedTrace(t);
		}

		for (Transition t : net.getTransitions()) {
			XEventClass clazz = mapping.get(t);
			short id = label2short.putIfAbsent(clazz, c);
			if (id == label2short.getNoEntryValue()) {
				short2label.put(c, clazz);
				id = c;
				c++;
			}
		}

		int t = 0;
		short[][] log = new short[tempLog.size()][];
		log2xLog = new Representative[tempLog.size()];
		TObjectByteIterator<Representative> it = tempLog.iterator();
		while (it.hasNext()) {
			it.advance();
			log[it.key().getNumber()] = it.key().getTrace().toArray();
			log2xLog[it.key().getNumber()] = it.key();
			t++;
		}

		AlignmentILPCalculator calculator = new AlignmentILPCalculator(net, initialMarking, finalMarking, label2short,
				short2label, mapping, log, mapTrans2Cost, mapEvClass2Cost, mapSync2Cost);
		calculator.setLPSolve();
		//		try {
		//						calculator.doExperiment(initialMarking, finalMarking);
		//		} catch (LPMatrixException e) {
		//			// TODO Auto-generated catch block
		//			e.printStackTrace();
		//		}
		calculator.setCutOffLength(7);
		calculator.setMinEvents(2);
		double minCost;
		try {
			TIntList moves = calculator.getAlignmentWithoutTrace(initialMarking, finalMarking);
			minCost = calculator.getCost(moves, new short[0]);
		} catch (LPMatrixException e) {
			minCost = 0.0;
		}
		context.getProgress().inc();

		List<SyncReplayResult> results = new ArrayList<>(log.length);
		for (int tr = 0; tr < log.length; tr++) {
			try {
				long start = System.currentTimeMillis();
				TIntList moves = calculator.getAlignment(initialMarking, finalMarking, tr);
				boolean reliable = calculator.checkFiringSequence(moves, initialMarking, finalMarking);
				long end = System.currentTimeMillis();

				Representative rep = log2xLog[tr];
				//				XTrace trace = xLog.get(rep.getRepresented().get(0));
				SyncReplayResult srr = getSyncReplayResult(calculator, moves, rep.getRepresented().get(0), minCost,
						(int) (end - start), reliable);
				for (int xt = 1; xt < rep.getRepresented().size(); xt++) {
					srr.addNewCase(rep.getRepresented().get(xt));

				}

				results.add(srr);
			} catch (LPMatrixException e) {
			}
			context.getProgress().inc();

		}

		PNRepResult result = new PNRepResult(results);

		return result;
	}

	protected SyncReplayResult getSyncReplayResult(AlignmentILPCalculator calculator, TIntList moves, int xTrace,
			double minCostMoveModel, int time, boolean reliable) {
		List<StepTypes> stepTypes = new ArrayList<StepTypes>(moves.size());
		List<Object> nodeInstance = new ArrayList<Object>(moves.size());

		double lmCost = 0;
		double mmCost = 0;
		double lmUpper = 0;
		double mmUpper = 0;
		double smCost = 0;
		for (int i = 0; i < moves.size(); i++) {
			Pair<Transition, Short> p = calculator.toPair(moves.get(i));

			if (p.getSecond() == null) {
				// Model Move
				nodeInstance.add(p.getFirst());
				if (p.getFirst().isInvisible()) {
					stepTypes.add(StepTypes.MINVI);
				} else {
					stepTypes.add(StepTypes.MREAL);
				}
				mmCost += calculator.getCost(p.getFirst(), null);
				mmUpper += calculator.getCost(p.getFirst(), null);
			} else if (p.getFirst() == null) {
				// Log Move
				XEventClass clazz = short2label.get(p.getSecond());

				stepTypes.add(StepTypes.L);
				nodeInstance.add(clazz);
				lmCost += calculator.getCost(null, clazz);
				lmUpper += calculator.getCost(null, clazz);

			} else {
				// Sync move
				XEventClass clazz = short2label.get(p.getSecond());

				stepTypes.add(StepTypes.LMGOOD);
				nodeInstance.add(p.getFirst());
				smCost += calculator.getCost(p.getFirst(), clazz);
				mmUpper += calculator.getCost(p.getFirst(), null);
				lmUpper += calculator.getCost(null, clazz);

			}

		}

		SyncReplayResult res = new SyncReplayResult(nodeInstance, stepTypes, xTrace);

		res.setReliable(reliable);
		Map<String, Double> info = new HashMap<String, Double>();
		info.put(PNRepResult.RAWFITNESSCOST, (mmCost + lmCost + smCost));

		if (lmCost > 0) {
			info.put(PNRepResult.MOVELOGFITNESS, 1 - (lmCost / lmUpper));
		} else {
			info.put(PNRepResult.MOVELOGFITNESS, 1.0);
		}

		if (mmCost > 0) {
			info.put(PNRepResult.MOVEMODELFITNESS, 1 - (mmCost / mmUpper));
		} else {
			info.put(PNRepResult.MOVEMODELFITNESS, 1.0);
		}
		info.put(PNRepResult.NUMSTATEGENERATED, (double) calculator.steps);
		info.put(PNRepResult.QUEUEDSTATE, 0.0);

		// set info fitness
		if (mmCost > 0 || lmCost > 0 || smCost > 0) {
			info.put(PNRepResult.TRACEFITNESS, 1 - ((mmCost + lmCost + smCost) / (lmUpper + minCostMoveModel)));
		} else {
			info.put(PNRepResult.TRACEFITNESS, 1.0);
		}
		info.put(PNRepResult.TIME, (double) time);
		info.put(PNRepResult.ORIGTRACELENGTH, (double) xLog.get(xTrace).size());
		res.setInfo(info);
		return res;

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
