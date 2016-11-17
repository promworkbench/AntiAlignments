package org.processmining.antialignments.ilp.alignment;

import gnu.trove.list.TIntList;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nl.tue.astar.AStarException;
import nl.tue.astar.util.ilp.LPMatrixException;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.classification.XEventClasses;
import org.deckfour.xes.classification.XEventClassifier;
import org.deckfour.xes.extension.std.XConceptExtension;
import org.deckfour.xes.factory.XFactoryRegistry;
import org.deckfour.xes.info.XLogInfo;
import org.deckfour.xes.info.XLogInfoFactory;
import org.deckfour.xes.model.XLog;
import org.deckfour.xes.model.XTrace;
import org.processmining.antialignments.ilp.util.AbstractHeuristicILPReplayer;
import org.processmining.antialignments.ilp.util.Representative;
import org.processmining.antialignments.test.AlignmentSetup;
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
import org.processmining.plugins.petrinet.replayer.annotations.PNReplayAlgorithm;
import org.processmining.plugins.petrinet.replayresult.PNRepResult;
import org.processmining.plugins.petrinet.replayresult.PNRepResultImpl;
import org.processmining.plugins.petrinet.replayresult.StepTypes;
import org.processmining.plugins.replayer.replayresult.SyncReplayResult;

@KeepInProMCache
@PNReplayAlgorithm
public class HeuristicPNetReplayerExperiment extends AbstractHeuristicILPReplayer<Petrinet> implements
		IPNReplayAlgorithm {

	static final String SOLVER = "ILP Solver";
	static final String TIME = "Total time (ms)";
	static final String EXPECTEDMOVES = "User expected model moves";

	static final String CUTOFF = "Parameter: window size";
	static final String MINEVENT = "Parameter: Minimal events in window";
	static final String BACKTRACKBOUND = "Parameter: Backtrack bound";
	static final String BACKTRACKTHRESHOLD = "Parameter: Backtrack threshold";

	static final String MINMODELMOVECOST = "Minimal model move cost";

	private final static int MAXTRACE = 100;

	private Map<Transition, Integer> mapTrans2Cost;
	private Map<XEventClass, Integer> mapEvClass2Cost;
	private Map<Transition, Integer> mapSync2Cost;

	private int expectedModelMovesParameter;
	private int backtrackLimitParameter;
	private double backtrackThresholdParameter;
	private int estRows;
	private int estColumns;

	public PNRepResultImpl replayLog(PluginContext context, PetrinetGraph net, XLog xLog, TransEvClassMapping mapping,
			IPNReplayParameter parameters) throws AStarException {

		importParameters((HeuristicParameters) parameters);

		setUpDataStructures((Petrinet) net, parameters.getInitialMarking(),
				((HeuristicParameters) parameters).getFinalMarkings()[0], xLog, mapping);

		context.getProgress().setMaximum(log.length + 1);

		// Set parameters just over the bounds for the ILP's
		int cutOffEvent = expectedModelMovesParameter + 1;
		int minEvent = 1;

		PrintStream out;
		try {
			out = new PrintStream(new File("d:/temp/antialignment/log.csv"));
		} catch (FileNotFoundException e1) {
			out = System.out;
		}
		String sep = ";";
		PNRepResultExportPlugin export = new PNRepResultExportPlugin();
		export.printResult(out, null, -1, -1, sep);

		context.log("Starting replay with " + cutOffEvent + " exact variables containing at least " + minEvent
				+ " labeled moves.");

		PNRepResultImpl result = null;

		context.getProgress().setMaximum(
				cutOffEvent * ((int) (Math.log(cutOffEvent) / Math.log(2) + 0.5)) * log.length + cutOffEvent);

		long startWhole = System.currentTimeMillis();

		// For BenchMarking Purposes, go through the log one trace at the time
		PNRepResult[][] realAlignments = new PNRepResult[Math.min(MAXTRACE, log.length)][2];
		XLog newLog = XFactoryRegistry.instance().currentDefault().createLog();
		newLog.add(XFactoryRegistry.instance().currentDefault().createTrace());

		// Setup the alignmentAlgorithm
		AlignmentSetup alignmentAlgorithm = new AlignmentSetup((Petrinet) net, newLog, mapping, mapTrans2Cost,
				mapEvClass2Cost);

		// Compute Alignments trace by trace
		context.log("Starting real replay on empty log with ILP.");
		PNRepResult res = alignmentAlgorithm.getAlignment(context, initialMarking, finalMarking, true, 0);
		int minCost = (int) Math.round(((Double) res.getInfo().get(PNRepResult.RAWFITNESSCOST)));

		for (int t = 0; t < realAlignments.length; t++) {

			newLog = XFactoryRegistry.instance().currentDefault().createLog();
			int tr = log2xLog[t].getNumber();
			newLog.add((XTrace) xLog.get(tr).clone());

			// Setup the alignmentAlgorithm
			alignmentAlgorithm = new AlignmentSetup((Petrinet) net, newLog, mapping, mapTrans2Cost, mapEvClass2Cost);

			// Compute Alignments trace by trace
			context.log("Starting real replay on trace " + t + " without ILP.");
			realAlignments[t][0] = alignmentAlgorithm.getAlignment(context, initialMarking, finalMarking, false,
					minCost);
			realAlignments[t][0].iterator().next().getTraceIndex().clear();
			realAlignments[t][0].iterator().next().getTraceIndex().add(tr);

			context.log("Starting real replay on trace " + t + " with ILP.");
			realAlignments[t][1] = alignmentAlgorithm
					.getAlignment(context, initialMarking, finalMarking, true, minCost);
			realAlignments[t][1].iterator().next().getTraceIndex().clear();
			realAlignments[t][1].iterator().next().getTraceIndex().add(tr);

			realAlignments[t][0].getInfo().put(HeuristicPNetReplayerAlgorithm.SOLVER, "A*");
			export.printResult(out, realAlignments[t][0], -1, -1, sep);
			realAlignments[t][1].getInfo().put(HeuristicPNetReplayerAlgorithm.SOLVER, "A*_ILP");
			export.printResult(out, realAlignments[t][1], -1, -1, sep);

			out.flush();

			// Starting full-blown experiment
			AlignmentILPCalculator calculator = new AlignmentILPCalculator(this.net, initialMarking, finalMarking,
					label2short, short2label, mapping, new short[][] { log[t] }, mapTrans2Cost, mapEvClass2Cost,
					mapSync2Cost);
			//DEBUGCODE
			calculator.VERBOSE = false;
			calculator.NAMES = false;

			for (int c = 1; c < cutOffEvent; c++) {
				for (int b = c; b > 0; b /= 2) {
					//			for (int b = 1; b <= c; b++) {
					for (int solver = 0; solver < 2; solver++) {

						if (solver == 0) {
							calculator.setLPSolve();
						} else {
							calculator.setGurobi();
						}

						result = computeAlignments(context, calculator, minCost, c, 0, b, backtrackThresholdParameter,
								t, t + 1);
						export.printResult(out, result, estRows, estColumns, sep);
						for (int e = 1; e < c; e *= 2) {
							System.out.print(".");

							result = computeAlignments(context, calculator, minCost, c, e, b,
									backtrackThresholdParameter, t, t + 1);
							export.printResult(out, result, estRows, estColumns, sep);
						}
						System.out.print(".");
						result = computeAlignments(context, calculator, minCost, c, c, b, backtrackThresholdParameter,
								t, t + 1);
						export.printResult(out, result, estRows, estColumns, sep);
						System.out.println("!");
					}
				}

			}

		}

		out.close();

		long endWhole = System.currentTimeMillis();

		result.addInfo(PNRepResult.VISTITLE, "Heuristic Alignments of "
				+ XConceptExtension.instance().extractName(xLog) + " on " + net.getLabel());
		//		result.addInfo(CBOUND, Integer.toString(cBound));
		//		result.addInfo(RBOUND, Integer.toString(rBound));
		result.addInfo(TIME, Integer.toString((int) (endWhole - startWhole)));
		return result;

	}

	private boolean switchToGurobi(PluginContext context, AlignmentILPCalculator calculator, int cutOffEvent,
			int minEvent) {
		calculator.setLPSolve();
		boolean gurobi = false;

		// estimate number of columns and rows (This is not exact!)
		estColumns = (net.getTransitions().size() + 2 * minEvent) * cutOffEvent + 2 * net.getTransitions().size()
				+ label2short.size();
		estRows = (net.getPlaces().size() + minEvent) * cutOffEvent + net.getPlaces().size() + label2short.size()
				+ cutOffEvent + 2;

		//		// the estimated minimum number of rows and colums are known.
		//		while (rows < rBound && columns < cBound) {
		//			// Try to increase without loosing too much CPU time in LpSolve
		//			cutOffEvent += expectedModelMoves;
		//			minEvent++;
		//			columns = (net.getTransitions().size() + 2 * minEvent) * cutOffEvent + 2 * net.getTransitions().size()
		//					+ label2short.size();
		//			rows = (net.getPlaces().size() + 2 * minEvent) * cutOffEvent + net.getPlaces().size() + label2short.size();
		//		}

		// But what if that fails?

		if (calculator.setGurobi()) {
			gurobi = true;
			context.log("Solver set to Gurobi, because of problem size");
		} else if (minEvent > 1) {
			context.log("Failed to load Gurobi...");
			// cannot setup gurobi. Nothing we can do, but resort to smallest case.
			minEvent = 1;
			cutOffEvent = 1 + expectedModelMovesParameter;
		}

		context.log("Estimated ILP size: " + estColumns + " columns and " + estRows + " rows.");

		return gurobi;
	}

	public PNRepResultImpl computeAlignments(PluginContext context, AlignmentILPCalculator calculator, double minCost,
			int cutOffEvent, int minEvent, int backtrackLimit, double backtrackThreshold, int trFrom, int trTo) {
		List<SyncReplayResult> results = new ArrayList<>(trTo - trFrom);
		for (int tr = trFrom; tr < trTo && !context.getProgress().isCancelled(); tr++) {
			//		for (int tr = log.length; tr-- > 0 && !context.getProgress().isCancelled();) {
			try {
				calculator.setMinEvents(Math.min(minEvent, log[tr].length));
				calculator.setCutOffLength(cutOffEvent);
				calculator.setBacktrackLimit(backtrackLimit);
				calculator.setBacktrackThreshold(backtrackThreshold);

				context.log("Starting replay of trace " + tr + "/" + log.length + " with " + cutOffEvent
						+ " exact variables containing at least " + minEvent + " labeled moves.");

				long start = System.currentTimeMillis();
				TIntList moves = calculator.getAlignment(context.getProgress(), initialMarking, finalMarking, tr
						- trFrom);
				boolean reliable = calculator.checkAndReorderFiringSequence(moves, initialMarking, finalMarking, true);
				reliable &= calculator.checkTrace(moves, log[tr], true);
				long end = System.currentTimeMillis();

				Representative rep = log2xLog[tr];
				SyncReplayResult srr = getSyncReplayResult(calculator, moves, rep.getRepresented().get(0), minCost,
						(int) (end - start), reliable);
				for (int xt = 1; xt < rep.getRepresented().size(); xt++) {
					srr.addNewCase(rep.getRepresented().get(xt));

				}

				results.add(srr);
			} catch (LPMatrixException _) {
			}
			context.getProgress().inc();

		}

		PNRepResultImpl result = new PNRepResultImpl(results);

		result.addInfo(EXPECTEDMOVES, Integer.toString(expectedModelMovesParameter));
		result.addInfo(CUTOFF, Integer.toString(cutOffEvent));
		result.addInfo(MINEVENT, Integer.toString(minEvent));
		result.addInfo(BACKTRACKBOUND, Integer.toString(backtrackLimit));
		result.addInfo(BACKTRACKTHRESHOLD, Double.toString(backtrackThreshold));
		result.addInfo(SOLVER, calculator.isGurobi() ? "Gurobi" : "LpSolve");

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
		return new HeuristicParameterProvider(context, net, log, mapping);
	}

	/**
	 * Return true if all replay inputs are correct
	 */
	public boolean isAllReqSatisfied(PluginContext context, PetrinetGraph net, XLog log, TransEvClassMapping mapping,
			IPNReplayParameter parameter) {

		if (isReqWOParameterSatisfied(context, net, log, mapping)) {
			if (isParameterReqCorrect(net, log, mapping, parameter)) {
				Marking[] finalMarking = ((HeuristicParameters) parameter).getFinalMarkings();
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
		if (parameter instanceof HeuristicParameters) {
			HeuristicParameters param = (HeuristicParameters) parameter;
			if ((param.getMapTrans2Cost() != null) && (param.getMapEvClass2Cost() != null)
					&& (param.getInitialMarking() != null) && (param.getFinalMarkings() != null)) {
				// check all transitions are indeed mapped to cost
				if ((param.getMapTrans2Cost().keySet().containsAll(net.getTransitions()))) {
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
		return "A Heuristic Cost-based Fitness Experiment, assuming at most " + Short.MAX_VALUE
				+ " tokens in each place.";
	}

	protected void importParameters(HeuristicParameters parameters) {
		// replay parameters
		mapTrans2Cost = parameters.getMapTrans2Cost();
		//		maxNumOfStates = parameters.getMaxNumOfStates();
		mapEvClass2Cost = parameters.getMapEvClass2Cost();
		mapSync2Cost = parameters.getMapSync2Cost();

		expectedModelMovesParameter = parameters.getExpecteModelMoves();
		backtrackLimitParameter = parameters.getBacktrackLimit();
		backtrackThresholdParameter = parameters.getBacktrackThreshold();
	}

}
