package org.processmining.antialignments;

import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gnu.trove.map.hash.TShortObjectHashMap;
import gnu.trove.set.TShortSet;
import gnu.trove.set.hash.TShortHashSet;
import gurobi.GRBException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Vector;

import lpsolve.LpSolveException;
import nl.tue.astar.AStarException;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.model.XLog;
import org.processmining.antialignments.algorithm.AntiAlignmentILPCalculator;
import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.contexts.uitopia.UIPluginContext;
import org.processmining.contexts.uitopia.annotations.UITopiaVariant;
import org.processmining.framework.connections.ConnectionCannotBeObtained;
import org.processmining.framework.connections.ConnectionManager;
import org.processmining.framework.plugin.annotations.Plugin;
import org.processmining.framework.plugin.annotations.PluginLevel;
import org.processmining.framework.plugin.annotations.PluginVariant;
import org.processmining.models.connections.petrinets.EvClassLogPetrinetConnection;
import org.processmining.models.connections.petrinets.behavioral.FinalMarkingConnection;
import org.processmining.models.connections.petrinets.behavioral.InitialMarkingConnection;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;
import org.processmining.plugins.petrinet.replayer.PNLogReplayer;
import org.processmining.plugins.petrinet.replayresult.PNRepResult;
import org.processmining.plugins.petrinet.replayresult.StepTypes;
import org.processmining.plugins.replayer.replayresult.SyncReplayResult;

@Plugin(name = "Anti-Alignment Precision/Generalization", level = PluginLevel.NightlyBuild, //
returnLabels = { "Anti-alignment-based Precision/Generalization" }, returnTypes = { AntiAlignmentValues.class },//
parameterLabels = { "Petri net", "Event Log", "Alignment" }, //
help = "Measure precision/generalization using anti alignments.", userAccessible = true)
public class AntiAlignmentPlugin {

	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	pack = "AntiAlignments")
	@PluginVariant(variantLabel = "Without alignments", requiredParameterLabels = { 0, 1 })
	public AntiAlignmentValues measurePrecisionWithoutAlignments(final UIPluginContext context, Petrinet net, XLog log) {
		// replay log on model (or obtain existing replay result)
		PNRepResult alignments;
		try {
			PNLogReplayer replayer = new PNLogReplayer();
			alignments = replayer.replayLog(context, net, log);

			return measurePrecision(context, net, log, alignments);

		} catch (ConnectionCannotBeObtained e) {
			context.log("No connection between the given net and log. For computing anti-alignment based precision, the plugin"
					+ " needs a Petri net with an initial and final marking!");
		} catch (AStarException e) {
			context.log("Error in replay algorithm: " + e.getMessage());
		}

		return null;
	}

	@UITopiaVariant(affiliation = UITopiaVariant.EHV, author = "Boudewijn van Dongen", email = "b.f.v.dongen@tue.nl", //
	pack = "AntiAlignments")
	@PluginVariant(variantLabel = "With alignments", requiredParameterLabels = { 0, 1, 2 })
	public AntiAlignmentValues measurePrecision(final UIPluginContext context, Petrinet net, XLog log,
			PNRepResult alignments) {
		// retrieve mapping between process model to log
		try {
			ConnectionManager connManager = context.getConnectionManager();
			EvClassLogPetrinetConnection conn = connManager.getFirstConnection(EvClassLogPetrinetConnection.class,
					context, net, log);
			TransEvClassMapping mapping = (TransEvClassMapping) conn
					.getObjectWithRole(EvClassLogPetrinetConnection.TRANS2EVCLASSMAPPING);

			// get marking
			InitialMarkingConnection initMarkingConn = connManager.getFirstConnection(InitialMarkingConnection.class,
					context, net);
			Marking initialMarking = initMarkingConn.getObjectWithRole(InitialMarkingConnection.MARKING);

			FinalMarkingConnection finalMarkingConn = connManager.getFirstConnection(FinalMarkingConnection.class,
					context, net);
			Marking finalMarking = finalMarkingConn.getObjectWithRole(FinalMarkingConnection.MARKING);

			return basicCodeStructureWithAlignments(net, initialMarking, finalMarking, log, alignments, mapping);

		} catch (ConnectionCannotBeObtained noConnection) {
			context.log("No connection between the given net and log. For computing anti-alignment based precision, the plugin"
					+ " needs a Petri net with an initial and final marking!");
		}

		return null;
	}

	public static AntiAlignmentValues basicCodeStructureWithoutAlignments(Petrinet net, Marking initialMarking,
			Marking finalMarking, XLog xLog) {

		// Setup the alignmentAlgorithm
		AlignmentSetup alignmentAlgorithm = new AlignmentSetup(net, xLog);

		// Compute Alignments
		PNRepResult alignments = alignmentAlgorithm.getAlignment(initialMarking, finalMarking, false);

		return basicCodeStructureWithAlignments(net, initialMarking, finalMarking, xLog, alignments,
				alignmentAlgorithm.mapping);

	}

	public static AntiAlignmentValues basicCodeStructureWithAlignments(Petrinet net, Marking initialMarking,
			Marking finalMarking, XLog xLog, PNRepResult alignments, TransEvClassMapping mapping) {

		// Obtain the fitness:
		double traceFitness = (Double) alignments.getInfo().get(PNRepResult.TRACEFITNESS);

		// Build the aligned event log for anti-alignment computations
		int[] frequencies = new int[alignments.size()];
		short[][] alignedLog = new short[alignments.size()][];
		List[] firingSequences = new List[alignments.size()];

		// We need to rebuild the log using the alignments. The frequencies and alignedlog objects
		// are filled with fitting traces and firing sequences (including invisible transitions) and 
		// the label2short map is built for future use. Max indicates the trace length.
		TObjectShortMap<String> label2short = new TObjectShortHashMap<>();
		TShortObjectMap<String> short2label = new TShortObjectHashMap<>();

		int max = rebuildLogFromAlignments(alignments, mapping, frequencies, alignedLog, firingSequences, label2short,
				short2label);

		//		AntiAlignmentCalculator calculator = new AntiAlignmentCalculator(net, initialMarking, finalMarking, label2short);

		// Start anti-alignment computation
		int maxFactor = 1;
		max *= 2 * maxFactor;

		AntiAlignmentILPCalculator calculator2 = null;
		try {
			calculator2 = new AntiAlignmentILPCalculator(net, initialMarking, finalMarking, label2short, short2label,
					alignedLog, max, maxFactor);
		} catch (LpSolveException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		long start = System.nanoTime();

		AntiAlignments aa3 = null;
		try {
			aa3 = calculator2.getAntiAlignments();
		} catch (LpSolveException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		long mid = System.nanoTime();
		System.out.println("new: " + (mid - start) / 1000000.0 + " ms");

		mid = System.nanoTime();
		//		AntiAlignments aa2 = calculator.getAntiAlignments(alignedLog, max, maxFactor);
		AntiAlignments aa = new DepthFirstTraceSearch(net, initialMarking, finalMarking, label2short)
				.getAntiAlignments(alignedLog, max, maxFactor, new DistanceMetric.Hamming());

		long end = System.nanoTime();

		System.out.println("old: " + (end - mid) / 1000000.0);

		System.out.println("old: " + Arrays.toString(aa.getMaxMinDistances()));
		System.out.println("new: " + Arrays.toString(aa3.getMaxMinDistances()));
		for (int i = 0; i < aa.getMaxMinDistances().length; i++) {
			if (aa.getMaxMinDistances()[i] != aa3.getMaxMinDistances()[i]) {
				System.out.print("OLD: ");
				System.out.println(toString(aa.getAntiAlignments()[i], short2label));
				System.out.print("NEW: ");
				System.out.println(toString(aa3.getAntiAlignments()[i], short2label));
			}
		}

		double[] weightedPrecision = computePrecision(aa, alignedLog, frequencies, max, maxFactor);
		double[] unweightedPrecision = computePrecision(aa, alignedLog, null, max, maxFactor);

		// now for the generalization part.
		// compute the states that were visited by the aligned log, using the
		// firing sequences in the aligned log
		Map<Marking, TShortSet> statesVisitedPerTrace = getStatesVisitedPerTrace(net, initialMarking, firingSequences);

		int[] count = new int[alignedLog.length + 1];
		double[] recDist = new double[alignedLog.length + 1];

		int[] cAndRD = countNewStatesAndRecoveryDistance(net, initialMarking, aa.getAAFiringSequenceForLog(),
				statesVisitedPerTrace, (short) -1);
		count[alignedLog.length] = cAndRD[0];
		recDist[alignedLog.length] = cAndRD[1];

		for (short tr = 0; tr < alignedLog.length; tr++) {

			cAndRD = countNewStatesAndRecoveryDistance(net, initialMarking,
					aa.getAAFiringSequenceForLogWithoutTrace(tr), statesVisitedPerTrace, tr);
			count[tr] = cAndRD[0];
			recDist[tr] = cAndRD[1];

		}

		double[] weightedGeneralization = computeGeneralization(aa, alignedLog, frequencies, recDist, count, max,
				maxFactor);
		double[] unweightedGeneralization = computeGeneralization(aa, alignedLog, null, recDist, count, max, maxFactor);

		//		printAntiAlignments("print", aa, alignedLog, frequencies, count, recDist, short2label, true);
		return new AntiAlignmentValues(0.5, traceFitness, unweightedPrecision[0], unweightedPrecision[1],
				weightedGeneralization[0], weightedGeneralization[1]);

	}

	private static int rebuildLogFromAlignments(PNRepResult alignments, TransEvClassMapping mapping, int[] frequencies,
			short[][] alignedLog, List[] firingSequences, TObjectShortMap<String> label2short,
			TShortObjectMap<String> short2label) {
		// move through alignments and rebuild event log to an aligned version with short labels.
		// We need to build both an aligned log and aligned firing sequences. The latter includes
		// invisible transitions which are required to compute the visited states later.

		// What is needed here is a consistent mapping from labels in the log to shorts and 
		// labels in the model to shorts. This mapping has been built internally by the replayer
		// but is currently not accessible...
		// Ugly, but map is from both XEventClass and Transition to short

		// add all mapped objects first
		short mapped = 1;
		for (Entry<Transition, XEventClass> tc : mapping.entrySet()) {
			label2short.put(tc.getKey().getLabel(), mapped);
			label2short.put(tc.getValue().toString(), mapped);
			short2label.put(mapped, tc.getKey().getLabel());
			mapped++;
		}

		int max = 0;
		int i = 0;
		for (SyncReplayResult alignment : alignments) {
			frequencies[i] = alignment.getTraceIndex().size();

			firingSequences[i] = new ArrayList<Transition>();
			TShortList modelSeq = new TShortArrayList();
			for (int s = 0; s < alignment.getStepTypes().size(); s++) {
				Transition t = null;
				XEventClass c = null;
				short m = label2short.getNoEntryValue();
				if (alignment.getStepTypes().get(s) == StepTypes.LMGOOD) {
					// synchronous move
					// Corresponding nodeStep is a transition.
					t = (Transition) alignment.getNodeInstance().get(s);
					c = mapping.get(t);
					m = label2short.get(t.getLabel());
				} else if (alignment.getStepTypes().get(s) == StepTypes.L) {
					// log move
					// Corresponding nodeStep is an event class.
					c = (XEventClass) alignment.getNodeInstance().get(s);
					if (label2short.get(c.toString()) == label2short.getNoEntryValue()) {
						label2short.put(c.toString(), mapped);
						short2label.put(mapped, c.toString());
						mapped++;
					}
				} else if (alignment.getStepTypes().get(s) == StepTypes.MREAL) {
					// Model move on visible transition
					// Corresponding nodeStep is a transition.
					t = (Transition) alignment.getNodeInstance().get(s);
					m = label2short.get(t.getLabel());
					if (m == label2short.getNoEntryValue()) {
						label2short.put(t.getLabel(), mapped);
						short2label.put(mapped, t.getLabel());
						m = mapped;
						mapped++;
					}
				} else if (alignment.getStepTypes().get(s) == StepTypes.MINVI) {
					// model move on invisible transition
					// Corresponding nodeStep is a transition.
					t = (Transition) alignment.getNodeInstance().get(s);
					if (label2short.get(t.getLabel()) == label2short.getNoEntryValue()) {
						label2short.put(t.getLabel(), mapped);
						short2label.put(mapped, t.getLabel());
						mapped++;
					}
				} else {
					System.out.println("error");
				}
				if (t != null) {
					firingSequences[i].add(t);

					if (m != label2short.getNoEntryValue()) {
						modelSeq.add(m);
					}
				}
			}
			alignedLog[i] = modelSeq.toArray();
			if (alignedLog[i].length > max) {
				max = alignedLog[i].length;
			}
			i++;
		}
		// all elements are mapped to shorts and we have constructed the aligned log, both in terms of
		// only mapped sync move elements and in terms of firing sequences.
		return max;
	}

	private static double[] computeGeneralization(AntiAlignments aa, short[][] log, int[] frequencies,
			double[] recDistances, int[] newStateCounts, int max, double maxFactor) {

		double sum = 0;
		int sumFreq = 0;
		int f = 1;
		for (int t = 0; t < aa.getLogLength(); t++) {
			if (frequencies != null) {
				f = frequencies[t];
			}
			// sum += f
			// * Math.max(0, aa.getAADistanceForLogWithoutTrace(t)
			// - recDistances[t])
			// / (double) aa.getMaxDistanceForTrace(t);
			double d = Math.max(aa.getAAForLogWithoutTrace(t).length, aa.getMaxDistanceForTrace(t));
			double aaDistance = d < 1 ? 1 : aa.getAADistanceForLogWithoutTrace(t) / d;
			assert 0 <= aaDistance && aaDistance <= 1;
			// double recDistance = newStateCounts[t] == 0 ? 0 : recDistances[t]
			// / (double) (newStateCounts[t]);
			d = (aa.getAAFiringSequenceForLogWithoutTrace(t).size() - 1);
			double recDistance = d < 1 ? 1 : recDistances[t] / d;
			assert 0 <= recDistance && recDistance <= 1;

			// double d = Math.max(0, aaDistance - recDistance);
			// sum += f * Math.sqrt(2 * d - d * d);
			// sum += f * Math.abs(aaDistance - recDistance);
			// sum += f * (1 - Math.abs(aaDistance + recDistance - 1));

			d = (1 - aaDistance) * (1 - aaDistance);
			d += (recDistance) * (recDistance);
			d = Math.sqrt(d);
			d = Math.min(1, d);
			sum += f * (1 - d);
			// sum += f * d;

			sumFreq += f;

		}

		double d = Math.min(aa.getAAForLog().length, max * maxFactor);
		double aaDistance = d < 1 ? 1 : aa.getAADistanceForLog() / d;
		assert 0 <= aaDistance && aaDistance <= 1;
		// double recDistance = newStateCounts[t] == 0 ? 0 : recDistances[t]
		// / (double) (newStateCounts[t]);
		d = aa.getAAFiringSequenceForLog().size() - 1;
		double recDistance = d < 1 ? 1 : recDistances[log.length] / d;
		d = (1 - aaDistance) * (1 - aaDistance);
		d += recDistance * recDistance;
		d = Math.sqrt(d);
		d = Math.min(1, d);

		return new double[] { (sum / (double) sumFreq), 1 - d };

	}

	private static double[] computePrecision(AntiAlignments aa, short[][] log, int[] frequencies, int max,
			double maxFactor) {

		double sum = 0;
		int sumFreq = 0;
		int f = 1;
		for (int t = 0; t < aa.getLogLength(); t++) {
			if (frequencies != null) {
				f = frequencies[t];
			}
			// the distance to the removed trace is at most the
			// length of the removed trace
			assert (int) (log[t].length * maxFactor) >= aa.getAADistanceToTrace(t);
			// we normalize on the length of the trace * maxFactor as this is
			// the given maximum
			// to the anti-alignment algorithm.
			// sum += f
			// * (aa.getAADistanceToTrace(t) / (double) ((int) (log[t].length *
			// maxFactor)));

			// We can also normalize on the length of the returned
			// anti-alignment!
			double d = Math.min(aa.getAAForLogWithoutTrace(t).length, log[t].length * maxFactor);
			sum += d < 1 ? f : f * (aa.getAADistanceToTrace(t) / d);
			sumFreq += f;

		}
		double d = Math.min(aa.getAAForLog().length, max * maxFactor);
		return new double[] { 1 - sum / (double) sumFreq, d < 1 ? 0 : 1 - aa.getAADistanceForLog() / d };

	}

	public static int[] countNewStatesAndRecoveryDistance(Petrinet net, Marking initialMarking,
			Vector<Transition> firingSequence, Map<Marking, TShortSet> statesVisitedPerTrace, short indexToIgnore) {

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);
		BreadthFirstTraceSearch breadthFirstSearch = new BreadthFirstTraceSearch(net, initialMarking,
				statesVisitedPerTrace);

		int length = firingSequence.size();
		int count = 0;
		int maxRecDist = 0;
		int sumRecDist = 0;
		for (Transition t : firingSequence) {
			try {
				semantics.executeExecutableTransition(t);
				length--;
			} catch (IllegalTransitionException e) {
				e.printStackTrace();
			}
			TShortSet set = statesVisitedPerTrace.get(semantics.getCurrentState());
			if (set == null || (set.size() == 1 && set.contains(indexToIgnore))) {
				count++;
				// state that is not reached in the log. Let's compute recovery
				// distance.
				int recDist = breadthFirstSearch
						.getShortestDistance(semantics.getCurrentState(), length, indexToIgnore);
				if (recDist > maxRecDist) {
					maxRecDist = recDist;
				}
				sumRecDist += recDist;
			}
		}
		return new int[] { count, maxRecDist, sumRecDist, firingSequence.size() };
	}

	/**
	 * Returns a map from markings to a set of trace indices, such that each
	 * trace in the set has visited this particular state.
	 * 
	 * @param net
	 * @param initialMarking
	 * @param firingSequences
	 * @return
	 */
	public static Map<Marking, TShortSet> getStatesVisitedPerTrace(Petrinet net, Marking initialMarking,
			List<Transition>[] firingSequences) {

		Map<Marking, TShortSet> map = new HashMap<>();

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		for (short i = 0; i < firingSequences.length; i++) {
			semantics.setCurrentState(initialMarking);
			if (!map.containsKey(initialMarking)) {
				map.put(initialMarking, new TShortHashSet(firingSequences.length, 0.75f, Short.MIN_VALUE));
			}
			map.get(initialMarking).add(i);
			for (Transition t : firingSequences[i]) {
				try {
					semantics.executeExecutableTransition(t);
				} catch (IllegalTransitionException e) {
					e.printStackTrace();
				}
				if (!map.containsKey(semantics.getCurrentState())) {
					map.put(semantics.getCurrentState(), new TShortHashSet(firingSequences.length, 0.75f,
							Short.MIN_VALUE));
				}
				map.get(semantics.getCurrentState()).add(i);
			}
		}

		return map;
	}

	public static void printAntiAlignments(String model, AntiAlignments aa, short[][] log, int[] frequencies,
			int[] newStateCount, double[] recDistances, TShortObjectMap<String> short2label, boolean printHeader) {
		// Depth first search concluded. Let's report:
		if (printHeader) {
			System.out
					.println("Model      \tLog           \tRemoved Trace\tLength\tFreq\tAnti Alignment\tDlog\tDrem\tnewS\tDrec\tFiring Sequence");
		}

		System.out.print(model);
		System.out.print("\t");
		System.out.print("Whole log\t\t\t\t");
		System.out.print(toString(aa.getAAForLog(), short2label));
		System.out.print("\t");
		System.out.print(aa.getAADistanceForLog());
		System.out.print("\t");
		System.out.print("0\t");
		System.out.print(newStateCount[log.length]);
		System.out.print("\t");
		System.out.print(recDistances[log.length]);
		System.out.print("\t");
		System.out.print(aa.getAAFiringSequenceForLog());
		System.out.println();

		for (int t = 0; t < aa.getLogLength(); t++) {
			System.out.print(model);
			System.out.print("\t");
			System.out.print("Trace " + t + " removed\t");
			System.out.print(toString(log[t], short2label));
			System.out.print("\t");
			System.out.print(log[t].length);
			System.out.print("\t");
			System.out.print(frequencies[t]);
			System.out.print("\t");
			if (aa.getAAForLogWithoutTrace(t) == null) {
				System.out.print("none");
			} else {
				System.out.print(toString(aa.getAAForLogWithoutTrace(t), short2label));
			}
			System.out.print("\t");
			System.out.print(aa.getAADistanceForLogWithoutTrace(t));
			System.out.print("\t");
			System.out.print(aa.getAADistanceToTrace(t));
			System.out.print("\t");
			System.out.print(newStateCount[t]);
			System.out.print("\t");
			System.out.print(recDistances[t]);
			System.out.print("\t");
			if (aa.getAAFiringSequenceForLogWithoutTrace(t) == null) {
				System.out.print("none");
			} else {
				System.out.print(aa.getAAFiringSequenceForLogWithoutTrace(t).toString());
			}
			System.out.println();

		}

	}

	public static String toString(short[] sequence, TShortObjectMap<String> short2label) {
		String s = "<";
		int i;
		for (i = 0; i < sequence.length; i++) {
			s += short2label.get(sequence[i]);
			s += i == sequence.length - 1 ? "" : ",";
		}
		s += ">";

		return s;
	}

}
