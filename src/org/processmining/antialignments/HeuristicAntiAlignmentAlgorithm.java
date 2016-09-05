package org.processmining.antialignments;

import gnu.trove.map.TShortObjectMap;
import gnu.trove.set.TShortSet;
import gnu.trove.set.hash.TShortHashSet;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import nl.tue.astar.util.LPMatrix.LPMatrixException;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.model.XLog;
import org.processmining.antialignments.algorithm.AntiAlignmentILPCalculator;
import org.processmining.antialignments.base.AbstractHeuristicILPReplayer;
import org.processmining.antialignments.base.AlignedRepresentative;
import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.IllegalTransitionException;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;
import org.processmining.plugins.petrinet.replayresult.PNRepResult;
import org.processmining.plugins.petrinet.replayresult.StepTypes;
import org.processmining.plugins.replayer.replayresult.SyncReplayResult;

public class HeuristicAntiAlignmentAlgorithm extends AbstractHeuristicILPReplayer<Petrinet> {

	public static final String HAMMINGDISTANCETOLOG = "Minimal Hamming Distance to Log";
	public static final String HAMMINGDISTANCETOREMOVED = "Hamming Distance to Removed Trace";
	private static final String PRECISION = "Precision";
	private static final String GENERALIZATION = "Generalization";
	private static final String LOGPRECISION = "Log-based Precision";
	private static final String LOGGENERALIZATION = "Log-based Generalization";
	private static final String TRACEPRECISION = "Trace-based Precision";
	private static final String TRACEGENERALIZATION = "Trace-based Generalization";

	private Double traceFitness;

	public HeuristicAntiAlignmentAlgorithm(Petrinet net, Marking initialMarking, Marking finalMarking, XLog xLog,
			PNRepResult alignments, TransEvClassMapping mapping) {
		setUpDataStructures(net, initialMarking, finalMarking, xLog, alignments, mapping);

		// Obtain the fitness:
		traceFitness = (Double) alignments.getInfo().get(PNRepResult.TRACEFITNESS);

		// Build the aligned event log for anti-alignment computations
	}

	public AntiAlignments computeAntiAlignments() {
		// Start anti-alignment computation
		int maxFactor = 1;
		int max = maxTraceLength * 2 * maxFactor;

		AntiAlignmentILPCalculator calculator2 = null;
		calculator2 = new AntiAlignmentILPCalculator(net, initialMarking, finalMarking, label2short, short2label,
				mapping, log, max, maxFactor);

		AntiAlignments aa;

		long start = System.nanoTime();

		try {
			aa = calculator2.getAntiAlignments(initialMarking, finalMarking);
		} catch (LPMatrixException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			aa = null;
		}

		long mid = System.nanoTime();
		System.out.println("Anti alignment computation: " + (mid - start) / 1000000.0 + " ms");

		return aa;
	}

	public AntiAlignmentValues computePrecisionAndGeneralization(AntiAlignments aa) {

		//		System.out.println("old: " + (end - mid) / 1000000.0);
		//
		//		System.out.println("old: " + Arrays.toString(aa.getMaxMinDistances()));
		//		System.out.println("new: " + Arrays.toString(aa3.getMaxMinDistances()));
		//		for (int i = 0; i < aa.getMaxMinDistances().length; i++) {
		//			if (aa.getMaxMinDistances()[i] != aa3.getMaxMinDistances()[i]) {
		//				System.out.print("OLD: ");
		//				System.out.println(toString(aa.getAntiAlignments()[i], short2label));
		//				System.out.print("NEW: ");
		//				System.out.println(toString(aa3.getAntiAlignments()[i], short2label));
		//			}
		//		}

		//		double[] weightedPrecision = computePrecision(aa, log, true, max, maxFactor);
		double[] unweightedPrecision = computePrecision(aa, log, false);

		// now for the generalization part.
		// compute the states that were visited by the aligned log, using the
		// firing sequences in the aligned log
		Map<Marking, TShortSet> statesVisitedPerTrace = getStatesVisitedPerTrace();

		int[] count = new int[log.length + 1];
		double[] recDist = new double[log.length + 1];

		int[] cAndRD = countNewStatesAndRecoveryDistance(net, initialMarking, aa.getAAFiringSequenceForLog(),
				statesVisitedPerTrace, (short) -1);
		count[log.length] = cAndRD[0];
		recDist[log.length] = cAndRD[1];

		for (short tr = 0; tr < log.length; tr++) {

			cAndRD = countNewStatesAndRecoveryDistance(net, initialMarking,
					aa.getAAFiringSequenceForLogWithoutTrace(tr), statesVisitedPerTrace, tr);
			count[tr] = cAndRD[0];
			recDist[tr] = cAndRD[1];

		}

		double[] weightedGeneralization = computeGeneralization(aa, log, true, recDist, count);
		//		double[] unweightedGeneralization = computeGeneralization(aa, log, false, recDist, count, max, maxFactor);

		//		printAntiAlignments("print", aa, alignedLog, frequencies, count, recDist, short2label, true);
		return new AntiAlignmentValues(0.5, traceFitness, unweightedPrecision[0], unweightedPrecision[1],
				weightedGeneralization[0], weightedGeneralization[1]);

	}

	public PNRepResult getPNRepResult(AntiAlignments aa, AntiAlignmentValues values) {
		Collection<SyncReplayResult> collection = new ArrayList<>();
		for (int t = 0; t <= log.length; t++) {
			collection.add(getSyncReplayResult(aa, t));
		}
		PNRepResult res = new PNRepResult(collection);

		res.addInfo(PRECISION, Double.toString(values.getPrecision()));
		res.addInfo(GENERALIZATION, Double.toString(values.getGeneralization()));
		res.addInfo(LOGPRECISION, Double.toString(values.getLogPrecision()));
		res.addInfo(LOGGENERALIZATION, Double.toString(values.getLogGeneralization()));
		res.addInfo(TRACEPRECISION, Double.toString(values.getTracePrecision()));
		res.addInfo(TRACEGENERALIZATION, Double.toString(values.getTraceGeneralization()));
		res.addInfo(PNRepResult.TRACEFITNESS, Double.toString(traceFitness));

		return res;
	}

	private double[] computeGeneralization(AntiAlignments aa, short[][] log, boolean weighted, double[] recDistances,
			int[] newStateCounts) {

		double sum = 0;
		int sumFreq = 0;
		int f = 1;
		for (int t = 0; t < log.length; t++) {

			f = weighted ? log2xLog[t].getRepresented().size() : 1;

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

		double d = Math.max(aa.getAAForLog().length, aa.getMaxDistanceForLog());
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

	private double[] computePrecision(AntiAlignments aa, short[][] log, boolean weighted) {

		double sum = 0;
		int sumFreq = 0;
		int f = 1;
		for (int t = 0; t < aa.getLogLength(); t++) {
			// Determine the number of represented xTraces by this trace in log
			f = weighted ? log2xLog[t].getRepresented().size() : 1;

			// the distance to the removed trace is at most the
			// length of the removed trace
			assert (int) (log[t].length * aa.getMaxFactorForTraces()) >= aa.getAADistanceToTrace(t);
			// we normalize on the length of the trace * maxFactor as this is
			// the given maximum
			// to the anti-alignment algorithm.
			// sum += f
			// * (aa.getAADistanceToTrace(t) / (double) ((int) (log[t].length *
			// maxFactor)));

			// We can also normalize on the length of the returned
			// anti-alignment!
			double d = Math.min(aa.getAAForLogWithoutTrace(t).length, log[t].length * aa.getMaxFactorForTraces());
			sum += d < 1 ? f : f * (aa.getAADistanceToTrace(t) / d);
			sumFreq += f;

		}
		double d = Math.min(aa.getAAForLog().length, aa.getMaxAntiAlignmentLength());
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
		for (int t_i = 0; t_i < firingSequence.size(); t_i++) {
			Transition t = firingSequence.get(t_i);
			try {
				semantics.executeExecutableTransition(t);
				length--;
			} catch (IllegalTransitionException e) {
				// so this transition was not enabled.
				assert (t.isInvisible());
				// push forward to first visible transition
				int j;
				for (j = t_i + 1; j < firingSequence.size(); j++) {
					if (firingSequence.get(j).isInvisible()) {
						firingSequence.set(j - 1, firingSequence.get(j));
					} else {
						firingSequence.set(j - 1, t);
						break;
					}
				}
				if (j == firingSequence.size()) {
					firingSequence.set(j - 1, t);
				}
				t_i--;
				continue;

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
	public Map<Marking, TShortSet> getStatesVisitedPerTrace() {

		Map<Marking, TShortSet> map = new HashMap<>();

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		for (short i = 0; i < log.length; i++) {
			AlignedRepresentative rep = (AlignedRepresentative) log2xLog[i];

			semantics.setCurrentState(initialMarking);
			if (!map.containsKey(initialMarking)) {
				map.put(initialMarking, new TShortHashSet(log.length, 0.75f, Short.MIN_VALUE));
			}
			map.get(initialMarking).add(i);
			for (Transition t : rep.getFiringSequence()) {
				try {
					semantics.executeExecutableTransition(t);
				} catch (IllegalTransitionException e) {
					e.printStackTrace();
				}
				if (!map.containsKey(semantics.getCurrentState())) {
					map.put(semantics.getCurrentState(), new TShortHashSet(log.length, 0.75f, Short.MIN_VALUE));
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

	protected SyncReplayResult getSyncReplayResult(AntiAlignments antiAlignments, int trace) {

		Vector<Transition> moves;
		if (trace >= 0 && trace < log.length) {
			moves = antiAlignments.getAAFiringSequenceForLogWithoutTrace(trace);
		} else {
			moves = antiAlignments.getAAFiringSequenceForLog();
		}

		List<StepTypes> stepTypes = new ArrayList<StepTypes>(moves.size());
		List<Object> nodeInstance = new ArrayList<Object>(moves.size());

		int e = 0;

		for (int i = 0; i < moves.size(); i++) {

			Transition t = moves.get(i);
			nodeInstance.add(t);
			if (t.isInvisible()) {
				// Move on invisible transtion
				stepTypes.add(StepTypes.MINVI);

			} else if (trace >= 0 && trace < log.length && e < log[trace].length) {
				// lookup the corresponding event in the ignored trace
				XEventClass c = short2label.get(log[trace][e]);
				e++;

				if (mapping.get(t).equals(c)) {
					// no difference between aa and trace
					stepTypes.add(StepTypes.LMGOOD);
				} else {
					// difference
					stepTypes.add(StepTypes.MREAL);
				}
			} else if (trace < 0 || trace >= log.length) {
				// check if, any trace has (at this location) an event with the correct label;
				XEventClass aaC = mapping.get(t);
				boolean found = false;
				for (int tr = 0; !found && tr < log.length; tr++) {
					if (e < log[tr].length && aaC.equals(short2label.get(log[tr][e]))) {
						stepTypes.add(StepTypes.LMGOOD);
						found = true;
					}
				}
				if (!found) {
					// difference
					stepTypes.add(StepTypes.MREAL);
				}
				e++;
			} else {
				// difference
				stepTypes.add(StepTypes.MREAL);
				e++;
			}
		}

		if (trace >= 0 && trace < log.length) {
			for (; e < log[trace].length; e++) {
				// add trailing logmoves
				XEventClass c = short2label.get(log[trace][e]);
				stepTypes.add(StepTypes.L);
				assert (c != null);
				nodeInstance.add(c);
			}
		}

		SyncReplayResult res;
		Map<String, Double> info = new HashMap<String, Double>();

		if (trace >= 0 && trace < log.length) {
			int xTrace = log2xLog[trace].getRepresented().get(0);
			res = new SyncReplayResult(nodeInstance, stepTypes, xTrace);

			info.put(HAMMINGDISTANCETOLOG, antiAlignments.getAADistanceForLogWithoutTrace(trace));
			info.put(HAMMINGDISTANCETOREMOVED, (double) antiAlignments.getAADistanceToTrace(trace));
			info.put(PNRepResult.ORIGTRACELENGTH, (double) xLog.get(xTrace).size());
			info.put(PNRepResult.TIME, antiAlignments.getTimeForTrace(trace));
		} else {
			// slightly hacking. I need to specify a trace (but I don't want to)
			res = new SyncReplayResult(nodeInstance, stepTypes, 1);
			// so remove it again
			res.getTraceIndex().remove(1);

			info.put(HAMMINGDISTANCETOLOG, antiAlignments.getAADistanceForLog());
			info.put(PNRepResult.TIME, antiAlignments.getTimeForLog());
		}

		res.setReliable(true);

		res.setInfo(info);
		return res;

	}
}
