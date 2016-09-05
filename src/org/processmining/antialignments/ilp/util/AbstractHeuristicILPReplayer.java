package org.processmining.antialignments.ilp.util;

import gnu.trove.iterator.TObjectByteIterator;
import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gnu.trove.map.hash.TShortObjectHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.deckfour.xes.classification.XEventClass;
import org.deckfour.xes.classification.XEventClasses;
import org.deckfour.xes.classification.XEventClassifier;
import org.deckfour.xes.info.XLogInfo;
import org.deckfour.xes.info.XLogInfoFactory;
import org.deckfour.xes.model.XEvent;
import org.deckfour.xes.model.XLog;
import org.deckfour.xes.model.XTrace;
import org.processmining.models.graphbased.directed.petrinet.PetrinetGraph;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;
import org.processmining.plugins.petrinet.replayresult.PNRepResult;
import org.processmining.plugins.petrinet.replayresult.StepTypes;
import org.processmining.plugins.replayer.replayresult.SyncReplayResult;

public abstract class AbstractHeuristicILPReplayer<P extends PetrinetGraph> {

	/**
	 * NONE OF THE FIELDS CAN BE ACCESSED WITHOUT CALLING setUpDataStructures
	 * first!
	 */
	protected Marking initialMarking;
	protected Marking finalMarking;
	protected XEventClassifier classifier;
	protected TShortObjectMap<XEventClass> short2label;
	protected TObjectShortMap<XEventClass> label2short;
	protected TransEvClassMapping mapping;
	protected XLog xLog;
	protected XEventClasses classes;
	protected Representative[] log2xLog;
	protected short[][] log;
	protected int maxTraceLength;
	protected P net;

	public AbstractHeuristicILPReplayer() {
		// Explicit constructor without parameters is essential for the use as a @PNReplayAlgorithm
	}

	public void setUpDataStructures(P net, Marking initialMarking, Marking finalMarking, XLog xLog, TransEvClassMapping mapping) {

		this.net = net;
		this.initialMarking = initialMarking;
		this.finalMarking = finalMarking;
		this.xLog = xLog;

		this.mapping = mapping;
		classifier = mapping.getEventClassifier();
		final XLogInfo summary = XLogInfoFactory.createLogInfo(xLog, classifier);
		this.classes = summary.getEventClasses();

		short2label = new TShortObjectHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);
		label2short = new TObjectShortHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);

		KeyLookupHashMap<Representative> tempLog = new KeyLookupHashMap<Representative>(xLog.size());
		short c = 0;
		maxTraceLength = 0;
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
			if (list.size() > maxTraceLength) {
				maxTraceLength = list.size();
			}
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

		log = new short[tempLog.size()][];
		log2xLog = new Representative[tempLog.size()];
		TObjectByteIterator<Representative> it = tempLog.iterator();
		while (it.hasNext()) {
			it.advance();
			log[it.key().getNumber()] = it.key().getTrace().toArray();
			log2xLog[it.key().getNumber()] = it.key();
		}

		assert allTracesUnique(log);

	}

	public void setUpDataStructures(P net, Marking initialMarking, Marking finalMarking, XLog xLog, PNRepResult alignments, TransEvClassMapping mapping) {

		this.net = net;
		this.initialMarking = initialMarking;
		this.finalMarking = finalMarking;
		this.xLog = xLog;

		this.mapping = mapping;
		classifier = mapping.getEventClassifier();
		final XLogInfo summary = XLogInfoFactory.createLogInfo(xLog, classifier);
		this.classes = summary.getEventClasses();

		short2label = new TShortObjectHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);
		label2short = new TObjectShortHashMap<>(net.getTransitions().size(), 0.7f, (short) -1);

		KeyLookupHashMap<AlignedRepresentative> tempLog = new KeyLookupHashMap<AlignedRepresentative>(xLog.size());
		short c = 0;
		maxTraceLength = 0;
		int number = 0;
		for (SyncReplayResult res : alignments) {

			List<Transition> firingSeq = new ArrayList<Transition>();
			TShortList list = new TShortArrayList(res.getStepTypes().size());

			for (int s = 0; s < res.getStepTypes().size(); s++) {
				Transition trans = null;
				XEventClass clazz = null;
				short m = label2short.getNoEntryValue();
				if (res.getStepTypes().get(s) == StepTypes.LMGOOD) {
					// synchronous move
					// Corresponding nodeStep is a transition.
					trans = (Transition) res.getNodeInstance().get(s);
					clazz = mapping.get(trans);
					m = label2short.get(clazz);
					if (m == label2short.getNoEntryValue()) {
						label2short.put(clazz, c);
						short2label.put(c, clazz);
						m = c;
						c++;
					}
				} else if (res.getStepTypes().get(s) == StepTypes.L) {
					// log move
					// Corresponding nodeStep is an event class.
					clazz = (XEventClass) res.getNodeInstance().get(s);
					if (label2short.get(clazz) == label2short.getNoEntryValue()) {
						label2short.put(clazz, c);
						short2label.put(c, clazz);
						c++;
					}
				} else if (res.getStepTypes().get(s) == StepTypes.MREAL) {
					// Model move on visible transition
					// Corresponding nodeStep is a transition.
					trans = (Transition) res.getNodeInstance().get(s);
					m = label2short.get(mapping.get(trans));
					if (m == label2short.getNoEntryValue()) {
						label2short.put(mapping.get(trans), c);
						short2label.put(c, mapping.get(trans));
						m = c;
						c++;
					}
				} else if (res.getStepTypes().get(s) == StepTypes.MINVI) {
					// model move on invisible transition
					// Corresponding nodeStep is a transition.
					trans = (Transition) res.getNodeInstance().get(s);
					if (label2short.get(mapping.get(trans)) == label2short.getNoEntryValue()) {
						label2short.put(mapping.get(trans), c);
						short2label.put(c, mapping.get(trans));
						c++;
					}
				} else {
					System.out.println("error");
				}

				if (trans != null) {
					// Append firing sequence with transition
					firingSeq.add(trans);
					// Only if label is visible, append trace
					if (m != label2short.getNoEntryValue()) {
						list.add(m);
					}
				}

			}
			AlignedRepresentative rep = new AlignedRepresentative(list, number, firingSeq);
			AlignedRepresentative existing = tempLog.getKeyIfPresent(rep);
			if (null == existing) {
				tempLog.put(rep, (byte) 1);
				number++;
			} else {
				rep = existing;
			}
			rep.addRepresentedTrace(res.getTraceIndex());
			if (list.size() > maxTraceLength) {
				maxTraceLength = list.size();
			}
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

		log = new short[tempLog.size()][];
		log2xLog = new AlignedRepresentative[tempLog.size()];
		TObjectByteIterator<AlignedRepresentative> it = tempLog.iterator();
		while (it.hasNext()) {
			it.advance();
			log[it.key().getNumber()] = it.key().getTrace().toArray();
			log2xLog[it.key().getNumber()] = it.key();
		}

		assert allTracesUnique(log);
	}

	private static boolean allTracesUnique(short[][] alignedLog) {
		for (int i = 0; i < alignedLog.length; i++) {
			for (int j = i + 1; j < alignedLog.length; j++) {
				if (Arrays.equals(alignedLog[i], alignedLog[j])) {
					return false;
				}

			}
		}
		return true;
	}

	//	private static int rebuildLogFromAlignments(PNRepResult alignments, TransEvClassMapping mapping,
	//			TObjectIntMap<AlignedRepresentative> alignedTraces, TObjectShortMap<XEventClass> label2short,
	//			TShortObjectMap<XEventClass> short2label) {
	//		// move through alignments and rebuild event log to an aligned version with short labels.
	//		// We need to build both an aligned log and aligned firing sequences. The latter includes
	//		// invisible transitions which are required to compute the visited states later.
	//
	//		// What is needed here is a consistent mapping from labels in the log to shorts and 
	//		// labels in the model to shorts. This mapping has been built internally by the replayer
	//		// but is currently not accessible...
	//		// Ugly, but map is from both XEventClass.toString() and Transition.getLabel() to short
	//
	//		// add all mapped objects first
	//		short mapped = 1;
	//		for (XEventClass ec : mapping.values()) {
	//			label2short.put(ec, mapped);
	//			short2label.put(mapped, ec);
	//			mapped++;
	//		}
	//
	//		int max = 0;
	//		int i = 0;
	//		for (SyncReplayResult alignment : alignments) {
	//			//			frequencies[i] = alignment.getTraceIndex().size();
	//
	//			ArrayList<Transition> firingSeq = new ArrayList<Transition>();
	//			//			firingSequences[i] = new ArrayList<Transition>();
	//			TShortList modelSeq = new TShortArrayList();
	//			for (int s = 0; s < alignment.getStepTypes().size(); s++) {
	//				Transition t = null;
	//				XEventClass c = null;
	//				short m = label2short.getNoEntryValue();
	//				if (alignment.getStepTypes().get(s) == StepTypes.LMGOOD) {
	//					// synchronous move
	//					// Corresponding nodeStep is a transition.
	//					t = (Transition) alignment.getNodeInstance().get(s);
	//					c = mapping.get(t);
	//					m = label2short.get(mapping.get(t));
	//				} else if (alignment.getStepTypes().get(s) == StepTypes.L) {
	//					// log move
	//					// Corresponding nodeStep is an event class.
	//					c = (XEventClass) alignment.getNodeInstance().get(s);
	//					if (label2short.get(c) == label2short.getNoEntryValue()) {
	//						label2short.put(c, mapped);
	//						short2label.put(mapped, c);
	//						mapped++;
	//					}
	//				} else if (alignment.getStepTypes().get(s) == StepTypes.MREAL) {
	//					// Model move on visible transition
	//					// Corresponding nodeStep is a transition.
	//					t = (Transition) alignment.getNodeInstance().get(s);
	//					m = label2short.get(mapping.get(t));
	//					if (m == label2short.getNoEntryValue()) {
	//						label2short.put(mapping.get(t), mapped);
	//						short2label.put(mapped, mapping.get(t));
	//						m = mapped;
	//						mapped++;
	//					}
	//				} else if (alignment.getStepTypes().get(s) == StepTypes.MINVI) {
	//					// model move on invisible transition
	//					// Corresponding nodeStep is a transition.
	//					t = (Transition) alignment.getNodeInstance().get(s);
	//					if (label2short.get(mapping.get(t)) == label2short.getNoEntryValue()) {
	//						label2short.put(mapping.get(t), mapped);
	//						short2label.put(mapped, mapping.get(t));
	//						mapped++;
	//					}
	//				} else {
	//					System.out.println("error");
	//				}
	//				if (t != null) {
	//					firingSeq.add(t);
	//
	//					if (m != label2short.getNoEntryValue()) {
	//						modelSeq.add(m);
	//					}
	//				}
	//			}
	//
	//			AlignedRepresentative at = new AlignedRepresentative(modelSeq, firingSeq);
	//			alignedTraces.adjustOrPutValue(at, alignment.getTraceIndex().size(), alignment.getTraceIndex().size());
	//			//			alignedLog[i] = modelSeq.toArray();
	//			if (modelSeq.size() > max) {
	//				max = modelSeq.size();
	//			}
	//			i++;
	//		}
	//		// all elements are mapped to shorts and we have constructed the aligned log, both in terms of
	//		// only mapped sync move elements and in terms of firing sequences.
	//		return max;
	//	}
	//
}
