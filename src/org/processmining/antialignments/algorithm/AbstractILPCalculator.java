package org.processmining.antialignments.algorithm;

import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import nl.tue.astar.util.LPMatrix;

import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetEdge;
import org.processmining.models.graphbased.directed.petrinet.PetrinetGraph;
import org.processmining.models.graphbased.directed.petrinet.PetrinetNode;
import org.processmining.models.graphbased.directed.petrinet.elements.Arc;
import org.processmining.models.graphbased.directed.petrinet.elements.InhibitorArc;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public abstract class AbstractILPCalculator {

	public static boolean VERBOSE = true;

	protected static final int MODE_LPSOLVE = 1;
	protected static final int MODE_GUROBI = 2;

	protected final PetrinetGraph net;
	protected final PetrinetSemantics semantics;
	protected final TObjectShortMap<String> label2short;
	protected final TShortObjectMap<String> short2label;
	protected final short transitions;
	protected final short places;
	protected final short labels;
	protected final Transition[] short2trans;
	protected final Place[] short2place;
	protected final short[] trans2label;
	protected final short[][] log;
	protected final TObjectShortHashMap<Object> trans2int;
	protected final TObjectShortHashMap<Object> place2int;
	protected final short invisibleTransitions;

	protected GRBEnv gbEnv;

	protected int mode;
	protected int useHY = 1;
	protected int cutOffLength = 5;
	protected double backtrackThreshold = 2.0;

	public AbstractILPCalculator(PetrinetGraph net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short, TShortObjectMap<String> short2label, short[][] log) {
		this.net = net;
		//		this.initialMarking = initialMarking;
		//		this.finalMarking = finalMarking;
		this.label2short = label2short;
		this.short2label = short2label;
		this.log = log;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

		mode = MODE_LPSOLVE;
		try {
			gbEnv = new GRBEnv();
			gbEnv.set(GRB.IntParam.OutputFlag, 0);
			mode = MODE_GUROBI;
		} catch (GRBException _) {
			mode = MODE_LPSOLVE;
		}

		// replay log on model (or obtain existing replay result)
		short transitions = 0;
		short places = 0;

		short2trans = new Transition[net.getTransitions().size()];
		trans2int = new TObjectShortHashMap<>(net.getTransitions().size() / 2 * 3, 0.7f, (short) 0);
		for (Transition t : net.getTransitions()) {
			if (t.isInvisible()) {
				trans2int.put(t, transitions);
				short2trans[transitions] = t;
				transitions++;
			}
		}
		invisibleTransitions = transitions;
		for (Transition t : net.getTransitions()) {
			if (!t.isInvisible()) {
				trans2int.put(t, transitions);
				short2trans[transitions] = t;
				transitions++;
			}
		}
		this.transitions = transitions;
		short2place = new Place[net.getPlaces().size()];
		place2int = new TObjectShortHashMap<>(net.getPlaces().size() / 2 * 3, 0.7f, (short) 0);
		for (Place p : net.getPlaces()) {
			place2int.put(p, places);
			short2place[places] = p;
			places++;
		}
		this.places = places;

		trans2label = new short[transitions];

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short t;
			Transition trans;
			if (e instanceof Arc) {
				if (e.getSource() instanceof Place) {
					trans = (Transition) e.getTarget();
					t = trans2int.get(trans);
				} else {
					trans = (Transition) e.getSource();
					t = trans2int.get(trans);
				}
			} else if (e instanceof InhibitorArc) {
				trans = (Transition) e.getTarget();
				t = trans2int.get(trans);
			} else {
				continue;
			}
			if (trans.isInvisible()) {
				trans2label[t] = -1;
			} else {
				trans2label[t] = label2short.get(trans.getLabel());
			}
		}

		this.labels = (short) label2short.size();
	}

	protected boolean equalLabel(short transition, short event) {
		return trans2label[transition] == event;
	}

	protected LPMatrix<?> setupMatrix(int rows, int columns) {
		if (mode == MODE_GUROBI) {
			return new LPMatrix.GUROBI(gbEnv, rows, columns);
		} else {
			return new LPMatrix.LPSOLVE(rows, columns);
		}
	}

}