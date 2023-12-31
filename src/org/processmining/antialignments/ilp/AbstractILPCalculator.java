package org.processmining.antialignments.ilp;

import org.deckfour.xes.classification.XEventClass;
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
import org.processmining.plugins.connectionfactories.logpetrinet.TransEvClassMapping;

import gnu.trove.iterator.TIntShortIterator;
import gnu.trove.map.TIntShortMap;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TIntShortHashMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import nl.tue.astar.util.ilp.LPMatrix;
import nl.tue.astar.util.ilp.LPMatrixException;

public abstract class AbstractILPCalculator {

	public static boolean VERBOSE = true;
	public static boolean NAMES = true;

	protected static final int MODE_LPSOLVE = 1;
	protected static final int MODE_GUROBI = 2;

	protected final PetrinetGraph net;
	protected final PetrinetSemantics semantics;
	protected final TObjectShortMap<XEventClass> label2short;
	protected final TShortObjectMap<XEventClass> short2label;
	protected final short transitions;
	protected final short places;
	protected final short labels;
	protected final Transition[] short2trans;
	protected final Place[] short2place;
	protected final short[] trans2label;
	protected final short[][] log;
	protected final TObjectShortHashMap<Object> trans2short;
	protected final TObjectShortHashMap<Object> place2int;
	protected final short invisibleTransitions;

	protected static GRBEnv gbEnv;

	protected int mode;
	protected int cutOffLength = 5;
	protected double backtrackThreshold = 2.0;
	protected int maxBackTrackDepth = 1;

	// Represents the matrix of consumption of the original net.
	protected final Matrix matrixAMin;

	// Represents the incidence matrix of the original net.
	protected final Matrix matrixA;
	protected TransEvClassMapping mapping;

	protected static class Matrix {
		public final short[] rows;
		public final short[] columns;
		public final short[] values;

		public Matrix(int size) {
			rows = new short[size];
			columns = new short[size];
			values = new short[size];
		}

		public int size() {
			return rows.length;
		}

		public void copyIntoMatrix(LPMatrix<?> lp, int lpRow, int lpCol) {
			for (int i = 0; i < rows.length; i++) {
				lp.setMat(lpRow + rows[i], lpCol + columns[i], values[i]);
			}
		}

		public void copyIntoMatrixFromColumn(int minCol, LPMatrix<?> lp, int lpRow, int lpCol) {
			for (int i = 0; i < rows.length; i++) {
				if (columns[i] >= minCol) {
					lp.setMat(lpRow + rows[i], lpCol + columns[i] - minCol, values[i]);
				}
			}
		}

		public void copyColumnIntoMatrix(int colM, LPMatrix<?> lp, int fromRow, int lpCol) {
			for (int i = 0; i < rows.length; i++) {
				if (columns[i] == colM) {
					lp.setMat(fromRow + rows[i], lpCol, values[i]);
				}
			}
		}

	}

	public AbstractILPCalculator(PetrinetGraph net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<XEventClass> label2short, TShortObjectMap<XEventClass> short2label,
			TransEvClassMapping mapping, short[][] log) {
		this.net = net;
		//		this.initialMarking = initialMarking;
		//		this.finalMarking = finalMarking;
		this.label2short = label2short;
		this.short2label = short2label;
		this.mapping = mapping;
		this.log = log;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

		mode = MODE_LPSOLVE;
		setGurobi();

		// replay log on model (or obtain existing replay result)
		short transitions = 0;
		short places = 0;

		short2trans = new Transition[net.getTransitions().size()];
		trans2short = new TObjectShortHashMap<>(net.getTransitions().size() / 2 * 3, 0.7f, (short) 0);
		for (Transition t : net.getTransitions()) {
			if (t.isInvisible()) {
				trans2short.put(t, transitions);
				short2trans[transitions] = t;
				transitions++;
			}
		}
		invisibleTransitions = transitions;
		for (Transition t : net.getTransitions()) {
			if (!t.isInvisible()) {
				trans2short.put(t, transitions);
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
					t = trans2short.get(trans);
				} else {
					trans = (Transition) e.getSource();
					t = trans2short.get(trans);
				}
			} else if (e instanceof InhibitorArc) {
				trans = (Transition) e.getTarget();
				t = trans2short.get(trans);
			} else {
				continue;
			}
			if (trans.isInvisible()) {
				trans2label[t] = -1;
			} else {
				trans2label[t] = label2short.get(mapping.get(trans));
			}
		}

		this.labels = (short) label2short.size();

		TIntShortMap aMinusList = new TIntShortHashMap(net.getEdges().size(), 0.7f, (Short.MIN_VALUE << 16)
				| Short.MIN_VALUE, (short) 0);
		TIntShortMap aMatrixList = new TIntShortHashMap(net.getEdges().size(), 0.7f, (Short.MIN_VALUE << 16)
				| Short.MIN_VALUE, (short) 0);

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			int p, t;
			short dir;
			if (e instanceof Arc) {
				if (e.getSource() instanceof Place) {
					p = place2int.get(e.getSource());
					t = trans2short.get(e.getTarget());
					dir = (short) -((Arc) e).getWeight();
				} else {
					t = trans2short.get(e.getSource());
					p = place2int.get(e.getTarget());
					dir = (short) ((Arc) e).getWeight();
				}
				int m = (p << 16) | t;
				if ((dir < 0)  || (trans2label[t] < 0)) {//{//BVD: Don't treat invisible steps in a special way
					aMinusList.adjustOrPutValue(m, dir, dir);
				}
				aMatrixList.adjustOrPutValue(m, dir, dir);
			}
		}

		matrixAMin = new Matrix(aMinusList.size());
		int i = 0;
		TIntShortIterator it = aMinusList.iterator();
		while (it.hasNext()) {
			it.advance();
			short row = (short) (it.key() >>> 16);
			short col = (short) (it.key() & 0x0000FFFF);
			matrixAMin.rows[i] = row;
			matrixAMin.columns[i] = col;
			matrixAMin.values[i] = it.value();
			i++;
		}

		matrixA = new Matrix(aMatrixList.size());
		i = 0;
		it = aMatrixList.iterator();
		while (it.hasNext()) {
			it.advance();
			short row = (short) (it.key() >>> 16);
			short col = (short) (it.key() & 0x0000FFFF);
			matrixA.rows[i] = row;
			matrixA.columns[i] = col;
			matrixA.values[i] = it.value();
			i++;
		}

	}

	protected boolean equalLabel(short transition, short event) {
		return trans2label[transition] == event;
	}

	protected LPMatrix<?> setupMatrix(int rows, int columns) throws LPMatrixException {
		if (mode == MODE_GUROBI) {
			return new LPMatrix.SPARSE.GUROBI(gbEnv, rows, columns);
		} else {
			return new LPMatrix.SPARSE.LPSOLVE(rows, columns);
		}
	}

	public void setCutOffLength(int cutOffLength) {
		this.cutOffLength = cutOffLength;
	}

	public boolean setLPSolve() {
		try {
			mode = MODE_LPSOLVE;
			// try if LpSolve works by creating a simple LP.
			LpSolve.makeLp(1, 1).deleteLp();
			return true;
		} catch (LpSolveException | UnsatisfiedLinkError | NoClassDefFoundError e) {
			return false;
		}
	}

	public boolean setGurobi() {
		try {
			if (gbEnv == null) {
				gbEnv = new GRBEnv();
				gbEnv.set(GRB.IntParam.OutputFlag, 0);
				gbEnv.set(GRB.IntParam.Threads, 1);
			}
			mode = MODE_GUROBI;
			return true;
		} catch (GRBException | UnsatisfiedLinkError | NoClassDefFoundError e) {
			System.err.println(e.getMessage());
			mode = MODE_LPSOLVE;
			return false;
		}
	}

	protected void unfire(Transition last, Marking marking) {
		for (PetrinetEdge<?, ?> e : net.getInEdges(last)) {
			if (e instanceof Arc) {
				Arc a = (Arc) e;
				Place p = (Place) a.getSource();
				int w = a.getWeight();
				marking.add(p, w);
			}
		}
		for (PetrinetEdge<?, ?> e : net.getOutEdges(last)) {
			if (e instanceof Arc) {
				Arc a = (Arc) e;
				Place p = (Place) a.getTarget();
				int w = a.getWeight();
				while (w > 0) {
					marking.remove(p);
					w--;
				}
			}
		}
	}

	public void setBacktrackLimit(int backtrackLimit) {
		this.maxBackTrackDepth = backtrackLimit;
	}

	public void setBacktrackThreshold(double backtrackThreshold) {
		this.backtrackThreshold = backtrackThreshold;
	}

	//	public void disposeGurobi() {
	//		finalize();
	//	}
	//
	//	protected void finalize() {
	//		if (gbEnv != null) {
	//			try {
	//				gbEnv.dispose();
	//			} catch (GRBException e) {
	//				// Nothing to do here. Forget about it.
	//			}
	//		}
	//
	//	}
}