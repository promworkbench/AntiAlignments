package org.processmining.antialignments.algorithm;

import gnu.trove.list.TShortList;
import gnu.trove.list.array.TShortArrayList;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.map.hash.TObjectShortHashMap;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBModel;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import lpsolve.LpSolve;
import lpsolve.LpSolveException;

import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetEdge;
import org.processmining.models.graphbased.directed.petrinet.PetrinetNode;
import org.processmining.models.graphbased.directed.petrinet.elements.Arc;
import org.processmining.models.graphbased.directed.petrinet.elements.InhibitorArc;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class AntiAlignmentILPCalculator {

	private final Petrinet net;
	//	private final Marking initialMarking;
	private final PetrinetSemantics semantics;
	//	private final Marking finalMarking;
	private final TObjectShortMap<String> label2short;
	private TShortObjectMap<String> short2label;
	private short transitions;
	private short places;
	private Transition[] short2trans;
	private Place[] short2place;

	private short[] trans2label;
	private final short[][] log;
	private final int maxFactor;

	private final GRBEnv gbEnv;
	private TObjectShortHashMap<Object> trans2int;
	private TObjectShortHashMap<Object> place2int;
	private final int maxLength;

	public AntiAlignmentILPCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short, TShortObjectMap<String> short2label, short[][] log, int maxLength,
			int maxFactor) throws LpSolveException, GRBException {
		this.net = net;
		//		this.initialMarking = initialMarking;
		//		this.finalMarking = finalMarking;
		this.label2short = label2short;
		this.short2label = short2label;
		this.log = log;
		this.maxLength = maxLength;
		this.maxFactor = maxFactor;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

		setupDataStructures();

		gbEnv = new GRBEnv();
		gbEnv.set(GRB.IntParam.OutputFlag, 0);

		List<Transition> firingSequence = new ArrayList<Transition>(maxLength * maxFactor);
		solveByDrillingDown(maxLength * maxFactor, initialMarking, finalMarking, firingSequence, 0);
		System.out.println(firingSequence);

	}

	protected void solveByDrillingDown(int maxLength, Marking initialMarking, Marking finalMarking,
			List<Transition> firingSequence, int startTracesAt) throws LpSolveException, GRBException {

		System.out
				.println("Solving maxlength " + maxLength + " from marking " + initialMarking + " to " + finalMarking);
		if (maxLength > 10) {

			int maxLengthX = maxLength / 2 + maxLength % 2;
			int maxLengthY = maxLength - maxLengthX;

			LPMatrix matrix = setupLpForSplit(maxLengthX, maxLengthY, true, initialMarking, finalMarking, startTracesAt);
			double[] vars = new double[matrix.getNcolumns()];

			Marking intermediate = determineSplitMarking(matrix, vars);

			System.out.println("Intermediate marking: " + intermediate);
			if (!intermediate.equals(initialMarking)) {
				solveByDrillingDown(maxLength / 2 + maxLength % 2, initialMarking, intermediate, firingSequence,
						startTracesAt);
				if (!intermediate.equals(finalMarking)) {
					solveByDrillingDown(maxLength / 2, intermediate, finalMarking, firingSequence, startTracesAt
							+ maxLengthX);
				}
			} else {
				// Initial Marking is kept intact. Possible spurious trace, or transition invariant detected
				// Now do 1+ rest search.
				matrix = setupLpForSplit(1, maxLength - 1, true, initialMarking, finalMarking, startTracesAt);
				intermediate = determineSplitMarking(matrix, vars);
				for (int t = 0; t < transitions; t++) {
					if (vars[t] > 0.5) {
						assert ((int) (vars[t] + 0.5)) == 1 || short2trans[t].isInvisible();
						do {
							firingSequence.add(short2trans[t % transitions]);
							vars[t] = vars[t] - 1;
						} while (vars[t] > 0.5);
					}
				}

				if (!intermediate.equals(finalMarking)) {
					solveByDrillingDown(maxLength - 1, intermediate, finalMarking, firingSequence, startTracesAt + 1);
				}
			}

		} else {
			LPMatrix matrix = setupLpForFullSequence(maxLength, true, initialMarking, finalMarking, startTracesAt);
			//			matrix.printLp();
			double[] vars = new double[matrix.getNcolumns()];
			solveForFullSequence(matrix, -1, maxLength, null, vars);
			for (int t = 0; t < vars.length - 2; t++) {
				if (vars[t] > 0.5) {
					assert ((int) (vars[t] + 0.5)) == 1;
					firingSequence.add(short2trans[t % transitions]);
				}
			}
			//			System.out.println(Arrays.toString(vars));
		}

	}

	protected Marking determineSplitMarking(LPMatrix matrix, double[] vars) throws LpSolveException, GRBException {
		//			matrix.toLpSolve().printLp();

		//			matrix.printLp();

		solveForSplit(matrix, null, vars);

		// get the intermediate marking:
		Marking intermediate = new Marking();
		// Iterate over the X vector
		for (int p = 0; p < places; p++) {
			double effect = -matrix.getRh(p + places);
			for (int t = 0; t < transitions; t++) {
				effect += vars[t] * matrix.getMat(p, t);
			}
			assert effect >= 0;
			if (effect > 0) {
				intermediate.add(short2place[p], (int) (effect + 0.5));
			}
		}
		return intermediate;
	}

	public AntiAlignments getAntiAlignments(Marking initialMarking, Marking finalMarking) throws LpSolveException,
			GRBException {

		System.out.println("Setting up LP");
		LPMatrix lpMatrix = setupLpForFullSequence(maxLength * maxFactor, true, initialMarking, finalMarking, 0);

		System.out.println("IP with " + lpMatrix.getNcolumns() + " columns and " + lpMatrix.getNrows() + " rows.");

		final AntiAlignments antiAlignments = new AntiAlignments(log.length);

		System.out.println("-----------------------------------------------");
		System.out.println("Whole log. ");
		System.out.println("Maxlength: " + maxFactor * maxLength);
		System.out.flush();

		//		lp.printLp();
		double[] result = new double[lpMatrix.getNcolumns()];
		solveForFullSequence(lpMatrix, -1, maxFactor * maxLength, null, result);// "D:/temp/antialignment/ilp_instance_log.mps");

		antiAlignments.getMaxMinDistances()[log.length] = (int) (result[result.length - 2] + 0.5);
		antiAlignments.getAntiAlignments()[log.length] = getAntiAlignment(result);
		antiAlignments.getTraces()[log.length] = getFiringSequence(result);
		antiAlignments.getMaxDistances()[log.length] = Math.max(maxLength * maxFactor,
				antiAlignments.getAntiAlignments()[log.length].length);

		//		int[] basis = new int[lp.getNrows() + lp.getNcolumns() + 1];
		//		lp.getBasis(basis, true);

		for (int t = 0; t < log.length; t++) {
			System.out.println("-----------------------------------------------");
			System.out.println("Removed Trace: " + toString(log[t]));
			System.out.println("Maxlength: " + maxFactor * log[t].length);
			System.out.flush();
			//			lp.setBasis(basis, true);
			solveForFullSequence(lpMatrix, t, maxFactor * log[t].length, null, result);//, "D:/temp/antialignment/ilp_instance_t" + t + ".mps");

			antiAlignments.getMaxMinDistances()[t] = (int) (result[result.length - 2] + 0.5);
			antiAlignments.getAntiAlignments()[t] = getAntiAlignment(result);
			antiAlignments.getMaxDistances()[t] = Math.max(log[t].length, antiAlignments.getAntiAlignments()[t].length);
			antiAlignments.getTraces()[t] = getFiringSequence(result);
			antiAlignments.getTraceDistances()[t] = (int) (result[result.length - 1] + 0.5);

		}
		System.out.println("-----------------------------------------------");

		//		System.out.println(Arrays.toString(vars));

		return antiAlignments;

	}

	private Vector<?> getFiringSequence(double[] result) {
		// TODO Auto-generated method stub
		return null;
	}

	private short[] getAntiAlignment(double[] result) {
		TShortList list = new TShortArrayList();
		for (int i = 0; i < result.length - 2; i++) {
			if (result[i] > 0.5) {
				// fired transition.
				short trans = (short) (i % trans2label.length);
				if (trans2label[trans] >= 0) {
					list.add(trans2label[trans]);
				}
			}
		}
		return list.toArray();
	}

	protected void solveForFullSequence(LPMatrix lpMatrix, int traceToIgnore, int maxLength, String filename,
			double[] vars) throws LpSolveException, GRBException {

		int row;

		row = lpMatrix.getNrows() - log.length + traceToIgnore - 1;

		int oldRH = (int) (lpMatrix.getRh(row) + 0.5);
		if (traceToIgnore >= 0) {
			// h doesn't matter
			lpMatrix.setMat(row, lpMatrix.getNcolumns() - 2, 0);
			// minus g for the removed trace.
			lpMatrix.setMat(row, lpMatrix.getNcolumns() - 1, -1);
			lpMatrix.setConstrType(row, LpSolve.EQ);
		}
		lpMatrix.setRh(lpMatrix.getNrows() - 1, maxLength);

		//		lp.printLp();
		//		System.out.print("Pushing to Solver .... ");

		//		LpSolve lp = lpMatrix.toLpSolve();
		GRBModel grbModel = lpMatrix.toGurobi(gbEnv);

		//		System.out.println("Done");
		try {
			int result = -1;
			if (filename != null) {
				//				lp.writeMps(filename);
			} else {

				//				lp.printLp();

				//				result = lp.solve();

				grbModel.optimize();

				//				if (result == lp.INFEASIBLE) {
				//					System.out.println("Result: INFEASIBLE");
				//				} else if (result == lpMatrix.OPTIMAL) {
				//					System.out.println("Result: OPTIMAL");
				//				} else {
				//					System.out.println("Result: " + result);
				//				}
				//
				//				lp.getVariables(vars);

				if (grbModel.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {

					for (int j = 0; j < vars.length; j++)
						vars[j] = grbModel.getVars()[j].get(GRB.DoubleAttr.X);

				} else if (grbModel.get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE) {

					lpMatrix.toLpSolve().writeLp("D:/temp/antialignment/test.lp");
					lpMatrix.toLpSolve().writeMps("D:/temp/antialignment/test.mps");
					System.out.println("INFEASIBLE");
				}

				for (int i = 0; i < vars.length - 2; i++) {
					if (vars[i] > 0.5) {
						System.out.print(lpMatrix.getColName(i));
						System.out.print("|");
					}
				}
				System.out.print("Dlog: ");
				System.out.print(((int) (vars[vars.length - 2] + .5)));

				System.out.print(" Drem: ");
				System.out.print(((int) (vars[vars.length - 1] + .5)));

				if (traceToIgnore >= 0) {
					// minimize h
					lpMatrix.setMat(row, lpMatrix.getNcolumns() - 2, -1);
					// ignore g.
					lpMatrix.setMat(row, lpMatrix.getNcolumns() - 1, 0);
					lpMatrix.setConstrType(row, LpSolve.GE);
				}
				//		lp.printLp();
				System.out.println();
			}
		} finally {
			grbModel.dispose();
			//			lp.deleteAndRemoveLp();
		}

	}

	protected void setupDataStructures() {
		// replay log on model (or obtain existing replay result)
		transitions = 0;
		places = 0;

		int[][] costs = new int[label2short.size() + 1][maxLength];
		for (int t = 0; t < log.length; t++) {
			for (int pos = 0; pos < log[t].length && pos < maxLength; pos++) {
				short label = log[t][pos];
				costs[label][pos]++;
			}
		}

		short2trans = new Transition[net.getTransitions().size()];
		trans2int = new TObjectShortHashMap<>(net.getTransitions().size() / 2 * 3, 0.7f, (short) 0);
		for (Transition t : net.getTransitions()) {
			trans2int.put(t, transitions);
			short2trans[transitions] = t;
			transitions++;
		}
		short2place = new Place[net.getPlaces().size()];
		place2int = new TObjectShortHashMap<>(net.getPlaces().size() / 2 * 3, 0.7f, (short) 0);
		for (Place p : net.getPlaces()) {
			place2int.put(p, places);
			short2place[places] = p;
			places++;
		}

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
	}

	protected LPMatrix setupLpForFullSequence(int maxLength, boolean integerVariables, Marking initialMarking,
			Marking finalMarking, int startTracesAt) throws LpSolveException {

		//		LpSolve lp = LpSolve.makeLp((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 2);
		LPMatrix lp = new LPMatrix((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 2);

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short p, t;
			int dir;
			int type = LpSolve.GE;
			Transition trans;
			if (e instanceof Arc) {
				if (e.getSource() instanceof Place) {
					p = place2int.get(e.getSource());
					trans = (Transition) e.getTarget();
					t = trans2int.get(trans);
					dir = -((Arc) e).getWeight();
				} else {
					trans = (Transition) e.getSource();
					t = trans2int.get(trans);
					p = place2int.get(e.getTarget());
					dir = ((Arc) e).getWeight();
				}
			} else if (e instanceof InhibitorArc) {
				p = place2int.get(e.getSource());
				trans = (Transition) e.getTarget();
				t = trans2int.get(trans);
				dir = 0;
				type = LpSolve.EQ;
			} else {
				continue;
			}
			for (int i = p; i < (maxLength + 1) * places; i += places) {
				for (int idx = t; idx < (i / places + 1) * transitions && idx < maxLength * transitions; idx += transitions) {
					// set up all the A matrixes.
					if (idx / transitions < i / places || trans.isInvisible()) {
						// This is one of the matrixes before the last.
						lp.setMat(i, idx, lp.getMat(i, idx) + dir);
					} else if (dir < 0 || i >= maxLength * places) {
						// This is the last matrix, only count the negative values
						lp.setMat(i, idx, lp.getMat(i, idx) + dir);
					}
					lp.setInt(idx, integerVariables);
					lp.setUpbo(idx, 1.0);
					// target function
					//						if (label2short.get(trans.getLabel()) >= 0) {
					//							lp.setMat(0, idx, maxLength
					//									* costs[label2short.get(trans.getLabel())][(idx - 1) / transitions] - 1);
					//						}
					lp.setColName(idx, trans.getLabel().replace("\\n$invisible$", "") + "_" + idx / transitions);
				}
				lp.setRowName(i, "A" + i / places + "_" + (i % places));
				if (i < maxLength * places) {
					// at least -current marking
					// or equal to - current marking in case of an inhibitorArc
					lp.setConstrType(i, type);
				} else {
					// equal to final marking - output of xi?
					lp.setConstrType(i, LpSolve.EQ);
				}
			}

		}
		int row = places * (maxLength + 1);
		// set up the comparisons and sums
		for (int i = 0; i < maxLength - 1; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			for (int t = (i + 1) * transitions; t < (i + 2) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : -1);
			}
			lp.setRowName(row + i, "X" + i + "X" + (i + 1));
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// sos constraints, sum max 1
		row = places * (maxLength + 1) + maxLength - 1;
		for (int i = 0; i < maxLength; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t, trans2label[t % transitions] < 0 ? 0 : 1);
			}
			lp.setRowName(row + i, "X" + i + ".1");
			lp.setConstrType(row + i, LpSolve.LE);
		}

		row = places * (maxLength + 1) + 2 * maxLength - 1;
		for (int t = 0; t < log.length; t++) {
			// what's the hamming distance of the x's to trace t?
			for (int e = startTracesAt; e < maxLength; e++) {
				for (int tr = 0; tr < trans2label.length; tr++) {
					if (e >= log[t].length && trans2label[tr] > 0) {
						lp.setMat(row + t, e * transitions + tr, 1);
					} //					if ((e >= log[t].length || trans2label[tr] != log[t][e]) && trans2label[tr] > 0) {
					else if (e < log[t].length && trans2label[tr] != log[t][e] && trans2label[tr] >= 0) {
						//						lp.setMat(row + t, e * transitions + tr + 1, 1);
						lp.setMat(row + t, e * transitions + tr, 0);
					} else if (trans2label[tr] < 0) {
						// tau steps do not count
						lp.setMat(row + t, e * transitions + tr, 0);
					} else {
						// based on the index and the fact that this transition matches that label,
						// set the value of this trace's row to 0;
						//						lp.setMat(row + t, e * transitions + tr + 1, 0);
						lp.setMat(row + t, e * transitions + tr, -1);
					}
				}
			}
			lp.setRowName(row + t, "trace_" + t);
			// minus h
			lp.setMat(row + t, transitions * maxLength, -1);
			// minus g for the removed trace.
			lp.setMat(row + t, transitions * maxLength + 1, 0);
			lp.setConstrType(row + t, LpSolve.GE);
		}

		// Setup maxLength constraint
		row = (maxLength + 1) * places + 2 * maxLength + log.length - 1;
		for (int t = 0; t < transitions * maxLength; t++) {
			lp.setMat(row, t, trans2label[t % transitions] < 0 ? 0 : 1);
			lp.setObjective(t, -1);
		}
		lp.setRowName(row, "SumX");
		lp.setConstrType(row, LpSolve.LE);

		// maximize h
		lp.setInt(transitions * maxLength, integerVariables);
		lp.setObjective(transitions * maxLength, maxLength * maxLength);
		// minimuze g
		lp.setInt(transitions * maxLength + 1, integerVariables);
		lp.setObjective(transitions * maxLength + 1, -maxLength);

		lp.setMaxim();

		double[] rhs = new double[(maxLength + 1) * places + 2 * maxLength + log.length];

		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row < places * (maxLength + 1); row += places) {
				rhs[row] = -initialMarking.occurrences(p);
			}
		}
		for (Place p : finalMarking.baseSet()) {
			row = maxLength * places + place2int.get(p);
			rhs[row] += finalMarking.occurrences(p);
		}
		for (row = (maxLength + 1) * places + maxLength; row < (maxLength + 1) * places + 2 * maxLength; row++) {
			rhs[row - 1] = 1;
		}

		row = lp.getNrows() - log.length - 1;
		for (int t = 0; t < log.length; t++) {
			rhs[row + t] = -log[t].length;
		}

		row = (maxLength + 1) * places + 2 * maxLength + log.length - 1;
		rhs[row] = maxLength;

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
	}

	private double getCost(Transition trans, int i) {
		return 1;
	}

	private String toString(short[] a) {
		if (a == null)
			return "null";
		int iMax = a.length - 1;
		if (iMax == -1)
			return "[]";

		StringBuilder b = new StringBuilder();
		b.append('[');
		for (int i = 0;; i++) {
			b.append(short2label.get(a[i]));
			if (i == iMax)
				return b.append(']').toString();
			b.append(", ");
		}
	}

	protected LPMatrix setupLpForSplit(int maxLengthX, int maxLengthY, boolean integerVariables,
			Marking initialMarking, Marking finalMarking, int startTracesAt) throws LpSolveException {

		//		LpSolve lp = LpSolve.makeLp((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 2);
		LPMatrix lp = new LPMatrix(2 * places + 2 * log.length + 2, transitions * 2 + 2);

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			short p, t;
			int dir;
			int type = LpSolve.GE;
			Transition trans;
			if (e instanceof Arc) {
				if (e.getSource() instanceof Place) {
					p = place2int.get(e.getSource());
					trans = (Transition) e.getTarget();
					t = trans2int.get(trans);
					dir = -((Arc) e).getWeight();
				} else {
					trans = (Transition) e.getSource();
					t = trans2int.get(trans);
					p = place2int.get(e.getTarget());
					dir = ((Arc) e).getWeight();
				}
			} else if (e instanceof InhibitorArc) {
				p = place2int.get(e.getSource());
				trans = (Transition) e.getTarget();
				t = trans2int.get(trans);
				dir = 0;
				type = LpSolve.EQ;
			} else {
				continue;
			}
			// set up all the two A matrixes.
			lp.setMat(p, t, lp.getMat(p, t) + dir);
			if (maxLengthX > 1 || trans.isInvisible() || dir < 0) {
				lp.setMat(p + places, t, lp.getMat(p + places, t) + dir);
			}
			lp.setMat(p, t + transitions, lp.getMat(p, t + transitions) + dir);

			lp.setInt(t, integerVariables);
			lp.setInt(t + transitions, integerVariables);

			lp.setColName(t, trans.getLabel().replace("\\n$invisible$", "") + ",X");
			lp.setColName(t + transitions, trans.getLabel().replace("\\n$invisible$", "") + ",Y");

			lp.setRowName(p, "A" + (maxLengthX == 1 ? "'1" : "1") + p);
			lp.setRowName(p + places, "A2" + p);

			lp.setConstrType(p, LpSolve.EQ);
			lp.setConstrType(p + places, LpSolve.GE);

			if (!trans.isInvisible()) {
				lp.setMat(2 * places, t, 1);
				lp.setMat(2 * places + 1, t + transitions, 1);

			} else {
				// where possible, transfer tau steps to Y instead of X
				lp.setObjective(transitions + t, 1);

			}
		}
		// The length of the X vector is equal to half the maxLength
		lp.setConstrType(2 * places, LpSolve.EQ);
		// The remainder is the Y vector
		lp.setConstrType(2 * places + 1, LpSolve.LE);

		int row = 2 * places + 2;

		for (int t = 0; t < log.length; t++) {
			// what's the hamming distance of the x's to trace t[0..maxLengthX]?
			for (int e = startTracesAt; e < startTracesAt + maxLengthX; e++) {
				for (int tr = 0; tr < trans2label.length; tr++) {
					if (e >= log[t].length && trans2label[tr] > 0) {
						lp.setMat(row + 2 * t, tr, lp.getMat(row + 2 * t, tr) + 1);
					} else if (e < log[t].length && trans2label[tr] != log[t][e] && trans2label[tr] >= 0) {

						//SKIP
						lp.setMat(row + 2 * t, tr, lp.getMat(row + 2 * t, tr) + 0);
					} else if (trans2label[tr] < 0) {
						// tau steps do not count
						//SKIP
						lp.setMat(row + 2 * t, tr, lp.getMat(row + 2 * t, tr) + 0);
					} else {
						// based on the index and the fact that this transition matches that label,
						// set the value of this trace's row to 0;
						//						lp.setMat(row + t, e * transitions + tr + 1, 0);
						lp.setMat(row + 2 * t, tr, lp.getMat(row + 2 * t, tr) - 1);
					}
				}
			}
			for (int e = startTracesAt + maxLengthX; e < maxLengthX + maxLengthY; e++) {
				for (int tr = 0; tr < trans2label.length; tr++) {
					if (e >= log[t].length && trans2label[tr] > 0) {
						lp.setMat(row + 2 * t + 1, transitions + tr, lp.getMat(row + 2 * t + 1, tr) + 1);
					} else if (e < log[t].length && trans2label[tr] != log[t][e] && trans2label[tr] >= 0) {
						//SKIP
						lp.setMat(row + 2 * t + 1, transitions + tr, lp.getMat(row + 2 * t + 1, tr) + 0);
					} else if (trans2label[tr] < 0) {
						// tau steps do not count
						//SKIP
						lp.setMat(row + 2 * t + 1, transitions + tr, lp.getMat(row + 2 * t + 1, tr) + 0);
					} else {
						// based on the index and the fact that this transition matches that label,
						// set the value of this trace's row to 0;
						//						lp.setMat(row + t, e * transitions + tr + 1, 0);
						lp.setMat(row + 2 * t + 1, transitions + tr, lp.getMat(row + 2 * t + 1, tr) - 1);
					}
				}
			}
			lp.setRowName(row + 2 * t, "t_" + t + "_X");
			lp.setRowName(row + 2 * t + 1, "t_" + t + "_Y");
			// minus h
			lp.setMat(row + 2 * t, 2 * transitions, -1);
			lp.setMat(row + 2 * t + 1, 2 * transitions + 1, -1);
			lp.setConstrType(row + 2 * t, LPMatrix.GE);
			lp.setConstrType(row + 2 * t + 1, LPMatrix.GE);
		}

		// maximize h
		lp.setColName(2 * transitions, "HX");
		lp.setObjective(2 * transitions, maxLengthX * maxLengthX);
		// maximize g
		lp.setColName(2 * transitions + 1, "HY");
		lp.setObjective(2 * transitions + 1, maxLengthY * maxLengthY);

		lp.setMaxim();

		double[] rhs = new double[2 * places + 2 * log.length + 2];

		for (Place p : initialMarking.baseSet()) {
			rhs[place2int.get(p)] = -initialMarking.occurrences(p);
			rhs[place2int.get(p) + places] = -initialMarking.occurrences(p);
		}
		for (Place p : finalMarking.baseSet()) {
			rhs[place2int.get(p)] += finalMarking.occurrences(p);
		}
		rhs[2 * places] = maxLengthX;
		rhs[2 * places + 1] = maxLengthY;

		row = 2 * places + 2;
		for (int t = 0; t < log.length; t++) {
			if (log[t].length < maxLengthX) {
				rhs[row + 2 * t] = -log[t].length;
			} else {
				rhs[row + 2 * t] = -maxLengthX;
			}
			if (log[t].length > maxLengthX) {
				rhs[row + 2 * t + 1] = -(log[t].length - maxLengthX);
			} else {
				rhs[row + 2 * t + 1] = 0;
			}

		}

		lp.setRhVec(rhs);

		//		lp.printLp();

		return lp;
	}

	protected void solveForSplit(LPMatrix lpMatrix, String filename, double[] vars) throws LpSolveException,
			GRBException {

		//		System.out.print("Pushing to Solver .... ");

		//		LpSolve lp = lpMatrix.toLpSolve();
		GRBModel grbModel = lpMatrix.toGurobi(gbEnv);

		//		System.out.println("Done");
		try {
			if (filename != null) {
				//				lp.writeMps(filename);
			} else {

				//				lp.printLp();

				//				result = lp.solve();

				grbModel.optimize();

				//				if (result == lp.INFEASIBLE) {
				//					System.out.println("Result: INFEASIBLE");
				//				} else if (result == lpMatrix.OPTIMAL) {
				//					System.out.println("Result: OPTIMAL");
				//				} else {
				//					System.out.println("Result: " + result);
				//				}
				//
				//				lp.getVariables(vars);

				if (grbModel.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {

					for (int j = 0; j < vars.length; j++)
						vars[j] = grbModel.getVars()[j].get(GRB.DoubleAttr.X);
				}
				//				System.out.println(Arrays.toString(vars));

			}
		} finally {
			grbModel.dispose();
			//			lp.deleteAndRemoveLp();
		}

	}
}
