package org.processmining.antialignments.algorithm;

import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.hash.TObjectShortHashMap;

import java.util.Arrays;

import lpsolve.LpSolve;
import lpsolve.LpSolveException;

import org.processmining.antialignments.pathfinder.AntiAlignments;
import org.processmining.models.graphbased.directed.petrinet.Petrinet;
import org.processmining.models.graphbased.directed.petrinet.PetrinetEdge;
import org.processmining.models.graphbased.directed.petrinet.PetrinetNode;
import org.processmining.models.graphbased.directed.petrinet.elements.Arc;
import org.processmining.models.graphbased.directed.petrinet.elements.Place;
import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;
import org.processmining.models.semantics.petrinet.PetrinetSemantics;
import org.processmining.models.semantics.petrinet.impl.PetrinetSemanticsFactory;

public class AntiAlignmentILPCalculator {

	private final Petrinet net;
	private final Marking initialMarking;
	private final PetrinetSemantics semantics;
	private final Marking finalMarking;
	private final TObjectShortMap<String> label2short;

	public AntiAlignmentILPCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short) {
		this.net = net;
		this.initialMarking = initialMarking;
		this.finalMarking = finalMarking;
		this.label2short = label2short;

		this.semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		this.semantics.initialize(net.getTransitions(), initialMarking);

	}

	public AntiAlignments getAntiAlignments(final short[][] log, int maxLength, int maxFactor) throws LpSolveException {

		LpSolve lp = setupLp(log, maxLength * maxFactor);

		System.out.println(label2short);

		System.out.println("-----------------------------------------------");
		System.out.println("Whole log. ");
		System.out.println("Maxlength: " + maxFactor * maxLength);
		solve(lp, log.length, -1, maxFactor * maxLength);

		int[] basis = new int[lp.getNrows() + lp.getNcolumns() + 1];
		lp.getBasis(basis, true);

		for (int t = 0; t < log.length; t++) {
			System.out.println("-----------------------------------------------");
			System.out.println("Removed Trace: " + Arrays.toString(log[t]));
			System.out.println("Maxlength: " + maxFactor * log[t].length);
			lp.setBasis(basis, true);
			solve(lp, log.length, t, maxFactor * log[t].length);
		}
		System.out.println("-----------------------------------------------");

		//		System.out.println(Arrays.toString(vars));
		lp.deleteAndRemoveLp();

		return null;

	}

	protected void solve(LpSolve lp, int logLength, int traceToIgnore, int maxLength) throws LpSolveException {

		int row = lp.getNrows() - logLength + traceToIgnore;

		int oldRH = (int) (lp.getRh(row) + 0.5);
		if (traceToIgnore >= 0) {
			lp.setRh(row, -maxLength);
		}
		lp.setRh(lp.getNrows(), maxLength);

		int result = lp.solve();
		if (result == lp.INFEASIBLE) {
			System.out.println("Result: INFEASIBLE");
		} else if (result == lp.OPTIMAL) {
			System.out.println("Result: OPTIMAL");
		} else {
			System.out.println("Result: " + result);
		}

		double[] vars = new double[lp.getNcolumns()];
		lp.getVariables(vars);

		for (int i = 0; i < vars.length - 1; i++) {
			if (vars[i] > 0.5) {
				System.out.print(lp.getColName(i + 1));
				System.out.print("|");
			}
		}
		System.out.print("HD: ");
		System.out.print(((int) vars[vars.length - 1]));
		System.out.println();

		if (traceToIgnore >= 0) {
			lp.setRh(row, oldRH);
		}

	}

	protected LpSolve setupLp(final short[][] log, int maxLength) throws LpSolveException {

		// replay log on model (or obtain existing replay result)
		short transitions = 0;
		short places = 0;

		int[][] costs = new int[label2short.size() + 1][maxLength];
		for (int t = 0; t < log.length; t++) {
			for (int pos = 0; pos < log[t].length && pos < maxLength; pos++) {
				short label = log[t][pos];
				costs[label][pos]++;
			}
		}

		TObjectShortMap<Transition> trans2int = new TObjectShortHashMap<>(net.getTransitions().size() / 2 * 3, 0.7f,
				(short) 0);
		for (Transition t : net.getTransitions()) {
			transitions++;
			trans2int.put(t, transitions);
		}
		TObjectShortMap<Place> place2int = new TObjectShortHashMap<>(net.getPlaces().size() / 2 * 3, 0.7f, (short) 0);
		for (Place p : net.getPlaces()) {
			places++;
			place2int.put(p, places);
		}

		LpSolve lp = LpSolve.makeLp((maxLength + 1) * places + 2 * maxLength + log.length, transitions * maxLength + 1);

		short[] trans2label = new short[transitions];

		for (PetrinetEdge<? extends PetrinetNode, ? extends PetrinetNode> e : net.getEdges()) {
			if (e instanceof Arc) {
				short p, t;
				int dir;
				Transition trans;
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
				double tmp = lp.getMat(p, t) + dir;
				if (trans.isInvisible()) {
					trans2label[t - 1] = 0;
				} else {
					trans2label[t - 1] = label2short.get(trans.getLabel());
				}
				for (int i = p; i <= (maxLength + 1) * places; i += places) {
					for (int idx = t; idx <= ((i - 1) / places + 1) * transitions && idx <= maxLength * transitions; idx += transitions) {
						// set up all the A matrixes.
						if ((idx - 1) / transitions < (i - 1) / places) {
							// This is one of the matrixes before the last.
							lp.setMat(i, idx, tmp);
						} else if (dir < 0 || i > maxLength * places) {
							// This is the last matrix, only count the negative values
							lp.setMat(i, idx, lp.getMat(i, idx) + dir);
						}
						lp.setInt(idx, true);
						lp.setUpbo(idx, 1.0);
						// target function
						//						if (label2short.get(trans.getLabel()) >= 0) {
						//							lp.setMat(0, idx, maxLength
						//									* costs[label2short.get(trans.getLabel())][(idx - 1) / transitions] - 1);
						//						}
						lp.setColName(idx, trans.getLabel() + "," + (idx - 1) / transitions);
					}
					lp.setRowName(i, "A " + (i - 1) / places + "," + ((i - 1) % places));
					if (i <= maxLength * places) {
						// at least -current marking
						lp.setConstrType(i, LpSolve.GE);
					} else {
						// equal to final marking - output of xi?
						lp.setConstrType(i, LpSolve.EQ);
					}
				}
			}
		}
		int row = places * (maxLength + 1) + 1;
		// set up the comparisons and sums
		for (int i = 0; i < maxLength - 1; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t + 1, trans2label[t % transitions] <= 0 ? 0 : 1);
			}
			for (int t = (i + 1) * transitions; t < (i + 2) * transitions; t++) {
				lp.setMat(row + i, t + 1, trans2label[t % transitions] <= 0 ? 0 : -1);
			}
			lp.setRowName(row + i, "X" + i + "-X" + (i + 1) + ">0");
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// sos constraints, sum max 1
		row = places * (maxLength + 1) + maxLength;
		for (int i = 0; i < maxLength; i++) {
			for (int t = i * transitions; t < (i + 1) * transitions; t++) {
				lp.setMat(row + i, t + 1, trans2label[t % transitions] <= 0 ? 0 : 1);
			}
			lp.setRowName(row + i, "X" + i + ".1 <=1");
			lp.setConstrType(row + i, LpSolve.LE);
		}

		row = places * (maxLength + 1) + 2 * maxLength;
		for (int t = 0; t < log.length; t++) {
			// what's the hamming distance of the x's to trace t?
			for (int e = 0; e < maxLength; e++) {
				for (int tr = 0; tr < trans2label.length; tr++) {
					if (e >= log[t].length || (trans2label[tr] != log[t][e] && trans2label[tr] > 0)) {
						lp.setMat(row + t, e * transitions + tr + 1, 1);
					} else {
						// based on the index and the fact that this transition matches that label,
						// set the value of this trace's row to 0;
						lp.setMat(row + t, e * transitions + tr + 1, 0);
					}
				}
			}
			lp.setRowName(row + t, "trace_" + t);
			lp.setMat(row + t, transitions * maxLength + 1, -1);
			lp.setConstrType(row + t, LpSolve.GE);
		}

		row = (maxLength + 1) * places + 2 * maxLength + log.length;
		for (int t = 1; t <= transitions * maxLength; t++) {
			lp.setMat(row, t, trans2label[t % transitions] <= 0 ? 0 : 1);
			lp.setMat(0, t, -1);
		}
		lp.setRowName(row, "SumX <= " + maxLength);
		lp.setConstrType(row, LpSolve.LE);

		// maximize h
		lp.setMat(0, transitions * maxLength + 1, maxLength * maxLength);

		lp.setMaxim();

		lp.setVerbose(1);

		lp.setScaling(LpSolve.SCALE_GEOMETRIC | LpSolve.SCALE_EQUILIBRATE | LpSolve.SCALE_INTEGERS);
		lp.setScalelimit(5);
		lp.setPivoting(LpSolve.PRICER_DEVEX | LpSolve.PRICE_ADAPTIVE);
		lp.setMaxpivot(250);
		lp.setBbFloorfirst(LpSolve.BRANCH_AUTOMATIC);
		lp.setBbRule(LpSolve.NODE_PSEUDONONINTSELECT | LpSolve.NODE_GREEDYMODE | LpSolve.NODE_DYNAMICMODE
				| LpSolve.NODE_RCOSTFIXING);
		lp.setBbDepthlimit(-50);
		lp.setAntiDegen(LpSolve.ANTIDEGEN_FIXEDVARS | LpSolve.ANTIDEGEN_STALLING);
		lp.setImprove(LpSolve.IMPROVE_DUALFEAS | LpSolve.IMPROVE_THETAGAP);
		lp.setBasiscrash(LpSolve.CRASH_NOTHING);
		lp.setSimplextype(LpSolve.SIMPLEX_DUAL_PRIMAL);

		double[] rhs = new double[(maxLength + 1) * places + 2 * maxLength + log.length + 1];

		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row <= places * (maxLength + 1); row += places) {
				rhs[row] = -initialMarking.occurrences(p);
			}
		}
		for (Place p : finalMarking.baseSet()) {
			row = maxLength * places + place2int.get(p);
			rhs[row] += finalMarking.occurrences(p);
		}
		for (row = (maxLength + 1) * places + maxLength; row < (maxLength + 1) * places + 2 * maxLength; row++) {
			rhs[row] = 1;
		}
		row = (maxLength + 1) * places + 2 * maxLength;
		for (int t = 0; t < log.length; t++) {
			if (log[t].length > maxLength) {
				// the remaining hamming distance is the part of the trace not covered by the 
				// maxLength
				rhs[row + t] = maxLength - log[t].length;
			}
		}
		row = (maxLength + 1) * places + 2 * maxLength + log.length;
		rhs[row] = maxLength;

		lp.setRhVec(rhs);

		lp.printLp();

		return lp;
	}

	private double getCost(Transition trans, int i) {
		return 1;
	}

}
