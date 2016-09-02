package org.processmining.antialignments.alignments;

import gnu.trove.iterator.TShortIterator;
import gnu.trove.map.TObjectShortMap;
import gnu.trove.map.TShortObjectMap;
import gnu.trove.set.TShortSet;
import gnu.trove.set.hash.TShortHashSet;

import java.util.Arrays;

import lpsolve.LpSolve;

import org.processmining.antialignments.algorithm.AbstractILPCalculator;
import org.processmining.antialignments.algorithm.ilp.LPMatrix;
import org.processmining.antialignments.algorithm.ilp.LPMatrix.LPMatrixException;
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

public class AlignmentILPCalculator extends AbstractILPCalculator {

	public AlignmentILPCalculator(Petrinet net, Marking initialMarking, Marking finalMarking,
			TObjectShortMap<String> label2short, TShortObjectMap<String> short2label, short[][] log) {
		super(net, initialMarking, finalMarking, label2short, short2label, log);
	}

	public void getAlignments(Marking initialMarking, Marking finalMarking) throws LPMatrixException {

		mode = MODE_LPSOLVE;

		LPMatrix<?> matrix = setupLpForHybrid(3, true, initialMarking, finalMarking, 0, 0);
		((LpSolve) matrix.toSolver()).printLp();

	}

	protected LPMatrix<?> setupLpForHybrid(int maxLengthX, boolean integerVariables, Marking initialMarking,
			Marking finalMarking, int traceToConsider, int startTraceAt) {

		// Find the transitions in the net labeled with an event in 
		// log[traceToConsider][startTraceAt] .. log[traceToConsider][startTraceAt+lengthX-1]
		short[] mappedTransitions;

		short[] mappedLabels = new short[Math.min(log[traceToConsider].length - startTraceAt, maxLengthX)];

		// map2post stores the relative position of the synchronous move version of a transition
		short[] mapped2pos = new short[transitions];
		short[] label2pos = new short[labels];
		{
			Arrays.fill(label2pos, (short) -1);
			TShortSet trs = new TShortHashSet();
			short e;
			for (e = (short) startTraceAt; e < log[traceToConsider].length && e < startTraceAt + maxLengthX; e++) {
				for (short t = invisibleTransitions; t < transitions; t++) {
					if (equalLabel(t, log[traceToConsider][e])) {
						trs.add(t);
					}
				}
				mappedLabels[e - startTraceAt] = log[traceToConsider][e];
				label2pos[log[traceToConsider][e]] = (short) (e - startTraceAt);
			}
			mappedTransitions = new short[trs.size()];

			// position the other labels
			for (int i = 0; i < label2pos.length; i++) {
				if (label2pos[i] == -1) {
					label2pos[i] = e;
					e++;
				}
			}

			short p = 0;
			TShortIterator it = trs.iterator();
			while (it.hasNext()) {
				short t = it.next();
				mappedTransitions[p] = t;

				mapped2pos[t] = (short) (transitions - t + p);
				p++;
			}
		}

		int visibleTransitions = mappedTransitions.length;
		int spCols = transitions + visibleTransitions + mappedLabels.length;
		int spRows = places + mappedLabels.length + 1;

		//-  
		//-                            AB   >= m0
		//-                            AAB  >= m0
		//-                            AAAA == mf        
		LPMatrix<?> lp = setupMatrix(spRows * maxLengthX + (places + labels) + //
				// Xi - Xi+1<=0        Xi <= 1  Y <= ly*X_lx
				(maxLengthX - 1) + maxLengthX + 1,//
				//  X part is small           Y part is full
				spCols * maxLengthX + 2 * transitions - invisibleTransitions + labels);

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
			// pos is the relative position of the column index of t if it is a synchronous move
			int pos = mapped2pos[t];
			// label is the relative position of the corresponding label row.
			int labelRow = places - p + label2pos[label2short.get(trans.getLabel())];

			for (int block = 0; block <= maxLengthX; block++) {
				// First the whole A matrices.
				int r = block * spRows + p;

				lp.setConstrType(r, type);
				lp.setRowName(r, "A" + block + "_" + p);
				for (int c = t; c < block * spCols; c += spCols) {
					// update the SP A matrix
					lp.setMat(r, c, lp.getMat(r, c) + dir);
					// update the SP A' matrix for mapped visible transitions
					if (pos > 0) {
						lp.setMat(r, c + pos, lp.getMat(r, c + pos) + dir);
						if (block < maxLengthX) {
							// token flow
							lp.setMat(r + labelRow, c + pos, -1);
							lp.setMat(r + labelRow + 1, c + pos, 1);
						} else {
							// count
							lp.setMat(r + labelRow, c + pos, 1);
						}
					}
				}

				// Then, the  A- matrix.
				int c = block * spCols + t;
				if ((dir < 0 || trans.isInvisible()) && block < maxLengthX) {
					lp.setColName(c, "M_" + trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
					lp.setInt(c, integerVariables);
					if (!trans.isInvisible()) {
						// set upper bound of 1 for visible transitions
						lp.setUpbo(c, 1.0);
					}
					if (pos > 0) {
						lp.setColName(c + pos, "S_" + trans.getLabel() + "-" + block);
						lp.setObjective(c + pos, 1);

						lp.setInt(c + pos, integerVariables);
						lp.setUpbo(c + pos, 1.0);
						lp.setMat(r, c + pos, lp.getMat(r, c + pos) + dir);
						lp.setMat(r + labelRow, c + pos, -1);
						// lp.setMat(r + label + 1, c + pos, 1);
					}
					// update the A matrix only for consumption and for invisible transitions
					// or in the last block

					// update the SP A matrix for consumption in case of visible,
					// or always in case of invisible transition
					lp.setMat(r, c, lp.getMat(r, c) + dir);
				}

				// Then, in the last block
				if (block == maxLengthX) {

					// In this Block, the positions of the labels are now different,
					// as are the positions of the synchronous transitions.

					// The relative position of the synchronous counterpart is at:
					pos = transitions - invisibleTransitions;

					if (finalMarking != null) {
						lp.setConstrType(r, LpSolve.EQ);
					}
					// update the SP A matrix
					lp.setMat(r, c, lp.getMat(r, c) + dir);
					lp.setColName(c, "M_" + trans.getLabel().replace("\\n$invisible$", "") + "-" + block);
					lp.setInt(c, integerVariables);

					lp.setColName(c + pos, "S_" + trans.getLabel() + "-" + block);
					lp.setObjective(c + pos, 1);

					lp.setInt(c + pos, integerVariables);
					lp.setMat(r, c + pos, lp.getMat(r, c + pos) + dir);
					lp.setMat(r + labelRow, c + pos, 1);

				}

			}
		}

		for (short l = 0; l < mappedLabels.length; l++) {
			for (int block = 0; block <= maxLengthX; block++) {
				// First the whole A matrices.
				int r = block * spRows + places + l;
				lp.setConstrType(r, LPMatrix.GE);
				lp.setRowName(r, "B" + block + "_" + l);
				if (l == mappedLabels.length - 1) {
					lp.setConstrType(r + 1, LPMatrix.GE);
					lp.setRowName(r + 1, "B" + block + "_" + (l + 1));
				}
				for (int c = transitions + visibleTransitions + l; c < block * spCols; c += spCols) {
					// update the SP B matrix
					if (block < maxLengthX) {
						// token flow
						lp.setMat(r, c, -1);
						lp.setMat(r + 1, c, 1);
					} else {
						// count
						lp.setMat(r, c, 1);
					}
					lp.setUpbo(c, 1);
					lp.setInt(c, integerVariables);
				}

				// Then, the  A- matrix.
				int c = block * spCols + transitions + visibleTransitions + l;
				if (block < maxLengthX) {
					lp.setColName(c, "L_" + short2label.get(mappedLabels[l]) + "-" + block);
					lp.setInt(c, integerVariables);
					lp.setUpbo(c, 1.0);
					lp.setMat(r, c, -1);
				}

			}
		}
		for (short l = 0; l < label2pos.length; l++) {
			int r = maxLengthX * spRows + places + label2pos[l];
			lp.setRowName(r, "C_" + short2label.get(l));

			int c = maxLengthX * spCols + 2 * transitions - invisibleTransitions + label2pos[l];

			// Then, in the last block
			lp.setConstrType(r, LpSolve.EQ);
			// count the number of B's
			lp.setMat(r, c, 1);
			lp.setColName(c, "L_" + short2label.get(l));
			lp.setInt(c, integerVariables);

		}

		int row = spRows * maxLengthX + places + labels;
		// set up the comparisons and sums
		for (int i = 0; i < maxLengthX - 1; i++) {
			for (int t = i * spCols + invisibleTransitions; t < (i + 1) * spCols; t++) {
				lp.setMat(row + i, t, 1);
			}
			for (int t = (i + 1) * spCols + invisibleTransitions; t < (i + 2) * spCols; t++) {
				lp.setMat(row + i, t, -1);
			}
			lp.setRowName(row + i, "X" + i + "~X" + (i + 1));
			lp.setConstrType(row + i, LpSolve.GE);
		}

		// Y < lm*l_lx
		row += maxLengthX - 1;
		for (int t = (maxLengthX - 1) * spCols + invisibleTransitions; t < maxLengthX * spCols; t++) {
			//TODO: Find some constant!
			lp.setMat(row, t, transitions * transitions);
		}
		for (int t = maxLengthX * spCols + invisibleTransitions; t < lp.getNcolumns(); t++) {
			lp.setMat(row, t, -1);
		}
		lp.setRowName(row, "Y~X" + maxLengthX);
		lp.setConstrType(row, LpSolve.GE);

		// sos constraints, sum max 1
		row++;
		for (int i = 0; i < maxLengthX; i++) {
			for (int t = i * spCols + invisibleTransitions; t < (i + 1) * spCols; t++) {
				lp.setMat(row + i, t, 1);
			}
			lp.setRowName(row + i, "X" + i + ".1");
			lp.setConstrType(row + i, LpSolve.LE);
		}

		//TODO: Cost Function!
		lp.setMaxim();

		double[] rhs = new double[lp.getNrows()];

		// first -initial Marking
		for (Place p : initialMarking.baseSet()) {
			for (row = place2int.get(p); row < spRows * maxLengthX + places; row += spRows) {
				rhs[row] = -initialMarking.occurrences(p);
			}
			if (log[traceToConsider].length > startTraceAt) {
				for (row = places + label2pos[log[traceToConsider][0]]; row < spRows * maxLengthX; row += spRows) {
					rhs[row] = -1;
				}
			}
		}
		// then add final marking
		if (finalMarking != null) {
			for (Place p : finalMarking.baseSet()) {
				row = maxLengthX * spRows + place2int.get(p);
				rhs[row] += finalMarking.occurrences(p);
			}
		}
		// then, the parikh vector of the log
		row = maxLengthX * spRows + places;
		// Count the number of times each label is in the trace from position startTraceAt
		for (int e = startTraceAt; e < log[traceToConsider].length; e++) {
			rhs[row + label2pos[log[traceToConsider][e]]]++;
		}

		// then X.i - X.i+1 >= 0
		// then Y   - ly * X.lx >= 0

		// then X.1 <= 1
		for (row = maxLengthX * spRows + places + labels + maxLengthX; row < rhs.length; row++) {
			rhs[row] = 1;
		}

		lp.setRhVec(rhs);

		return lp;
	}

	protected void solveSequential(Marking initialMarking, Marking finalMarking, final int traceToConsider)
			throws LPMatrixException {

		PetrinetSemantics semantics = PetrinetSemanticsFactory.regularPetrinetSemantics(Petrinet.class);
		semantics.initialize(net.getTransitions(), initialMarking);

		Marking marking = initialMarking;
		int startTracesAt = 0;

		LPMatrix<?> matrix;

		if (VERBOSE) {
			System.out.println("Trying to get from " + initialMarking + " to " + finalMarking + //
					" starting with  " + cutOffLength + " exact steps.");
		}

		matrix = setupLpForHybrid(cutOffLength, true, marking, finalMarking, traceToConsider, startTracesAt);

		double[] vars = new double[matrix.getNcolumns()];
		matrix.solve(vars);

		// translate vars into log, sync and model moves
		System.out.println(Arrays.toString(matrix.getColNames()));
		System.out.println(Arrays.toString(vars));

	}

}
