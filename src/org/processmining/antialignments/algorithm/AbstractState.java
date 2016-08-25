package org.processmining.antialignments.algorithm;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Vector;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

public abstract class AbstractState implements State {

	protected final Marking marking;
	protected final int[] vector;
	protected final AbstractState predecessor;
	protected final Transition executedTransition;
	protected final short executedLabel;

	/**
	 * Constructs the initial state
	 */
	public AbstractState(Marking marking, int logSize) {
		this.marking = marking;
		this.vector = new int[logSize + 2];
		this.predecessor = null;
		this.executedTransition = null;
		this.executedLabel = NOLABEL;
	}

	protected AbstractState(Marking marking, AbstractState predecessor, int[] vector, Transition executedTransition,
			short executedLabel) {
		this.marking = marking;
		this.predecessor = predecessor;
		this.vector = vector;
		this.executedTransition = executedTransition;
		this.executedLabel = executedLabel;
	}

	@Override
	public Marking getMarking() {
		return marking;
	}

	@Override
	public int getLength() {
		return vector[0];
	}

	public String toString() {
		return (executedTransition == null || executedTransition.isInvisible() ? "." : executedTransition.getLabel())
				+ ": " + marking.toString() + " " + Arrays.toString(vector);
	}

	@Override
	public String getAntiAlignmentString() {
		if (predecessor == null) {
			return "";
		}
		return predecessor.getAntiAlignmentString()
				+ (executedTransition.isInvisible() ? "" : ";" + executedTransition.getLabel());
	}

	@Override
	public short[] getAntiAlignment() {
		return buildAntiAlignment(0);
	}

	private short[] buildAntiAlignment(int length) {
		short[] aa;
		if (predecessor == null) {
			aa = new short[length];
		} else {
			if (executedTransition.isInvisible()) {
				aa = predecessor.buildAntiAlignment(length);
			} else {
				length++;
				aa = predecessor.buildAntiAlignment(length);
				aa[aa.length - length] = executedLabel;
			}
		}
		return aa;
	}

	@Override
	public Vector<Transition> getFiringSequence() {
		return buildFiringSequence();
	}

	private Vector<Transition> buildFiringSequence() {
		Vector<Transition> fs;
		if (predecessor == null) {
			fs = new Vector<Transition>();
		} else {
			fs = predecessor.buildFiringSequence();
			fs.add(executedTransition);
		}
		return fs;
	}

	@Override
	public void printMatrix(PrintStream out) {
		if (predecessor != null) {
			predecessor.printMatrix(out);
			if (executedTransition.isInvisible()) {
				out.print(".");
			} else {
				out.print(executedTransition.getLabel());
			}
		} else {
			out.print(" ");
		}
		out.print(": ");
		out.println(Arrays.toString(vector));

	}

	public int getPathDistance() {
		return vector[0] - vector[vector.length - 1];
	}

	public void setFinalMarkingReached(short[][] log, int traceToIgnore) {
		int min = Integer.MAX_VALUE;
		int maxDist = 0;
		for (int t = 0; t < log.length; t++) {
			int dist = log[t].length - vector[t + 1];
			if (dist > 0) {
				// some trailing tau steps are needed
				vector[t + 1] += dist;
			}
			if (dist > maxDist) {
				maxDist = dist;
			}
			if (t != traceToIgnore && vector[t + 1] < min) {
				min = vector[t + 1];
			}
		}
		// update the path length
		vector[0] += maxDist;
		// update the minimal value
		vector[vector.length - 1] = min;
	}

}