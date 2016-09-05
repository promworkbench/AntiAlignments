package org.processmining.antialignments.bruteforce;

import java.io.PrintStream;
import java.util.Arrays;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

/**
 * A state is a pair of a marking and a row in the edit-distance matrix.
 * 
 * A state has a predecessor state (or null if it is the intial state)
 * 
 * @author bfvdonge
 * 
 */
public class HammingState extends AbstractState {

	protected final short[] vector;
	protected short minimumDistance;
	protected short length;

	/**
	 * Constructs the initial state
	 */
	public HammingState(Marking marking, int logSize) {
		super(marking);
		this.vector = new short[logSize];
		this.minimumDistance = 0;
		this.length = 0;
	}

	private HammingState(Marking marking, HammingState predecessor, short[] vector, short minimumDistance,
			short length, Transition executedTransition, short executedLabel) {
		super(marking, predecessor, executedTransition, executedLabel);
		this.vector = vector;
		this.minimumDistance = minimumDistance;
		this.length = length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.processmining.antialignments.algorithm.State#getNextState(org.
	 * processmining.models.semantics.petrinet.Marking, short[][], int, short,
	 * org
	 * .processmining.models.graphbased.directed.petrinet.elements.Transition)
	 */
	@Override
	public HammingState getNextState(Marking newMarking, final short[][] log, final int traceToIgnore,
			short executedLabel, Transition executedTransition) {
		short[] newVector = new short[vector.length];
		short newLength = length;
		short newMinimumDistance = minimumDistance;

		if (executedLabel == NOLABEL) {
			// invisible transition, 
			for (int i = newVector.length; i-- > 0;) {
				newVector[i] = vector[i];
			}
		} else {
			// store the length of the backtrace
			newLength = (short) (length + 1);
			// store the minimum
			newMinimumDistance = Short.MAX_VALUE;
			for (int t = 0; t < log.length; t++) {
				if (length < log[t].length && log[t][length] == executedLabel) {
					// Match
					newVector[t] = vector[t];
				} else {
					// no match or passed trace t's length				
					newVector[t] = (short) (vector[t] + 1);
				}
				// Keep track of the minimum value in the last element of the vector.
				// over all traces nog equal to traceToIgnore.
				if (t != traceToIgnore && newVector[t] < newMinimumDistance) {
					newMinimumDistance = newVector[t];
				}
			}
		}

		return new HammingState(newMarking, this, newVector, newMinimumDistance, newLength, executedTransition,
				executedLabel);

	}

	public boolean equals(Object o) {
		if (o instanceof HammingState) {
			HammingState s = (HammingState) o;
			return s.marking.equals(marking) && length == s.length && Arrays.equals(s.vector, vector);
		} else {
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.processmining.antialignments.algorithm.State#hasSameMarkingAs(org
	 * .processmining.antialignments.algorithm.HammingState)
	 */
	@Override
	public boolean hasSameMarkingAs(State s) {
		return s instanceof HammingState && marking.equals(((HammingState) s).marking);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.processmining.antialignments.algorithm.State#isGreaterOrEqual(org
	 * .processmining.antialignments.algorithm.HammingState)
	 */
	@Override
	public boolean isGreaterOrEqual(State s) {
		if (!(s instanceof HammingState) || length != ((HammingState) s).length) {
			return false;
		}
		for (int i = 1; i < vector.length - 1; i++) {
			if (vector[i] < ((HammingState) s).vector[i]) {
				return false;
			}
		}
		return true;
	}

	public int hashCode() {
		return marking.hashCode() + 31 * Arrays.hashCode(vector) + 31 * 31 * length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.processmining.antialignments.algorithm.State#getMinimalDistance(short
	 * [][], int)
	 */
	@Override
	public int getMinimalDistance(short[][] log, int traceToIgnore) {
		return minimumDistance;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.processmining.antialignments.algorithm.State#getDistance(int)
	 */
	@Override
	public int getDistance(short[][] log, int trace) {
		return vector[trace] + (log[trace].length > length ? log[trace].length - length : 0);
	}

	public int getPathDistance() {
		return length - minimumDistance;
	}

	public void setFinalMarkingReached(short[][] log, int traceToIgnore) {
		short min = Short.MAX_VALUE;
		short maxDist = 0;
		for (int t = 0; t < log.length; t++) {
			short dist = (short) (log[t].length - vector[t]);
			if (dist > 0) {
				// some trailing tau steps are needed
				vector[t] += dist;
			}
			if (dist > maxDist) {
				maxDist = dist;
			}
			if (t != traceToIgnore && vector[t] < min) {
				min = vector[t];
			}
		}
		// update the path length
		length += maxDist;
		// update the minimal value
		minimumDistance = min;
	}

	@Override
	public int getLength() {
		return length;
	}

	public String toString() {
		return (executedTransition == null || executedTransition.isInvisible() ? "." : executedTransition.getLabel())
				+ ": " + marking.toString() + " " + Arrays.toString(vector) + " l:" + length + " d:"
				+ getPathDistance();
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

}
