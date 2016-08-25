package org.processmining.antialignments.algorithm;

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

	/**
	 * Constructs the initial state
	 */
	public HammingState(Marking marking, int logSize) {
		super(marking, logSize);
	}

	private HammingState(Marking marking, HammingState predecessor, int[] vector, Transition executedTransition,
			short executedLabel) {
		super(marking, predecessor, vector, executedTransition, executedLabel);
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
		int[] newVector = new int[vector.length];

		if (executedLabel == NOLABEL) {
			// invisible transition, 
			for (int i = newVector.length; i-- > 0;) {
				newVector[i] = vector[i];
			}
		} else {
			// store the length of the backtrace
			newVector[0] = vector[0] + 1;
			// store the minimum
			newVector[newVector.length - 1] = Integer.MAX_VALUE;
			for (int t = 0; t < log.length; t++) {
				if (vector[0] < log[t].length && log[t][vector[0]] == executedLabel) {
					// Match
					newVector[t + 1] = vector[t + 1];
				} else {
					// no match or passed trace t's length				
					newVector[t + 1] = vector[t + 1] + 1;
				}
				// Keep track of the minimum value in the last element of the vector.
				// over all traces nog equal to traceToIgnore.
				if (t != traceToIgnore && newVector[t + 1] < newVector[newVector.length - 1]) {
					newVector[newVector.length - 1] = newVector[t + 1];
				}
			}
		}

		return new HammingState(newMarking, this, newVector, executedTransition, executedLabel);

	}

	public boolean equals(Object o) {
		if (o instanceof HammingState) {
			HammingState s = (HammingState) o;
			return s.marking.equals(marking) && Arrays.equals(s.vector, vector);
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
		if (!(s instanceof HammingState) || vector[0] != ((HammingState) s).vector[0]) {
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
		return marking.hashCode() + 31 * Arrays.hashCode(vector);
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
		int min = Integer.MAX_VALUE;
		for (int t = 0; t < log.length; t++) {
			int dist = getDistance(log, t);
			if (t != traceToIgnore && dist < min) {
				min = dist;
			}
		}
		assert min == vector[vector.length - 1];
		return min;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.processmining.antialignments.algorithm.State#getDistance(int)
	 */
	@Override
	public int getDistance(short[][] log, int trace) {
		return vector[trace + 1] + (log[trace].length > vector[0] ? log[trace].length - vector[0] : 0);
	}

}
