package org.processmining.antialignments.pathfinder;

import java.util.Arrays;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

/**
 * A state is a combination of a marking and the arrays for the edit distance
 * 
 * Since the minimum edit distance is non-decreasing, should use this for the
 * sort order in the priority queue. This is not implemented as the natural sort
 * order as it goes against the hascode and equals which are based on the
 * marking and the length of the partial anti-alignment.
 * 
 * delta stores the length of the search path so far, i.e. the number of
 * non-edits introduced by the path to the origin. This is what we aim to
 * MINIMIZE
 * 
 * @author bfvdonge
 * 
 */
public class State {

	private final Marking marking;
	private final int[][] current;
	private State predecessor;
	private final Transition fired;
	private final int firingSequenceLength;
	private double pathLength;

	public State(Marking marking, State predecessor, Transition fired,
			int[][] current, int firingSequenceLength, double pathLength) {
		this.marking = marking;
		this.predecessor = predecessor;
		this.fired = fired;
		this.current = current;

		this.firingSequenceLength = firingSequenceLength;
		this.pathLength = pathLength;
	}

	public int hashCode() {
		return marking.hashCode() + 37 * getAntiAlignmentLength();
	}

	public boolean equals(Object o) {
		if (!(o instanceof State)) {
			return false;
		}
		State s = (State) o;
		return getFiringSequenceLength() == s.getAntiAlignmentLength()
				&& s.marking.equals(marking);
		// return s.marking.equals(marking);
	}

	public String toString() {
		return marking.toString() + "l: " + pathLength;
	}

	public Marking getMarking() {
		return marking;
	}

	public int[][] getCurrent() {
		return current;
	}

	public int getAntiAlignmentLength() {
		return current[0][0];
	}

	public int getFiringSequenceLength() {
		return firingSequenceLength;
	}

	public State getPredecessor() {
		return predecessor;
	}

	public Transition getFiredTransition() {
		return fired;
	}

	public void setPredecessor(State predecessor) {
		this.predecessor = predecessor;
	}

	// public boolean isFinal() {
	// return
	// }

	public int getDistance(int t) {
		return current[t][current[t].length - 1];
	}

	public double getPathLength() {
		return pathLength;
	}

	public int getMinMinDistance(int traceToIgnore) {
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < current.length; i++) {
			if (i != traceToIgnore && getDistance(i) < min) {
				min = getDistance(i);
			}
		}
		return min;
	}

}
