package org.processmining.antialignments.algorithm;

import java.io.PrintStream;
import java.util.Arrays;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

public class EditDistanceState extends AbstractState {

	protected final short[][] vector;
	protected short minimumDistance;
	protected short length;

	/**
	 * Constructs the initial state
	 */
	public EditDistanceState(Marking marking, short[][] log) {
		super(marking);
		this.vector = new short[log.length][];
		for (int t = 0; t < log.length; t++) {
			this.vector[t] = new short[log[t].length];
			for (int e = 0; e < log[t].length; e++) {
				this.vector[t][e] = (short) (e + 1);
			}
		}
		this.length = 0;
		this.minimumDistance = 0;
	}

	private EditDistanceState(Marking marking, EditDistanceState predecessor, short[][] vector, short minimumDistance,
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
	public EditDistanceState getNextState(Marking newMarking, final short[][] log, final int traceToIgnore,
			short executedLabel, Transition executedTransition) {
		short[][] newVector;

		short newLength = length;
		short newMinDistance = minimumDistance;
		if (executedLabel == NOLABEL) {
			// invisible transition,
			newVector = vector; // POINTERS are ok here as the previous state won't change.
		} else {
			newVector = new short[vector.length][];
			// store the length of the backtrace
			newLength = (short) (length + 1);
			// store the minimum
			newMinDistance = Short.MAX_VALUE;

			for (int t = 0; t < log.length; t++) {
				newVector[t] = new short[log[t].length];
				if (log[t][0] == executedLabel) {
					newVector[t][0] = length;
				} else {
					newVector[t][0] = newLength;
					if (vector[t][0] < length) {
						newVector[t][0] = (short) (vector[t][0] + 1);
					}
				}
				for (int e = 1; e < log[t].length; e++) {
					if (log[t][e] == executedLabel) {
						newVector[t][e] = vector[t][e - 1];
					} else {
						newVector[t][e] = (short) (newVector[t][e - 1] + 1);
						if (vector[t][e - 1] < newVector[t][e - 1]) {
							newVector[t][e] = (short) (vector[t][e - 1] + 1);
						}
						if (vector[t][e] < vector[t][e - 1] && vector[t][e] < newVector[t][e - 1]) {
							newVector[t][e] = (short) (vector[t][e - 1] + 1);

						}
					}
				}
				// Keep track of the minimum value
				// over all traces nog equal to traceToIgnore.
				short distTrace = newVector[t][Math.min(log[t].length - 1, length)];
				if (t != traceToIgnore && distTrace < newMinDistance) {
					newMinDistance = distTrace;
				}
			}

		}

		return new EditDistanceState(newMarking, this, newVector, newMinDistance, newLength, executedTransition,
				executedLabel);

	}

	public boolean equals(Object o) {
		if (o instanceof EditDistanceState) {
			EditDistanceState s = (EditDistanceState) o;
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
		return s instanceof EditDistanceState && marking.equals(((EditDistanceState) s).marking);
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
		if (!(s instanceof EditDistanceState) || length != ((EditDistanceState) s).length) {
			return false;
		}
		for (int t = 0; t < vector.length - 1; t++) {
			for (int e = 0; e < vector[t].length && e < length; e++) {
				if (vector[t][e] < ((EditDistanceState) s).vector[t][e]) {
					return false;
				}
			}
		}
		return true;
	}

	public int hashCode() {
		return marking.hashCode() + 31 * Arrays.hashCode(vector) + 31 * 31 * length + 31 * 31 * 31 * minimumDistance;
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
		return minimumDistance + (log[trace].length > length ? log[trace].length - length : 0);
	}

	public int getPathDistance() {
		return length - minimumDistance;
	}

	public void setFinalMarkingReached(short[][] log, int traceToIgnore) {
		short min = Short.MAX_VALUE;
		short maxLength = 0;
		for (int t = 0; t < log.length; t++) {

			if (t != traceToIgnore && log[t].length > maxLength) {
				maxLength = (short) log[t].length;
			}

			if (t != traceToIgnore && vector[t][log[t].length - 1] < min) {
				min = vector[t][log[t].length - 1];
			}
		}
		// update the path length
		length = maxLength;
		// update the minimal value
		minimumDistance = min;
	}

	@Override
	public int getLength() {
		return length;
	}

	public String toString() {
		return (executedTransition == null || executedTransition.isInvisible() ? "." : executedTransition.getLabel())
				+ ": " + marking.toString() + " " + Arrays.deepToString(vector) + " l:" + length + " d:"
				+ minimumDistance;
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
		out.println(Arrays.deepToString(vector));

	}
}