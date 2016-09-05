package org.processmining.antialignments.bruteforce;

import java.io.PrintStream;
import java.util.Vector;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

public interface State {

	public static final short NOLABEL = Short.MIN_VALUE;

	public abstract <S extends State> S getNextState(Marking newMarking, short[][] log, int traceToIgnore,
			short executedLabel, Transition executedTransition);

	/**
	 * Checks if the marking equals the marking of State s
	 * 
	 * @param s
	 * @return
	 */
	public abstract boolean hasSameMarkingAs(State s);

	/**
	 * returns true if and only if for all elements of the edit distance vector,
	 * the value of this is greater or equal to the value of s.
	 * 
	 * Note that s.isGreaterThan(s) returns true;
	 * 
	 * @param s
	 * @return
	 */
	public abstract boolean isGreaterOrEqual(State s);

	public abstract Marking getMarking();

	public abstract int getMinimalDistance(short[][] log, int traceToIgnore);

	public abstract int getDistance(short[][] log, int trace);

	public abstract int getLength();

	/**
	 * Returns the distance from the initial state to this state. This distance
	 * is equal to getLength() - the minimum hamming distance to any partial
	 * trace in the log upto length getLength()
	 * 
	 * @return
	 */
	public abstract int getPathDistance();

	public abstract String getAntiAlignmentString();

	public abstract short[] getAntiAlignment();

	public abstract Vector<Transition> getFiringSequence();

	public abstract void printMatrix(PrintStream out);

	public abstract void setFinalMarkingReached(short[][] log, int traceToIgnore);

}