package org.processmining.antialignments.ilp.util;

import gnu.trove.list.TIntList;
import gnu.trove.list.TShortList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;

import java.util.Collection;

/**
 * The class represents a trace from a TShortList
 * 
 * @author bfvdonge
 * 
 */
public class Representative {
	private final int number;
	private final TIntList represented = new TIntArrayList(10);
	private final TShortList trace;

	public Representative(TShortList trace, int number) {
		this.trace = trace;
		this.number = number;
	}

	public void addRepresentedTrace(int trace) {
		represented.add(trace);
	}

	public TIntList getRepresented() {
		return represented;
	}

	public TShortList getTrace() {
		return trace;
	}

	public boolean equals(Object o) {
		return o != null && (o instanceof Representative ? ((Representative) o).trace.equals(trace) : false);
	}

	public int hashCode() {
		return trace.hashCode();
	}

	public String toString() {
		return trace.toString() + " representing " + represented.toString();
	}

	public int getNumber() {
		return number;
	}

	public void addRepresentedTrace(Collection<? extends Integer> traceIndices) {
		represented.addAll(traceIndices);

	}

	public void addRepresentedTrace(TIntSet traceIndices) {
		represented.addAll(traceIndices);

	}

}