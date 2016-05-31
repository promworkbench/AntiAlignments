package org.processmining.antialignments.pathfinder;

import java.util.Arrays;
import java.util.Vector;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;

public class AntiAlignments {

	private final double[] maxMinDistances;
	private final int[] traceDistances;
	private final short[][] antiAlignments;
	private final Vector<?>[] traces;
	private final double[] maxDistances;

	public AntiAlignments(int logLength) {
		maxMinDistances = new double[logLength + 1];
		Arrays.fill(getMaxMinDistances(), -1);
		traceDistances = new int[logLength];
		antiAlignments = new short[logLength + 1][];
		traces = new Vector<?>[logLength + 1];
		maxDistances = new double[logLength + 1];
	}

	public int getLogLength() {
		return getTraceDistances().length;
	}

	public short[] getAAForLog() {
		return getAntiAlignments()[getAntiAlignments().length - 1];
	}

	public short[] getAAForLogWithoutTrace(int trace) {
		return getAntiAlignments()[trace];
	}

	public Vector<Transition> getAAFiringSequenceForLog() {
		return (Vector<Transition>) getTraces()[getAntiAlignments().length - 1];
	}

	public Vector<Transition> getAAFiringSequenceForLogWithoutTrace(int trace) {
		return (Vector<Transition>) getTraces()[trace];
	}

	public double getAADistanceForLog() {
		return getMaxMinDistances()[getMaxMinDistances().length - 1];
	}

	public double getAADistanceForLogWithoutTrace(int trace) {
		return getMaxMinDistances()[trace];
	}

	public int getAADistanceToTrace(int trace) {
		return getTraceDistances()[trace];
	}

	public double getMaxDistanceForTrace(int t) {
		return getMaxDistances()[t];
	}

	public double getMaxDistanceForLog() {
		return getMaxDistances()[getMaxDistances().length - 1];
	}

	public double[] getMaxDistances() {
		return maxDistances;
	}

	public double[] getMaxMinDistances() {
		return maxMinDistances;
	}

	public short[][] getAntiAlignments() {
		return antiAlignments;
	}

	public Vector<?>[] getTraces() {
		return traces;
	}

	public int[] getTraceDistances() {
		return traceDistances;
	}

}