package org.processmining.antialignments.ilp.antialignment;

public class AntiAlignmentParameters {

	private final int cutOffLength;

	private final double maxFactor;

	public AntiAlignmentParameters(int cutOffLength, double maxFactor) {
		this.cutOffLength = cutOffLength;
		this.maxFactor = maxFactor;

	}

	public int getCutOffLength() {
		return cutOffLength;
	}

	public double getMaxFactor() {
		return maxFactor;
	}

}
