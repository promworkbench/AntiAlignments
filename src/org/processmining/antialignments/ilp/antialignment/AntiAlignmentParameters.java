package org.processmining.antialignments.ilp.antialignment;

public class AntiAlignmentParameters {

	private final int cutOffLength;

	private final double maxFactor;

	private final int backtrackLimit;

	private final double backtrackThreshold;

	public AntiAlignmentParameters(int cutOffLength, double maxFactor, int backtrackLimit, double backtrackThreshold) {
		this.cutOffLength = cutOffLength;
		this.maxFactor = maxFactor;
		this.backtrackLimit = backtrackLimit;
		this.backtrackThreshold = backtrackThreshold;

	}

	public int getCutOffLength() {
		return cutOffLength;
	}

	public double getMaxFactor() {
		return maxFactor;
	}

	public int getBacktrackLimit() {
		return backtrackLimit;
	}

	public double getBacktrackThreshold() {
		return backtrackThreshold;
	}

}
