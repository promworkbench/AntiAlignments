package org.processmining.antialignments.ilp.util;

import org.processmining.models.semantics.petrinet.Marking;

public class HybridEquationResult {
	private Marking marking;
	private final int lengthX;
	private final int lengthY;

	public HybridEquationResult(Marking marking, int lengthX, int lengthY) {
		this.setMarking(marking);
		this.lengthX = lengthX;
		this.lengthY = lengthY;
	}

	public HybridEquationResult(int lengthX, int lengthY) {

		this.lengthX = lengthX;
		this.lengthY = lengthY;
	}

	public Marking getMarking() {
		return marking;
	}

	public int getLengthX() {
		return lengthX;
	}

	public int getLengthY() {
		return lengthY;
	}

	public boolean equals(Object o) {
		if (o != null && o instanceof HybridEquationResult) {
			HybridEquationResult h = (HybridEquationResult) o;
			return h.lengthX == lengthX && h.lengthY == lengthY && h.getMarking().equals(getMarking());
		}
		return false;
	}

	public int hashCode() {
		return getMarking().hashCode() + 31 * lengthX + 31 * 31 * lengthY;
	}

	public String toString() {
		return getMarking().toString() + " Lx: " + lengthX + " Ly: " + lengthY;
	}

	public void setMarking(Marking marking) {
		this.marking = marking;
	}

}