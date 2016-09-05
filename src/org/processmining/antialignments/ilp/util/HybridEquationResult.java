package org.processmining.antialignments.ilp.util;

import gnu.trove.list.TShortList;

import java.util.Stack;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

public class HybridEquationResult {
	private final Marking marking;
	final int lengthX;
	final int lengthY;

	public HybridEquationResult(Marking marking, int lengthX, int lengthY) {
		this.marking = marking;
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
			return h.lengthX == lengthX && h.lengthY == lengthY && h.getMarking().equals(marking);
		}
		return false;
	}

	public int hashCode() {
		return marking.hashCode() + 31 * lengthX + 31 * 31 * lengthY;
	}

	public String toString() {
		return marking.toString() + " Lx: " + lengthX + " Ly: " + lengthY;
	}

	public void undo(TShortList antiAlignment, Stack<Transition> firingSequence) {
		// remove lengthX elements from anitAlignment
		antiAlignment.remove(antiAlignment.size() - lengthX, lengthX);
		int popped = 0;
		while (firingSequence.peek().isInvisible() || popped < lengthX) {
			popped += firingSequence.pop().isInvisible() ? 0 : 1;
		}
	}
}