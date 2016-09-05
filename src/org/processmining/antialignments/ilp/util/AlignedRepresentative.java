package org.processmining.antialignments.ilp.util;

import gnu.trove.list.TShortList;

import java.util.List;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;

public class AlignedRepresentative extends Representative {

	private final List<Transition> firingSequence;

	public AlignedRepresentative(TShortList modelSequence, int number, List<Transition> firingSequence) {
		super(modelSequence, number);
		this.firingSequence = firingSequence;

	}

	public List<Transition> getFiringSequence() {
		return firingSequence;
	}

	public boolean equals(Object o) {
		return o instanceof AlignedRepresentative && ((AlignedRepresentative) o).firingSequence.equals(firingSequence);
	}

	public int hashCode() {
		return firingSequence.hashCode();
	}

}
