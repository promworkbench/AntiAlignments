package org.processmining.antialignments;

import gnu.trove.list.TShortList;

import java.util.List;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;

public class AlignedTrace {

	private final TShortList modelSequence;
	private final List<Transition> firingSequence;

	public AlignedTrace(TShortList modelSequence, List<Transition> firingSequence) {
		this.modelSequence = modelSequence;
		this.firingSequence = firingSequence;

	}

	public List<Transition> getFiringSequence() {
		return firingSequence;
	}

	public TShortList getModelSequence() {
		return modelSequence;
	}

	public boolean equals(Object o) {
		return o instanceof AlignedTrace && ((AlignedTrace) o).firingSequence.equals(firingSequence);
	}

	public int hashCode() {
		return firingSequence.hashCode();
	}

}
