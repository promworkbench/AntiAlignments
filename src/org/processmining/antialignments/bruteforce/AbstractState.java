package org.processmining.antialignments.bruteforce;

import java.util.Vector;

import org.processmining.models.graphbased.directed.petrinet.elements.Transition;
import org.processmining.models.semantics.petrinet.Marking;

public abstract class AbstractState implements State {

	protected final Marking marking;
	protected final AbstractState predecessor;
	protected final Transition executedTransition;
	protected final short executedLabel;

	/**
	 * Constructs the initial state
	 */
	public AbstractState(Marking marking) {
		this.marking = marking;
		this.predecessor = null;
		this.executedTransition = null;
		this.executedLabel = NOLABEL;
	}

	protected AbstractState(Marking marking, AbstractState predecessor, Transition executedTransition,
			short executedLabel) {
		this.marking = marking;
		this.predecessor = predecessor;
		this.executedTransition = executedTransition;
		this.executedLabel = executedLabel;
	}

	@Override
	public Marking getMarking() {
		return marking;
	}

	@Override
	public String getAntiAlignmentString() {
		if (predecessor == null) {
			return "";
		}
		return predecessor.getAntiAlignmentString()
				+ (executedTransition.isInvisible() ? "" : ";" + executedTransition.getLabel());
	}

	@Override
	public short[] getAntiAlignment() {
		return buildAntiAlignment(0);
	}

	private short[] buildAntiAlignment(int length) {
		short[] aa;
		if (predecessor == null) {
			aa = new short[length];
		} else {
			if (executedTransition.isInvisible()) {
				aa = predecessor.buildAntiAlignment(length);
			} else {
				length++;
				aa = predecessor.buildAntiAlignment(length);
				aa[aa.length - length] = executedLabel;
			}
		}
		return aa;
	}

	@Override
	public Vector<Transition> getFiringSequence() {
		return buildFiringSequence();
	}

	private Vector<Transition> buildFiringSequence() {
		Vector<Transition> fs;
		if (predecessor == null) {
			fs = new Vector<Transition>();
		} else {
			fs = predecessor.buildFiringSequence();
			fs.add(executedTransition);
		}
		return fs;
	}

}