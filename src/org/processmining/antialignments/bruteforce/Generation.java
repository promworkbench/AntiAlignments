package org.processmining.antialignments.bruteforce;

import java.util.Iterator;

import org.apache.commons.collections15.list.TreeList;

/**
 * A generation is a collection of states, such that: - For a state s, there is
 * no state s' in the collection with s.isGreaterOrEqual(s') true.
 * 
 * @author bfvdonge
 * 
 */

public class Generation<S extends State> implements Iterable<S> {

	//For now, just use a LinkedList of states, sorted by Marking.
	private TreeList<S> states = new TreeList<>();
	private S finalState;
	private int finalDistance;

	public Generation() {
		this(null, 0);
	}

	public Generation(S finalState, int minimalEditDistance) {
		this.finalState = finalState;
		this.finalDistance = minimalEditDistance;
	}

	public boolean setFinalState(S s, int minimalEditDistance) {
		//		System.out.println("Final state found: " + s.getAntiAlignmentString() + " distance: " + minimalEditDistance);
		if (finalState == null || minimalEditDistance > finalDistance) {
			finalState = s;
			finalDistance = minimalEditDistance;
			return true;
		}

		return false;
	}

	public boolean addNonFinalState(S s) {

		int key = s.getMarking().hashCode();

		int start = 0;
		int end = states.size() - 1;
		int mid = 0, midKey;
		while (start <= end) {
			mid = (start + end) / 2;
			midKey = states.get(mid).getMarking().hashCode();
			if (key == midKey) {
				// we are in a block with the same hashCode. Now look for the start of that block
				while (mid-- > 0 && (midKey = states.get(mid).getMarking().hashCode()) == key) {
				}

				if (midKey != key || mid < 0) {
					mid++;
				}
				assert (midKey = states.get(mid).getMarking().hashCode()) == key;

				// mid is the first element in states list of which the hashcode Matches the hashCode of the new state
				// now do a linear search.
				for (int m = mid; m < states.size();) {
					S stateAtM = states.get(m);
					if (stateAtM.getMarking().hashCode() != midKey) {
						// passed through the block. 
						break;
					}
					if (stateAtM.hasSameMarkingAs(s)) {
						if (s.equals(stateAtM)) {
							//							System.out.println("Not queueing " + s + " because it's already there.");
							// the new state is strictly smaller than an existing state, do not insert.
							return false;
						} else if (s.isGreaterOrEqual(stateAtM)) {
							// an existing state is strictly smaller than the new one. Remove!
							//							System.out.println("Removing " + stateAtM + " because of " + s);
							end--;
							states.remove(m);
						} else if (stateAtM.isGreaterOrEqual(s)) {
							//							System.out.println("Not queueing " + s + " because of " + stateAtM);
							// the new state is strictly smaller than an existing state, do not insert.
							return false;
						} else {
							m++;
						}

					} else {
						m++;
					}

				}
				// index mid is the insertion point as there is no reason to sort states within the
				// block in any order (there is no logical ordering aside from the Marking's hashCode.
				break;
			} else if (key < midKey) {
				end = mid - 1;
			} else {
				start = mid + 1;
			}
		}
		states.add(mid, s);
		return true;
	}

	public S getFirst() {
		return states.get(0);
	}

	public S pull() {
		return states.remove(0);
	}

	boolean isEmpty() {
		return states.isEmpty();
	}

	public Iterator<S> iterator() {
		return states.iterator();
	}

	public String toString() {
		return states.toString();
	}

	public S getFinalState() {
		return finalState;
	}

	public int getMinimalEditDistance() {
		return finalDistance;
	}

}
