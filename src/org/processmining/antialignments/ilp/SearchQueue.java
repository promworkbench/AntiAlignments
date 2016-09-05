package org.processmining.antialignments.ilp;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.collections15.list.TreeList;
import org.processmining.antialignments.bruteforce.State;

/**
 * This class implements the search queue for an anti alignment.
 * 
 * Internally, states are ordered by their marking, in order to allow for each
 * new state to be quickly compared to existing states (binary search on the
 * marking, followed by a linear search through identical markings)
 * 
 * On the side, an iteration order is kept for iterating based on the state's
 * "getPathDistance" value. This iteration order is kept consistent under
 * updates which is a linear operation.
 * 
 * @author bfvdonge
 * 
 * @param <S>
 */
public class SearchQueue<S extends State> implements Iterable<S> {

	private static final int INITIALSIZE = 1024;
	//For now, just use a LinkedList of states, sorted by Marking.
	private TreeList<S> states = new TreeList<>();
	// nextElement stores the relative index of the next element in an iteration order
	// by state.getPathDistance. The next element to iterate at index i, is "(i+nextElement[i])%states.size();"
	private int startAt;
	private int[] nextElement;
	private int[] pathDistances;

	public SearchQueue() {
		// initialize the nextElement array to some size;
		nextElement = new int[INITIALSIZE];
		pathDistances = new int[INITIALSIZE];
		startAt = -1;
		added = 0;
		removed = 0;
	}

	public boolean add(S s) {

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
							// Nothing changes, no need to update the iteration.
							return false;
						} else if (s.isGreaterOrEqual(stateAtM)) {
							// an existing state is strictly smaller than the new one. Remove!
							//							System.out.println("Removing " + stateAtM + " because of " + s);
							end--;
							states.remove(m);
							// removing state at index m, hence the nextElement array has to be updated.
							updateAfterRemoval(m);
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
		updateAfterAdding(mid, s.getPathDistance());
		return true;
	}

	private static int removed = 0;

	private void updateAfterRemoval(int removedAt) {
		removed++;
		int oldSize = states.size() + 1;
		int oldRemoveAtPointer = nextElement[removedAt];
		int oldStartAt = startAt;
		if (startAt == removedAt) {
			startAt = (removedAt + oldRemoveAtPointer) % oldSize;
			if (startAt >= removedAt) {
				startAt--;
			}
		} else if (startAt > removedAt) {
			startAt--;
		}
		for (int i = 0; i < removedAt; i++) {
			if (i + nextElement[i] > removedAt) {
				// points passed the removed element, decrease by 1
				nextElement[i]--;
			} else if (i + nextElement[i] == removedAt) {
				// points to the removed element, take new pointer
				nextElement[i] = (nextElement[i] + oldRemoveAtPointer - 1) % states.size();
				if (i + nextElement[i] > states.size() && pathDistances[i + nextElement[i]] == pathDistances[i]) {
					// swap [i] and [i+nextElement[i]%states.size()]
					int j = (i + nextElement[i]) % states.size();
					assert j < i;
					int tmp = nextElement[j];
					nextElement[j] = i - j;
					nextElement[i] = tmp - nextElement[j];
				}
			}
		}
		for (int i = removedAt + 1; i < oldSize; i++) {
			if (i + nextElement[i] > oldSize + removedAt) {
				// points passed the removed element, decrease by 1 and move up one index
				nextElement[i - 1] = nextElement[i] - 1;
			} else if (i + nextElement[i] == oldSize + removedAt) {
				// points to the removed element, take new pointer
				if (oldRemoveAtPointer == 0) {
					// we removed the last element in the queue
					nextElement[i - 1] = 0;
				} else {
					// we removed an intermediate element
					nextElement[i - 1] = oldRemoveAtPointer - 1;
				}
			} else {
				// simply move up one index
				nextElement[i - 1] = nextElement[i];
			}
			pathDistances[i - 1] = pathDistances[i];
		}
		// possibly shrink the array, making sure it has 33% free space after shrinking
		if (states.size() < nextElement.length / 3 && nextElement.length > INITIALSIZE) {
			// shrink array.
			int size = nextElement.length / 2;
			if (size < INITIALSIZE) {
				size = INITIALSIZE;
			}
			int[] newNextElement = new int[size];
			System.arraycopy(nextElement, 0, newNextElement, 0, states.size());
			nextElement = newNextElement;
			int[] newPathDistances = new int[size];
			System.arraycopy(pathDistances, 0, newPathDistances, 0, states.size());
			pathDistances = newPathDistances;
		}
		assert isConsistent();

	}

	static int added = 0;

	private void updateAfterAdding(int addedAt, int addedPathDistance) {
		added++;
		int oldSize = states.size() - 1;
		if (states.size() > nextElement.length) {
			// grow array by 50%
			int[] newNextElement = new int[(nextElement.length * 3) / 2];
			System.arraycopy(nextElement, 0, newNextElement, 0, nextElement.length);
			nextElement = newNextElement;
			int[] newPathDistances = new int[(pathDistances.length * 3) / 2];
			System.arraycopy(pathDistances, 0, newPathDistances, 0, pathDistances.length);
			pathDistances = newPathDistances;
		}

		int predecessor = -1;
		int nextPre = 0;
		for (int i = addedAt; i-- > 0;) {
			if (pathDistances[i] <= addedPathDistance
					&& (predecessor == -1 || pathDistances[i] > pathDistances[predecessor])) {
				// found closest predecessor with same path distance.
				// points to the removed element, take new pointer
				predecessor = i;
				nextPre = nextElement[i];
			}
			if (i + nextElement[i] >= addedAt) {
				// points passed the added element, increase by 1
				nextElement[i]++;
			}
		}
		for (int i = states.size(); i-- > addedAt + 1;) {
			if (i - 1 + nextElement[i - 1] >= oldSize + addedAt) {
				// points passed the removed element, increase by 1 and move down up index
				nextElement[i] = nextElement[i - 1] + 1;
			} else {
				// simply move up one index
				nextElement[i] = nextElement[i - 1];
			}
			pathDistances[i] = pathDistances[i - 1];
			if (pathDistances[i] < addedPathDistance && i > predecessor
					&& (predecessor == -1 || pathDistances[i] >= pathDistances[predecessor])) {
				// found closest predecessor with same path distance.
				// points to the removed element, take new pointer
				predecessor = i;
				nextPre = nextElement[i];
			}
		}
		// predecessor becomes my predecessor.
		pathDistances[addedAt] = addedPathDistance;
		if (startAt >= addedAt) {
			startAt++;
		}
		if (predecessor == -1) {
			// smallest element, i.e. head of the queue
			if (startAt == -1) {
				// size is now 1.
				nextElement[addedAt] = 0;
				startAt = 0;
			} else {
				// size is now > 1
				nextElement[addedAt] = startAt < addedAt ? states.size() - (addedAt - startAt) : startAt - addedAt;
				startAt = addedAt;
			}
		} else {
			if (nextPre == 0) {
				// predecessor was last element in the queue
				nextElement[addedAt] = 0;
				// now point the predecessor to addedAt
				nextElement[predecessor] = (states.size() + (addedAt - predecessor)) % states.size();
			} else {
				// point addedAt to nextPre.
				nextElement[addedAt] = (states.size() + (predecessor + nextElement[predecessor]) - addedAt)
						% states.size();
				// now point the predecessor to addedAt
				nextElement[predecessor] = (states.size() + (addedAt - predecessor)) % states.size();
			}

		}
		assert isConsistent();
	}

	public S getFirst() {
		return states.get(startAt);
	}

	public S pull() {
		S state = states.remove(startAt);
		updateAfterRemoval(startAt);
		return state;
	}

	public boolean isEmpty() {
		return states.isEmpty();
	}

	public Iterator<S> iterator() {
		return states.iterator();
	}

	public String toString() {
		return states.toString();
	}

	private static class SortedIterator<S extends State> implements Iterator<S> {

		private SearchQueue<S> queue;
		private int next;

		public SortedIterator(SearchQueue<S> queue) {
			this.queue = queue;
			this.next = queue.startAt;
		}

		public boolean hasNext() {
			return next >= 0;
		}

		public S next() {
			int oldNext = next;
			if (queue.nextElement[next] != 0) {
				assert queue.nextElement[next] > 0;
				next = ((next + queue.nextElement[next]) % queue.states.size());
			} else {
				next = -1;
			}
			return queue.states.get(oldNext);
		}

		public void remove() {
			throw new UnsupportedOperationException("Cannot remove from SearchQueue trough iterator");
		}

	}

	private boolean isConsistent() {
		boolean c = true;
		if (states.isEmpty()) {
			c &= startAt == -1;
			assert c;
		} else {
			c &= startAt >= 0;
			assert c;
			c &= startAt < states.size();
			assert c;
		}
		boolean oneZero = false;
		for (int i = 0; i < states.size(); i++) {
			c &= nextElement[i] > 0 || (!oneZero && nextElement[i] == 0);
			oneZero = nextElement[i] == 0;
			assert c;
		}
		for (int i = 0; i < states.size(); i++) {
			if (i + nextElement[i] < states.size()) {
				c &= pathDistances[i + nextElement[i]] >= pathDistances[i];
				assert c;
			} else {
				c &= pathDistances[i + nextElement[i] - states.size()] > pathDistances[i];
				assert c;
			}
		}
		Set<S> seen = new HashSet<S>();
		Iterator<S> it = iterator();
		while (c && it.hasNext()) {
			c &= seen.add(it.next());
			assert c;
		}
		c &= seen.size() == states.size();
		assert c;
		return c;

	}

}