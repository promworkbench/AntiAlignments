package org.processmining.antialignments.base;

import gnu.trove.map.hash.TObjectByteHashMap;

public class KeyLookupHashMap<S> extends TObjectByteHashMap<S> {

	public KeyLookupHashMap(int size) {
		super(size);
	}

	@SuppressWarnings("unchecked")
	public S getKeyIfPresent(S key) {
		int i = index(key);
		if (i == -1) {
			return null;
		} else {
			return (S) _set[i];
		}
	}
}