"""Merge disjoint source-region edits into one literal shared attribute bank."""
import numpy as np


def merge_source_updates(source_arrays, updates):
    arrays = {key: np.array(value, copy=True) for key, value in source_arrays.items()}
    written = {key: np.zeros(len(value), dtype=bool) for key, value in arrays.items()}
    for ids, values in updates:
        ids = np.asarray(ids, dtype=int)
        if ids.ndim != 1 or len(np.unique(ids)) != len(ids):
            raise ValueError('Source control IDs must be unique and one dimensional')
        for key, value in values.items():
            target = arrays[key]
            value = np.asarray(value)
            if len(ids) and (ids.min() < 0 or ids.max() >= len(target)):
                raise ValueError('Source control ID outside the shared bank')
            if value.dtype != target.dtype or value.shape != target[ids].shape:
                raise ValueError('Shared update must preserve source dtype and width')
            seen = written[key][ids]
            if not np.array_equal(target[ids[seen]], value[seen]):
                raise ValueError('Conflicting updates on a shared source control')
            target[ids] = value
            written[key][ids] = True
    return arrays
