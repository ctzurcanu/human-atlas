"""Replay face ancestry before reusing a source section bank.

Raw control IDs may change when unused controls are removed. The directed
source triangles must still induce one consistent, injective control map.
Newly constructed faces are deliberately not assigned guessed quad partners.
"""
import numpy as np


def remap_source_pairs(source_faces, current_faces, kept_source_ids, pairs):
    old = np.asarray(source_faces, dtype=np.int64)
    new = np.asarray(current_faces, dtype=np.int64)
    kept = np.asarray(kept_source_ids, dtype=np.int64)
    pairing = np.asarray(pairs, dtype=np.int64)
    if old.ndim != 2 or old.shape[1] != 3 or new.ndim != 2 or new.shape[1] != 3:
        raise ValueError('Directed source triangles required')
    if kept.ndim != 1 or not len(kept) or len(kept) > len(new):
        raise ValueError('Explicit retained source face ancestry required')
    if kept.min() < 0 or kept.max() >= len(old) or len(np.unique(kept)) != len(kept):
        raise ValueError('Invalid or repeated source face ancestry')
    if pairing.ndim != 2 or pairing.shape[1] != 2 or not len(pairing):
        raise ValueError('Explicit source triangle pairs required')
    if pairing.min() < 0 or pairing.max() >= len(old) or np.unique(pairing).size != pairing.size:
        raise ValueError('Invalid or repeated source pair ownership')
    forward, reverse = {}, {}
    for a, b in zip(old[kept].ravel(), new[:len(kept)].ravel()):
        a, b = int(a), int(b)
        if forward.setdefault(a, b) != b or reverse.setdefault(b, a) != a:
            raise ValueError('Retained directed faces do not replay one source control map')
    lookup = np.full(len(old), -1, dtype=np.int64)
    lookup[kept] = np.arange(len(kept))
    mapped = lookup[pairing]
    mapped = mapped[np.all(mapped >= 0, axis=1)]
    return mapped, dict(retainedFaces=len(kept), constructedFaces=len(new)-len(kept),
                        originalPairs=len(pairing), retainedPairs=len(mapped),
                        controlRenumberingVerified=True, constructedQuadPairingInferred=False)
