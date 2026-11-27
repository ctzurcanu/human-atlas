"""Load explicit section-pair evidence only when its directed faces replay."""
import numpy as np
from geometry_tools.course_io import sha


def load_source_pairing(binding, faces):
    if not binding:
        raise ValueError('Missing explicit source triangle-pair bank; rebind source profiles first')
    if sha(binding['path']) != binding['sha256']:
        raise ValueError('Source section bank changed')
    with np.load(binding['path'], allow_pickle=False) as z:
        if 'source_faces' not in z or not np.array_equal(z['source_faces'], faces):
            raise ValueError('Stale source section face ancestry; rebind source profiles first')
        return z[binding.get('array', 'retained_source_triangle_pairs')].copy()
