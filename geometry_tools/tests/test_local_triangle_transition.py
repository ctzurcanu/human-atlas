import numpy as np
import unittest

from geometry_tools.algorithms.local_triangle_transition import check_transition
from geometry_tools.model import Mesh


def mesh(v, f):
    return Mesh(np.asarray(v, float), np.asarray(f, int), 'test:metres')


def check_new_crossing_is_rejected_and_unmoved_pairs_cancel():
    a = mesh([[-1, -1, 1], [1, -1, 1], [0, 1, 1]], [[0, 1, 2]])
    b = mesh([[-1, -1, -1], [1, -1, 1], [0, 1, 1]], [[0, 1, 2]])
    wall = mesh([[-2, -2, 0], [2, -2, 0], [0, 2, 0]], [[0, 1, 2]])
    result = check_transition(a, b, ['artery'], [('skin', wall)])
    assert result['complete'] and not result['passed']
    assert result['rows'] == [dict(owners=['artery', 'fixed:skin'], before=0, after=1)]
    assert check_transition(b, a, ['artery'], [('skin', wall)])['passed']


def check_coherent_member_transport_and_literal_fixed_context():
    a = mesh([[-1, -1, 1], [1, -1, 1], [0, 1, 1],
              [-1, -1, 3], [1, -1, 3], [0, 1, 3]], [[0, 1, 2], [3, 4, 5]])
    b = mesh(a.vertices + [0, 0, .2], a.faces)
    assert check_transition(a, b, ['artery', 'vein'])['passed']
    # A moving artery may not collide with a fixed vein.
    vertices = a.vertices.copy()
    vertices[:3] += [0, 0, 2]
    c = mesh(vertices, a.faces)
    r = check_transition(a, c, ['artery', 'vein'])
    assert not r['passed'] and r['rows'][0]['after'] == 1
    with unittest.TestCase().assertRaises(ValueError):
        check_transition(a, b, ['artery'])


class LocalTriangleTransitionTests(unittest.TestCase):
    def test_new_crossing_is_rejected_and_unmoved_pairs_cancel(self):
        check_new_crossing_is_rejected_and_unmoved_pairs_cancel()

    def test_coherent_member_transport_and_literal_fixed_context(self):
        check_coherent_member_transport_and_literal_fixed_context()
