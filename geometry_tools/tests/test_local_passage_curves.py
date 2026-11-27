import unittest
import numpy as np
from geometry_tools.algorithms.local_passage_curves import compare_local_passage_curves


def loops(centers):
    points=[];edges=[];components=[]
    for c in centers:
        ids=list(range(len(points),len(points)+4));p=np.array([[0,0,0],[.001,0,0],[.001,.001,0],[0,.001,0]])+c
        points.extend(p.tolist());edges.extend([[ids[i],ids[(i+1)%4]]for i in range(4)]);components.append(dict(controls=ids,closedUnbranched=True,centroid_m=p.mean(0).tolist(),length_m=.004))
    return dict(points_m=points,edges=edges,components=components)


class LocalPassageTests(unittest.TestCase):
    def test_small_existing_passage_shift(self):
        r=compare_local_passage_curves(loops([[0,0,0]]),loops([[.0004,0,0]]),maximum_displacement_m=.001)
        self.assertTrue(r['passed']);self.assertLess(r['matched'][0]['hausdorffUpperBound_m'],.0005)

    def test_additional_loop_is_rejected(self):
        self.assertFalse(compare_local_passage_curves(loops([[0,0,0]]),loops([[0,0,0],[.0001,0,0]]),maximum_displacement_m=.001)['passed'])

    def test_moved_to_another_patch_is_rejected(self):
        self.assertFalse(compare_local_passage_curves(loops([[0,0,0]]),loops([[.01,0,0]]),maximum_displacement_m=.001)['passed'])

    def test_removed_unaccepted_contact_is_recorded(self):
        r=compare_local_passage_curves(loops([[0,0,0]]),loops([]),maximum_displacement_m=.001)
        self.assertTrue(r['passed']);self.assertEqual(r['removedOldLoops'],[0])


if __name__=='__main__':unittest.main()
