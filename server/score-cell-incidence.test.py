"""Dependent score planes must not create tissue on their contact line."""
import unittest
import tempfile
import numpy as np
from geometry_tools.score_cells import ScoreArrangement
from geometry_tools.implicit_cells import grid_tetrahedra
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

class ScoreIncidences(unittest.TestCase):
    def test_three_competitive_regions_have_the_correct_closed_volumes(self):
        p,t,_=grid_tetrahedra([0.,0.,0.],[1.,1.,1.],1.)
        scores=np.c_[p[:,0],p[:,1],np.full(len(p),.035)]
        def classify(values):
            winner=2
            for candidate in [0,1]:
                k,sign=arrangement.comparison_index[(winner,candidate)]
                if sign*values[k]<0:winner=candidate
            return['x','y','background'][winner]
        arrangement=ScoreArrangement(p,scores,[(2,0),(2,1),(0,1)],classify,background_index=2)
        v,s,c=arrangement.build(t);volumes={}
        for name,uses in c.items():
            f=np.concatenate([s[u['surface']]if u['side']==1 else s[u['surface']][:,[0,2,1]]for u in uses])
            proof=topology(Mesh(v,f,'test:metres'));self.assertTrue(proof['closed_oriented_manifold'],name)
            volumes[name]=proof['oriented_triangle_integral_m3']
        self.assertAlmostEqual(volumes['background'],.035**2)
        self.assertAlmostEqual(volumes['x'],(1-.035**2)/2)
        self.assertAlmostEqual(volumes['y'],(1-.035**2)/2)
        self.assertAlmostEqual(sum(volumes.values()),1.)
        self.assertLessEqual(arrangement.maximum_derived_projection,4e-12)

    def test_saved_stage_replays_the_identical_completed_complex(self):
        p,t,_=grid_tetrahedra([0.,0.,0.],[1.,1.,1.],1.)
        scores=np.c_[p[:,0],p[:,1],np.full(len(p),.035)]
        def make():
            def classify(values):
                winner=2
                for candidate in [0,1]:
                    k,sign=arrangement.comparison_index[(winner,candidate)]
                    if sign*values[k]<0:winner=candidate
                return['x','y','background'][winner]
            arrangement=ScoreArrangement(p,scores,[(2,0),(2,1),(0,1)],classify,background_index=2)
            return arrangement
        complete=make().build(t);first=make();first.build(t[:3])
        with tempfile.TemporaryDirectory()as tmp:
            first.save_checkpoint(tmp,dict(inputs={'fixture':'cube'},next_tetrahedron=3))
            resumed=make();metadata=resumed.restore_checkpoint(tmp,{'fixture':'cube'})
            self.assertEqual(metadata['next_tetrahedron'],3);actual=resumed.build(t[3:])
        self.assertTrue(np.array_equal(complete[0],actual[0]));self.assertEqual(complete[2],actual[2])
        self.assertEqual(set(complete[1]),set(actual[1]))
        for name in complete[1]:self.assertTrue(np.array_equal(complete[1][name],actual[1][name]),name)

if __name__=='__main__':unittest.main()
