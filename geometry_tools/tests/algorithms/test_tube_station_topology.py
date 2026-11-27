import unittest
import numpy as np
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence,StationTopologyError

def tube():
    n,k=12,8;a=np.arange(n)*2*np.pi/n;v=np.concatenate([np.c_[np.full(n,x),.002*np.cos(a)+.0003*np.sin(x*100),.002*np.sin(a)]for x in np.linspace(0,.02,k)]);f=[]
    for s in range(k-1):
        for j in range(n):
            x=s*n+j;y=s*n+(j+1)%n;f.extend([[x,y,y+n],[x,y+n,x+n]])
    return v,np.array(f)
class StationTests(unittest.TestCase):
    def test_scrambled_ring_cycles_recover_material_phase_from_source_faces(self):
        v,f=tube();rng=np.random.default_rng(313);p=rng.permutation(len(v));inv=np.argsort(p)
        vv=v[p];ff=inv[f[rng.permutation(len(f))]];saved_faces=ff.copy()
        group=recover_station_rings(vv,ff)[0]
        rings,proof=align_station_correspondence(vv,ff,group['rings'])
        normalized=[vv[q]-vv[q].mean(0)for q in rings]
        for q in normalized[1:]:np.testing.assert_allclose(q,normalized[0],rtol=0.,atol=1e-16)
        self.assertTrue(all(q['sourceFacesReplayedExactly']for q in proof))
        self.assertEqual(sum(len(q['sourceFaceIds'])for q in proof),len(ff))
        np.testing.assert_array_equal(vv,v[p]);np.testing.assert_array_equal(ff,saved_faces)
    def test_vertex_and_face_order_do_not_determine_sections(self):
        v,f=tube();rng=np.random.default_rng(71);permutation=rng.permutation(len(v));inverse=np.argsort(permutation);order=rng.permutation(len(f));r=recover_station_rings(v[permutation],inverse[f[order]])
        self.assertEqual(len(r),1);self.assertEqual(len(r[0]['rings']),8)
        for ring in r[0]['rings']:
            x=v[permutation][ring,0];self.assertEqual(len(x),12);self.assertEqual(np.ptp(x),0)
    def test_separate_source_tubes_keep_separate_stations(self):
        v,f=tube();r=recover_station_rings(np.r_[v,v+[0,.01,0]],np.r_[f,f+len(v)]);self.assertEqual(len(r),2);self.assertEqual([len(x['rings'])for x in r],[8,8])
    def test_missing_wall_strip_is_not_inferred_back_into_a_tube(self):
        v,f=tube()
        with self.assertRaises(StationTopologyError):recover_station_rings(v,f[2:])
    def test_closed_solid_does_not_supply_an_open_lumen(self):
        v,f=tube();points=list(v);ff=list(f)
        for s,sign in [(0,-1),(7,1)]:
            c=len(points);points.append(v[s*12:(s+1)*12].mean(0))
            for j in range(12):
                x=s*12+j;y=s*12+(j+1)%12;ff.append([y,x,c]if sign<0 else[x,y,c])
        with self.assertRaises(StationTopologyError):recover_station_rings(np.array(points),np.array(ff))
if __name__=='__main__':unittest.main()
