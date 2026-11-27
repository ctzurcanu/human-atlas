import unittest
from collections import Counter
import numpy as np
from geometry_tools.algorithms.quantized_surface_topology import conform_stored_surfaces
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces,contract_quantized_intersection_edges

class QuantizedSurfaceTopologyTests(unittest.TestCase):
    def test_decoded_chart_precision_is_not_destroyed_by_another_world_float32_cast(self):
        v=np.array([[0.,1.,0],[1e-9,1.,0],[0,1.+1e-9,0],[0,1.,1e-9]])
        stored=np.r_[v,v[[0]]];source=stored.copy();source[4,0]+=1e-12;before=stored.copy()
        surfaces={'first':np.array([[0,2,4],[4,2,1]]),
                  'second':np.array([[0,4,3],[4,1,3],[0,3,2],[1,2,3]])}
        out,_,proof=conform_stored_surfaces(stored,surfaces,source_positions=source)
        np.testing.assert_array_equal(stored,before)
        f=np.r_[out['first'],out['second']];q=stored[f]
        self.assertEqual(len(f),4)
        self.assertTrue(np.all(np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0))
        self.assertFalse(proof['storedCoordinatesChanged'])
        self.assertEqual(proof['zeroAreaFacesRemoved'],{'first':1,'second':1})

    def test_shared_tiny_edge_collapses_conformingly_without_moving_stored_points(self):
        v=np.array([[.04,.07,.03],[.041,.07,.03],[.04,.071,.03],[.04,.07,.031]])
        v=np.r_[v,(v[0]+1e-8*(v[1]-v[0]))[None]]
        surfaces={'first':np.array([[0,2,4],[4,2,1]]),'second':np.array([[0,4,3],[4,1,3],[0,3,2],[1,2,3]])}
        out,parents,proof=conform_quantized_surfaces(v,surfaces,[0,0,0]);f=np.r_[out['first'],out['second']]
        edges=Counter(tuple(sorted((int(a),int(b))))for t in f for a,b in zip(t,np.roll(t,-1)))
        self.assertTrue(all(n==2 for n in edges.values()));self.assertEqual(len(f),4)
        self.assertEqual(proof['zeroAreaFacesRemoved'],{'first':1,'second':1});self.assertFalse(proof['storedCoordinatesChanged']);self.assertLess(proof['maximumEquivalentSourcePointSeparationMeters'],2e-11)
        self.assertTrue(np.array_equal(v.astype('f4')[out['first']],v.astype('f4')[surfaces['first'][parents['first']]]))
    def test_whole_vanishing_surface_is_rejected(self):
        v=np.array([[1.,1.,1.],[1.+1e-10,1.,1.],[1.,1.+1e-10,1.]])
        with self.assertRaises(ValueError):conform_quantized_surfaces(v,{'tiny':np.array([[0,1,2]])},[0,0,0])
    def test_only_new_intersection_control_can_contract(self):
        p=np.array([-.027999998561722837,1.1449999982719425,.03800000013332115]);q=np.array([-.027799998561722838,1.1446999982719426,.03780000013332115]);r=np.array([-.027999998050680038,1.1449999964518582,.037999999531041315])
        v=np.array([p,q,p+[0,.0005,0],p+[0,0,.0005],r]);before=v.copy();f=np.array([[0,1,4],[1,2,4],[2,0,4],[0,3,1],[0,2,3],[1,3,2]])
        out,_,proof=contract_quantized_intersection_edges(v,{'shell':f},[0,1.1,0],4,5)
        self.assertTrue(np.array_equal(v,before));self.assertNotIn(4,out['shell']);self.assertEqual(len(out['shell']),4);self.assertFalse(proof['literalSourceControlsMoved']);self.assertLess(proof['maximumSourceEdgeMeters'],3e-9)
        untouched,_,protected=contract_quantized_intersection_edges(v,{'shell':f},[0,1.1,0],5,5)
        self.assertTrue(np.array_equal(untouched['shell'],f));self.assertEqual(protected['intersectionOnlyEdgeContractions'],[])
if __name__=='__main__':unittest.main()
