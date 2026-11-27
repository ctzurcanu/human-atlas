import unittest
import numpy as np
from geometry_tools.algorithms.axis_interface_retriangulation import retriangulate_axis_interfaces,edge_chain
from geometry_tools.distance import triangle_pairs


class AxisInterfaceTests(unittest.TestCase):
    def test_an_internal_fold_is_retriangulated_without_moving_its_boundary(self):
        v=np.array([[0.,0.,.25],[1.,0.,.25],[1.,1.,.25],[0.,1.,.25],[1.2,.5,.25]])
        f=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
        out,r=retriangulate_axis_interfaces(v,{'inside':f},['inside'])
        self.assertEqual(len(r['changes']),1);self.assertEqual(len(out['inside']),2)
        self.assertEqual(edge_chain(out['inside']),edge_chain(f));self.assertFalse(r['coordinates_changed'])
        protected,_=retriangulate_axis_interfaces(v,{'source':f},[])
        np.testing.assert_array_equal(protected['source'],f)

    def test_retriangulation_preserves_an_actual_hole(self):
        xy=np.array([[-2.,-2.],[2.,-2.],[2.,2.],[-2.,2.],[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.],[0.,0.]])
        v=np.c_[xy,np.full(9,.5)];f=[]
        for k in range(4):
            j=(k+1)%4;f.extend([[k,j,4+j],[k,4+j,4+k]])
        t=f.pop(0);f=np.array([[t[0],t[1],8],[t[1],t[2],8],[t[2],t[0],8],*f])
        out,r=retriangulate_axis_interfaces(v,{'inside':f},['inside']);new=out['inside']
        self.assertEqual(edge_chain(new),edge_chain(f));self.assertEqual(r['changes'][0]['boundary_loops'],2)
        hole=v[np.array([[4,5,6],[4,6,7]])];a=np.repeat(v[new],len(hole),axis=0);b=np.tile(hole,(len(new),1,1))
        self.assertFalse(np.any(np.isin(triangle_pairs(a,b)['relation'],['crossing','coplanar-overlap'])))


if __name__=='__main__':unittest.main()
