import unittest
import manifold3d as md
import numpy as np

from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.tube_rim_bridge import bridge_boundaries, end_annulus
from geometry_tools.algorithms.constrained_boundary_disk import constrained_boundary_disk
from geometry_tools.audit import audit
from geometry_tools.model import Mesh


class RimBridgeTests(unittest.TestCase):
    def test_flared_profile_cannot_pinch_a_cross_edge(self):
        a=np.arange(12)*2*np.pi/12;b=np.arange(7)*2*np.pi/7
        p=np.r_[np.c_[.001*np.cos(a),.001*np.sin(a),np.zeros(12)],
                np.c_[.01*np.cos(b)+.006,.01*np.sin(b),np.full(7,.001)]]
        f,_=bridge_boundaries(p,np.arange(12)[::-1],np.arange(12,19))
        edges=np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
        unique,count=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
        self.assertLessEqual(count.max(),2);self.assertEqual(len(np.unique(f))-len(unique)+len(f),0)

    def test_computational_disk_cannot_reuse_a_forbidden_source_diagonal(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
                    [2.,2.,1.],[3.,2.,1.],[2.,3.,1.]])
        f,_=constrained_boundary_disk(p,[0,1,2,3],np.array([[4,5,6]]),forbidden_edges=[(0,2)])
        edges={tuple(sorted((int(a),int(b))))for t in f for a,b in zip(t,np.roll(t,-1))}
        self.assertNotIn((0,2),edges);self.assertIn((1,3),edges)

    def test_nonplanar_unequal_pipe_rims_join_without_a_lumen_disk(self):
        pieces=[];boundaries=[];points=[];offset=0
        for n,z,upper in [(8,-.02,True),(13,.015,False)]:
            solid=md.Manifold.cylinder(.02,.01,.01,n)-md.Manifold.cylinder(.02,.008,.008,n)
            mesh=solid.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3].copy();f=np.asarray(mesh.tri_verts,int)
            v[:,2]+=z
            if not upper:v[:,2]+=.0005*np.sin(3*np.arctan2(v[:,1],v[:,0]))
            radius=np.linalg.norm(v[:,:2],axis=1)
            inner=f[np.all(abs(radius[f]-.008)<1e-12,axis=1)]
            mouths=boundary_loops(inner)
            mouth=(max if upper else min)(mouths,key=lambda q:v[q,2].mean())
            annulus,outer=end_annulus(f,inner,mouth)
            pieces.append(f[np.setdiff1d(np.arange(len(f)),annulus)]+offset)
            boundaries.append((outer+offset,mouth+offset));points.append(v);offset+=len(v)
        p=np.concatenate(points);literal=p.copy()
        for k in [0,1]:
            f,proof=bridge_boundaries(p,boundaries[0][k],boundaries[1][k])
            self.assertFalse(proof['lumen_disk_created']);self.assertTrue(proof['exact_opposed_boundary_chain'])
            pieces.append(f)
        result=audit(Mesh(p,np.concatenate(pieces),'fixture'),require_closed=True,expected_components=1,expected_genus=[1])
        self.assertTrue(result['passed'],result['failures']);self.assertTrue(np.array_equal(p,literal))


if __name__=='__main__':unittest.main()
