import copy
from pathlib import Path
import tempfile
import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.common import digest_file
from geometry_tools.reference_refit import verify_reference_refit


class ReferenceRefitTests(unittest.TestCase):
    def fixture(self, directory):
        v=np.array([[0.,0,0],[.01,0,0],[.01,.01,0],[0,.01,0]])
        f=np.array([[0,1,2],[0,2,3]])
        old=Mesh(v,f,'fixture',{'uvs':v[:,:2].astype('<f4')})
        rv=np.vstack([v,old.triangles.mean(1)])
        rf=np.array([[0,1,4],[1,2,4],[2,0,4],[0,2,5],[2,3,5],[3,0,5]])
        base=Mesh(rv,rf,'fixture',{'uvs':rv[:,:2].astype('<f4')})
        base.save(directory/'source.npz')
        np.savez(directory/'lineage.npz',parent_faces=np.repeat([0,1],3),original_vertices=np.arange(4))
        np.savez(directory/'controls.npz',pinned=np.array([True]*4+[False]*2))
        policy=dict(schema='open-source-z-refit/v1',source_geometry_sha256=old.geometry_sha,
            refined_source=str(directory/'source.npz'),lineage=str(directory/'lineage.npz'),controls=str(directory/'controls.npz'),
            encoding_tolerance_m=1e-8,maximum_z_shift_m=.001)
        for key in ['refined_source','lineage','controls']:policy[key+'_sha256']=digest_file(policy[key])
        cv=rv.copy();cv[4:,2]=.0003
        return old,Mesh(cv,rf,'fixture',base.attributes),policy

    def test_local_refit_retains_source_parents_pins_and_uvs(self):
        with tempfile.TemporaryDirectory()as td:
            old,new,policy=self.fixture(Path(td))
            result=verify_reference_refit(old,new,policy)
            self.assertTrue(result['passed']);self.assertFalse(result['new_finite_tissue_claim'])
            self.assertAlmostEqual(result['maximum_z_shift_m'],.0003)

    def test_stale_source_changed_uvs_and_unpinned_reference_are_rejected(self):
        with tempfile.TemporaryDirectory()as td:
            old,new,policy=self.fixture(Path(td))
            stale=copy.deepcopy(policy);stale['source_geometry_sha256']='0'*64
            with self.assertRaises(ValueError):verify_reference_refit(old,new,stale)
            uv=new.attributes['uvs'].copy();uv[0]=1
            with self.assertRaises(ValueError):verify_reference_refit(old,Mesh(new.vertices,new.faces,new.frame,{'uvs':uv}),policy)
            v=new.vertices.copy();v[0,2]=1e-5
            with self.assertRaises(ValueError):verify_reference_refit(old,Mesh(v,new.faces,new.frame,new.attributes),policy)


if __name__=='__main__':unittest.main()
