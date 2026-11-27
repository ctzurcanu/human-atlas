import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import Mesh
from geometry_tools.reference_refit import verify_reference_refit


class ReferenceApertureTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        vertices=np.array([[x,y,0.] for y in range(4) for x in range(4)])
        faces=[]
        for y in range(3):
            for x in range(3):
                a=y*4+x;b=a+1;c=a+4;d=c+1
                faces.extend([[a,b,d],[a,d,c]])
        faces=np.array(faces)
        attributes={'normals':np.tile([0,0,32767],(16,1)).astype('<i2'),
                    'uvs':vertices[:,:2].astype('<f4')}
        self.original=Mesh(vertices,faces,'fixture:metres',attributes)
        self.removed=np.array([8,9])
        self.keep=np.ones(len(faces),bool);self.keep[self.removed]=False
        self.candidate=Mesh(vertices,faces[self.keep],self.original.frame,attributes)
        outer=np.array([[1,1,0],[2,1,0],[2,2,0],[1,2,0]],float)
        inner=np.array([[1.2,1.2,0],[1.8,1.2,0],[1.8,1.8,0],[1.2,1.8,0]],float)
        v=np.r_[outer,inner,outer-[0,0,.1],inner-[0,0,.1]]
        f=[]
        for k in range(4):
            j=(k+1)%4
            f.extend([[k,j,4+j],[k,4+j,4+k],
                      [8+k,12+j,8+j],[8+k,12+k,12+j],
                      [k,8+k,8+j],[k,8+j,j],
                      [4+k,12+j,12+k],[4+k,4+j,12+j]])
        self.collar=Mesh(v,np.array(f),self.original.frame)
        self.collar_path=self.root/'collar.npz';self.collar.save(self.collar_path)
        self.evidence=self.root/'evidence.json'
        self.evidence.write_text(json.dumps({'schema':'open-source-aperture-evidence/v1',
            'source_geometry_sha256':self.original.geometry_sha,
            'removed_source_faces':self.removed.tolist(),
            'collar_geometry_sha256':self.collar.geometry_sha}))
        self.policy={'schema':'open-source-aperture/v1','source_geometry_sha256':self.original.geometry_sha,
            'evidence':str(self.evidence),'evidence_sha256':digest_file(self.evidence),
            'collar':str(self.collar_path),'collar_sha256':digest_file(self.collar_path)}

    def test_literal_window_preserves_old_mouth_and_finite_annular_collar(self):
        result=verify_reference_refit(self.original,self.candidate,self.policy)
        self.assertEqual(result['old_source_mouths_preserved'],1)
        self.assertEqual(result['new_source_boundaries'],1)
        self.assertTrue(result['remaining_parent_is_open_reference'])
        self.assertFalse(result['new_whole_parent_finite_tissue_claim'])

    def test_wrong_mask_or_stale_source_is_rejected(self):
        bad=copy.deepcopy(self.policy);bad['source_geometry_sha256']='stale'
        with self.assertRaisesRegex(ValueError,'source or frame changed'):
            verify_reference_refit(self.original,self.candidate,bad)
        extra=Mesh(self.candidate.vertices,self.candidate.faces[:-1],self.original.frame,self.candidate.attributes)
        with self.assertRaisesRegex(ValueError,'outside the aperture'):
            verify_reference_refit(self.original,extra,self.policy)

    def test_a_source_attribute_change_is_not_an_aperture(self):
        attrs={k:v.copy() for k,v in self.candidate.attributes.items()};attrs['uvs'][0,0]=.25
        bad=Mesh(self.candidate.vertices,self.candidate.faces,self.original.frame,attrs)
        with self.assertRaisesRegex(ValueError,'attribute corners'):
            verify_reference_refit(self.original,bad,self.policy)

    def test_solid_cap_cannot_replace_annular_collar(self):
        from geometry_tools.algorithms.communicating_vascular_network import tube_classifier
        v,f,_=tube_classifier([[1.5,1.5,0],[1.5,1.5,-.1]],.4,ring_size=12,extend_joints=False)
        cap=Mesh(v,f,self.original.frame);p=self.root/'solid-cap.npz';cap.save(p)
        evidence=json.loads(self.evidence.read_text());evidence['collar_geometry_sha256']=cap.geometry_sha
        e=self.root/'cap-evidence.json';e.write_text(json.dumps(evidence))
        bad={**self.policy,'collar':str(p),'collar_sha256':digest_file(p),
             'evidence':str(e),'evidence_sha256':digest_file(e)}
        with self.assertRaisesRegex(ValueError,'finite annular tissue'):
            verify_reference_refit(self.original,self.candidate,bad)

    def test_same_boundary_points_without_the_literal_seam_edges_is_rejected(self):
        from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank
        from geometry_tools.algorithms.shared_triangle_edges import refine_surfaces
        bank=CommonPositionBank(self.collar.vertices,tolerance=0.)
        plans={}
        for a in range(4):
            b=(a+1)%4
            midpoint=bank.insert((self.collar.vertices[a]+self.collar.vertices[b])/2)
            plans[tuple(sorted((a,b)))]={midpoint}
        surfaces,_,_=refine_surfaces(bank,{'wall':self.collar.faces},plans)
        remeshed=Mesh(bank.positions,surfaces['wall'],self.original.frame)
        from geometry_tools.audit import topology
        self.assertTrue(topology(remeshed)['closed_oriented_manifold'])
        self.assertTrue(all(any(np.array_equal(p,q) for q in remeshed.vertices) for p in self.collar.vertices))
        p=self.root/'wrong-seam.npz';remeshed.save(p)
        evidence=json.loads(self.evidence.read_text());evidence['collar_geometry_sha256']=remeshed.geometry_sha
        e=self.root/'wrong-seam-evidence.json';e.write_text(json.dumps(evidence))
        bad={**self.policy,'collar':str(p),'collar_sha256':digest_file(p),
             'evidence':str(e),'evidence_sha256':digest_file(e)}
        with self.assertRaisesRegex(ValueError,'literal shared collar edge'):
            verify_reference_refit(self.original,self.candidate,bad)


if __name__=='__main__':
    unittest.main()
