import json,tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.fit import landmark_transform,transform_mesh,constrained_translation,cage_mesh
from geometry_tools.tissue import layer_arrangement,pack,validate_cells
from geometry_tools.tubes import audit_fields,join,sweep,path
from geometry_tools.constraints import check_constraint
from geometry_tools.install import transaction,restore_transaction,catalogue_diff
from geometry_tools.common import digest_file
from test_core import tetra
from geometry_tools.algorithms.tube_continuity import loft_faces


def tube_fields(start=0.,stop=.1,stations=5):
    a=np.arange(12)*2*np.pi/12;c=np.c_[np.zeros(stations),np.zeros(stations),np.linspace(start,stop,stations)]
    q=np.array([c[:,None]+r*np.c_[np.cos(a),np.sin(a),np.zeros(12)][None]for r in [.01,.008,.006]])
    return q.reshape(3,-1,3),loft_faces(stations,12),c


class ConstructionTests(unittest.TestCase):
    def test_translation_keeps_literal_source_normal_words_and_uvs(self):
        m=tetra();words=np.tile([1000,-2000,32000],(4,1)).astype('<i2')
        m=Mesh(m.vertices,m.faces,m.frame,{'normals':words,'uvs':np.arange(8,dtype='<f4').reshape(4,2)})
        matrix=np.eye(4);matrix[:3,3]=[.013,0,-.004]
        q=transform_mesh(m,matrix)
        np.testing.assert_array_equal(q.attributes['normals'],words)
        np.testing.assert_array_equal(q.attributes['uvs'],m.attributes['uvs'])
        np.testing.assert_array_equal(q.faces,m.faces)
        np.testing.assert_array_equal(q.vertices,m.vertices+matrix[:3,3])

    def test_proper_similarity_and_shared_attributes(self):
        m=tetra();m=Mesh(m.vertices,m.faces,m.frame,{'normals':np.tile([0,0,32767],(4,1)).astype('<i2'),'uvs':np.zeros((4,2),'<f4'),'normalSigns':np.array([1,-1,1,1],'i1')})
        rotation=np.array([[0,-1.,0],[1,0,0],[0,0,1]]);target=.8*m.vertices@rotation.T+[2,-1,4]
        matrix,r=landmark_transform(m.vertices,target,similarity=True);np.testing.assert_allclose(transform_mesh(m,matrix).vertices,target,atol=2e-15);self.assertAlmostEqual(r['scale'],.8)
        out=transform_mesh(m,matrix);np.testing.assert_array_equal(out.attributes['uvs'],m.attributes['uvs']);np.testing.assert_array_equal(out.attributes['normalSigns'],m.attributes['normalSigns'])
        with self.assertRaises(ValueError):transform_mesh(m,np.diag([1.,1,-1,1]))

    def test_translation_feasibility_does_not_shrink_geometry(self):
        matrix,r=constrained_translation([0,0,0],[{'normal':[0,1,0],'minimum_dot_m':.01}]);self.assertAlmostEqual(matrix[1,3],.01)
        with self.assertRaises(ValueError):constrained_translation([0,0,0],[{'normal':[1,0,0],'minimum_dot_m':1},{'normal':[-1,0,0],'minimum_dot_m':0}])

    def test_cage_has_positive_jacobian_and_does_not_change_faces(self):
        m=tetra();spec=dict(source_bounds=[0,1],target_bounds=[0,1.2],fractions=[0,1],source_centers=[[0,0],[0,0]],target_centers=[[.1,0],[.2,0]],source_covariances=[[[1,0],[0,.5]]]*2,target_covariances=[[[1.1,0],[0,.6]]]*2)
        out,r=cage_mesh(m,spec);self.assertGreater(r['minimum_sampled_jacobian'],0);np.testing.assert_array_equal(out.faces,m.faces)

    def test_encoded_cells_close_only_wall_and_share_normals(self):
        fields,faces,c=tube_fields();positions,surfaces,cells,finite=layer_arrangement(fields,faces,['outer','inner'])
        with tempfile.TemporaryDirectory()as td:
            r=pack(positions,surfaces,cells,{'outer':'outer','inner':'inner'},'fixture',Path(td)/'tissue',construction=finite)
            self.assertTrue(r['passed']);self.assertTrue(r['positions_encoded_once']);self.assertGreater(r['checks']['cell_volumes_m3']['inner'],0)
            # Every annular-rim point stays at radius >= the inner ring.
            rim=np.concatenate([positions[f].reshape(-1,3)for sid,f in surfaces.items()if sid.startswith('rim:')]);self.assertGreaterEqual(np.linalg.norm(rim[:,:2],axis=1).min(),.006-1e-15)
            rows=r['parts'];raw=Path(r['bank']).read_bytes();a,b=rows['outer'],rows['inner'];self.assertEqual(a['positions'],b['positions']);self.assertEqual(a['normals'],b['normals'])
            sa=np.ndarray((a['vertexCount'],),'i1',buffer=raw,offset=a['normalSigns']);sb=np.ndarray((b['vertexCount'],),'i1',buffer=raw,offset=b['normalSigns'])
            shared=r['checks']['cells']['outer']['topology']['closed_oriented_manifold'];self.assertTrue(shared)
            self.assertTrue(any(sa!=sb))

    def test_volume_overlap_is_rejected(self):
        r=validate_cells({'a':tetra(),'b':tetra((.1,.1,.1))});self.assertFalse(r['passed']);self.assertTrue(any(x.startswith('volume_overlap')for x in r['failures']))

    def test_inverted_wall_and_obstructed_guide_are_rejected(self):
        fields,faces,c=tube_fields();good=audit_fields(fields,faces,'fixture',centres=c);self.assertTrue(good['passed'])
        self.assertFalse(audit_fields(fields[::-1],faces,'fixture')['passed'])
        r=audit_fields(fields,faces,'fixture',centres=np.array([[0,0,.04],[.02,0,.06]]));self.assertFalse(r['lumen_guide']['passed'])

    def test_join_uses_shared_open_mouth_not_plug(self):
        a,fa,ca=tube_fields(0,.1);b,fb,cb=tube_fields(.1,.2)
        with tempfile.TemporaryDirectory()as td:
            td=Path(td);np.savez(td/'a.npz',fields=a,faces=fa,frame=np.array('fixture'));np.savez(td/'b.npz',fields=b,faces=fb,frame=np.array('fixture'))
            r=join(dict(first=str(td/'a.npz'),second=str(td/'b.npz'),first_mouth=list(range(48,60)),second_mouth=list(range(12)),position_tolerance_m=1e-15),td/'joined.npz')
            self.assertTrue(r['passed']);self.assertFalse(r['central_lumen_capped']);self.assertEqual(r['checks']['carriers'][0]['mouths'],2)

    def test_rounded_sweep_replays_emitted_area_and_wall(self):
        fields,faces,c=tube_fields(0,.18,25)
        with tempfile.TemporaryDirectory()as td:
            td=Path(td);source=td/'source.npz';guide=td/'path.npz';target=td/'bent.npz';np.savez(source,fields=fields,faces=faces,frame=np.array('fixture'))
            r=path({'frame':'fixture','waypoints':[[0,0,0],[0,0,.09],[.09,0,.09]],'bend_radius_m':.04,'spacing_m':.002},guide);self.assertTrue(r['passed'])
            r=sweep({'source':str(source),'path':str(guide),'ring_size':12},target)
            self.assertTrue(r['passed']);self.assertLess(r['checks']['encoded_area_change_relative'],1e-5)
            with self.assertRaises(FileExistsError):sweep({'source':str(source),'path':str(guide),'ring_size':12},target)

    def test_mesh_constraints_require_geometry_and_keep_frames_explicit(self):
        a=tetra();b=tetra((2,0,0));r=check_constraint({'kind':'clearance','parameters':{'minimum_m':.5}},a,b);self.assertTrue(r['passed'])
        self.assertFalse(check_constraint({'kind':'fixed','parameters':{'geometry_sha256':a.geometry_sha}},b)['passed'])
        with self.assertRaises(ValueError):check_constraint({'kind':'clearance'},a,Mesh(b.vertices,b.faces,'other'))
        inner=Mesh(a.vertices*.1+.1,a.faces,a.frame);self.assertTrue(check_constraint({'kind':'nested','parameters':{}},inner,a)['passed'])
        self.assertFalse(check_constraint({'kind':'disjoint','parameters':{'occupancy':'solid'}},inner,a)['passed'])

    def test_conduit_constraint_rejects_a_disk_inside_matched_mouths(self):
        af,f,_=tube_fields(0,.1);bf,g,_=tube_fields(.1,.2);a=Mesh(af[-1],f,'fixture');b=Mesh(bf[-1],g,'fixture')
        spec={'kind':'open-conduit','parameters':{'first_ring':list(range(48,60)),'second_ring':list(range(12)),'position_tolerance_m':1e-15}}
        self.assertTrue(check_constraint(spec,a,b)['passed'])
        v=np.r_[b.vertices,b.vertices[:12].mean(0)[None]];plug=np.array([[len(b.vertices),i,(i+1)%12]for i in range(12)]);blocked=Mesh(v,np.r_[g,plug],'fixture')
        self.assertFalse(check_constraint(spec,a,blocked)['passed'])

    def test_shared_interface_checks_actual_normal_words(self):
        v=np.array([[0.,0,0],[1,0,0],[0,1,0]]);a=Mesh(v,np.array([[0,1,2]]),'fixture',{'normals':np.tile([0,0,32767],(3,1)).astype('<i2')})
        b=Mesh(v,np.array([[0,2,1]]),'fixture',{'normals':-a.attributes['normals']});spec={'kind':'shared-interface','parameters':{'first_faces':[0],'second_faces':[0]}}
        self.assertTrue(check_constraint(spec,a,b)['passed']);wrong=Mesh(v,b.faces,'fixture',a.attributes);self.assertFalse(check_constraint(spec,a,wrong)['passed'])

    def test_failed_file_transaction_restores_literal_previous_bytes(self):
        with tempfile.TemporaryDirectory()as td:
            td=Path(td);a=td/'a.json';b=td/'b.json';a.write_bytes(b'original');expected={str(a):digest_file(a),str(b):None}
            with self.assertRaises(RuntimeError):transaction({a:b'candidate',b:b'new'},td/'backup',expected=expected,after_write=lambda:(_ for _ in ()).throw(RuntimeError('failed derived rebuild')))
            self.assertEqual(a.read_bytes(),b'original');self.assertFalse(b.exists());self.assertTrue(list((td/'backup').glob('retained-*')))
            with self.assertRaises(ValueError):transaction({a:b'wrong'},td/'other',expected={str(a):'bad'})

    def test_rollback_refuses_later_user_edit(self):
        with tempfile.TemporaryDirectory()as td:
            td=Path(td);a=td/'a';a.write_bytes(b'original');r=transaction({a:b'candidate'},td/'backup',expected={str(a):digest_file(a)})
            a.write_bytes(b'user edit')
            with self.assertRaises(ValueError):restore_transaction(r['manifest'])
            self.assertEqual(a.read_bytes(),b'user edit')

    def test_catalogue_replacement_cannot_move_skeleton_or_other_parts(self):
        a={'parts':[{'id':'bone','system':'skeletal'},{'id':'organ','system':'digestive'}],'chunks':[]};b=json.loads(json.dumps(a));b['parts'][0]['positions']=0
        with self.assertRaises(ValueError):catalogue_diff(a,b,['bone'])
        with self.assertRaises(ValueError):catalogue_diff(a,b,['organ'])


if __name__=='__main__':unittest.main()
