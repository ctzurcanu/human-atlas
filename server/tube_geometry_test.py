"""Analytic geometry fixtures; no anatomical acceptance from synthetic tubes."""
from pathlib import Path
import sys, unittest, tempfile, subprocess, json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry.tube_sections import (TubeGeometryError,TubeTolerances,rotation_minimizing_frames,
    filleted_path,centripetal_curve,arclength,section_metrics,sweep_sections,resample_path,
    redirect_local_sections)
from geometry.tube_constraints import capped_translation_transition
from geometry.tube_prism_constraints import StationPrismConstraints
from geometry.tube_biarc import biarc_path
from geometry.tube_continuity import (loft_faces,boundary_loops,join_mouths,continuity_audit,
    finite_wall_cells,annular_rim)
from geometry.tube_tool import _full_mesh_checks
from geometry.audit_alimentary_atlas import audit_catalogue


def straight_sections(stations=81,controls=24,length=.2,radii=(.006,.004)):
    angle=np.arange(controls)*2*np.pi/controls
    xy=np.column_stack([radii[0]*np.cos(angle),radii[1]*np.sin(angle),np.zeros(controls)])
    center=np.column_stack([np.zeros(stations),np.zeros(stations),np.linspace(0,length,stations)])
    return np.array([center[:,None]+xy*f for f in np.linspace(1,.75,9)])


class SectionGeometryTests(unittest.TestCase):
    def test_local_redirection_preserves_fixed_ports_and_section_calibre(self):
        source=straight_sections(stations=101,length=.4)
        centres=source[0].mean(1);u=np.clip((centres[:,2]-.10)/.20,0,1)
        target=centres.copy();target[:,0]+=.008*np.sin(np.pi*u)**4
        target[u==1,0]=centres[u==1,0]
        result,rotation=redirect_local_sections(source,target,dtype=np.float64)
        np.testing.assert_array_equal(result[:,[0,-1]],source[:,[0,-1]])
        np.testing.assert_array_equal(result[:,:20],source[:,:20])
        np.testing.assert_allclose(np.linalg.det(rotation),1,atol=1e-14)
        for key in ('area','shape_distances','thickness'):
            np.testing.assert_allclose(section_metrics(result)[key],section_metrics(source)[key],atol=1e-14)
        self.assertTrue(finite_wall_cells(result.reshape(9,-1,3),loft_faces(101,24))['finite_wall_cells_pass'])
        checked,_,report=sweep_sections(source,target,dtype=np.float64,preserve_source_roll=True)
        self.assertTrue(report['section_transport_pass'],report)
        np.testing.assert_array_equal(checked,result)
        with self.assertRaises(TubeGeometryError):
            sweep_sections(source,target,tangents=np.tile([0,0,1],(101,1)),preserve_source_roll=True)

    def test_local_course_reversal_is_not_hidden_by_profile_flipping(self):
        source=straight_sections()
        with self.assertRaises(TubeGeometryError):
            redirect_local_sections(source,source[0].mean(1)[::-1])

    def test_constrained_translation_preserves_endpoints_and_hard_caps(self):
        desired=np.array([.01,.009,.005,.003,0.]);caps=np.array([.003,.001,.004,.004])
        result=capped_translation_transition(desired,caps)
        np.testing.assert_array_equal(result[[0,-1]],desired[[0,-1]])
        self.assertTrue(np.all(-np.diff(result)>=0))
        self.assertTrue(np.all(-np.diff(result)<=caps+1e-15))
        np.testing.assert_array_equal(capped_translation_transition(desired,np.full(4,.02)),desired)
        # Even a sub-picometre capacity shortfall must terminate as infeasible.
        with self.assertRaises(TubeGeometryError):
            capped_translation_transition([.01,0],[.01-1e-14])
        for values,limits in (([],[]),([1,2],[1]),([1,0],[-1]),([np.nan,0],[1])):
            with self.assertRaises(TubeGeometryError):capped_translation_transition(values,limits)

    def test_helix_frame_converges_to_analytic_bishop_frame(self):
        errors=[]
        for n in (25,49,97):
            phi=np.linspace(0,4,n);R,h=.1,.03;L=np.hypot(R,h)
            p=np.column_stack([R*np.cos(phi),R*np.sin(phi),h*phi])
            t=np.column_stack([-R*np.sin(phi),R*np.cos(phi),np.full(n,h)])/L
            normal=np.column_stack([-np.cos(phi),-np.sin(phi),np.zeros(n)])
            binormal=np.cross(t,normal);theta=h/L*phi
            exact=np.cos(theta)[:,None]*normal-np.sin(theta)[:,None]*binormal
            frames=rotation_minimizing_frames(p,t,normal[0])
            errors.append(np.linalg.norm(frames[:,:,0]-exact,axis=1).max())
            np.testing.assert_allclose(frames.transpose(0,2,1)@frames,np.broadcast_to(np.eye(3),frames.shape),atol=1e-13)
            np.testing.assert_allclose(np.linalg.det(frames),1,atol=1e-13)
        self.assertGreater(errors[0]/errors[1],14)
        self.assertGreater(errors[1]/errors[2],14)

    def test_ellipse_and_all_thicknesses_preserved_for_wide_and_tight_bends(self):
        source=straight_sections()
        for radius in (.025,.060):
            theta=np.linspace(0,np.pi/2,source.shape[1])
            centers=np.column_stack([radius*(1-np.cos(theta)),np.zeros(len(theta)),radius*np.sin(theta)])
            tangents=np.column_stack([np.sin(theta),np.zeros(len(theta)),np.cos(theta)])
            result,_,report=sweep_sections(source,centers,tangents,dtype=np.float64)
            self.assertTrue(report['section_transport_pass'],report)
            self.assertLess(report['maximum_relative_area_error'],1e-12)
            self.assertLess(report['maximum_shape_distance_error_m'],1e-14)
            self.assertFalse(report['radii_scaled'])
            wall=finite_wall_cells(result.reshape(9,-1,3),loft_faces(source.shape[1],source.shape[2]))
            self.assertTrue(wall['finite_wall_cells_pass'],wall)

    def test_impossible_tight_bend_is_reported_without_thinning(self):
        source=straight_sections();theta=np.linspace(0,np.pi/2,81);r=.003
        c=np.column_stack([r*(1-np.cos(theta)),np.zeros(81),r*np.sin(theta)])
        t=np.column_stack([np.sin(theta),np.zeros(81),np.cos(theta)])
        result,_,report=sweep_sections(source,c,t,dtype=np.float64)
        self.assertFalse(report['section_transport_pass'])
        self.assertIn('turn too tight for the actual section',report['failures'])
        np.testing.assert_allclose(section_metrics(result)['area'],section_metrics(source)['area'],rtol=1e-12)
        self.assertFalse(finite_wall_cells(result.reshape(9,-1,3),loft_faces(81,24))['finite_wall_cells_pass'])

    def test_straight_and_inflection_course_have_no_frame_flips(self):
        z=np.linspace(0,.3,201);x=.008*np.sin(z*2*np.pi/.3)
        p=np.column_stack([x,np.zeros(len(z)),z]);t=np.column_stack([.008*2*np.pi/.3*np.cos(z*2*np.pi/.3),np.zeros(len(z)),np.ones(len(z))])
        f=rotation_minimizing_frames(p,t,[0,1,0])
        np.testing.assert_allclose(f[:,:,0],np.broadcast_to([0,1,0],(len(z),3)),atol=1e-14)
        f=rotation_minimizing_frames(np.column_stack([0*z,0*z,z]))
        np.testing.assert_allclose(f,np.broadcast_to(f[0],f.shape),atol=1e-14)

    def test_fillet_radius_is_controlled_and_infeasible_radius_rejected(self):
        wp=np.array([[0,0,0],[.1,0,0],[.1,.1,0]])
        for r in (.02,.04):
            p,t=filleted_path(wp,r,max_step_m=.001)
            np.testing.assert_array_equal(p[[0,-1]],wp[[0,-1]])
            curved=(p[:,0]>.1-r+1e-10)&(p[:,1]<r-1e-10)
            np.testing.assert_allclose(np.linalg.norm(p[curved]-[.1-r,r,0],axis=1),r,atol=1e-14)
            self.assertLess(np.linalg.norm(np.diff(p,axis=0),axis=1).max(),.001000001)
            self.assertLess(np.arccos(np.clip((t[:-1]*t[1:]).sum(1),-1,1)).max(),np.deg2rad(5)+1e-10)
        with self.assertRaises(TubeGeometryError):filleted_path(wp,.11)

    def test_centripetal_curve_handles_nonuniform_steps(self):
        wp=np.array([[0,0,0],[.15,0,0],[.151,.008,0],[.24,.07,.02]])
        p,t=centripetal_curve(wp,chord_tolerance_m=1e-5,max_step_m=.002)
        np.testing.assert_array_equal(p[[0,-1]],wp[[0,-1]])
        self.assertLess(np.linalg.norm(np.diff(p,axis=0),axis=1).max(),.002000001)
        self.assertTrue(np.isfinite(t).all())
        q,qt=resample_path(p,101,t)
        np.testing.assert_array_equal(q[[0,-1]],wp[[0,-1]])
        np.testing.assert_allclose(np.linalg.norm(qt,axis=1),1,atol=1e-14)

    def test_duplicate_points_and_reversed_tangents_are_rejected(self):
        with self.assertRaises(TubeGeometryError):rotation_minimizing_frames([[0,0,0],[0,0,0]])
        with self.assertRaises(TubeGeometryError):rotation_minimizing_frames([[0,0,0],[0,0,1]],[[0,0,1],[0,0,-1]])
        with self.assertRaises(TubeGeometryError):TubeTolerances(curvature_support_limit=1.1)

    def test_real_triangle_checker_on_a_valid_shared_wall(self):
        p=straight_sections(stations=9,controls=12).reshape(9,-1,3)
        validator=Path(__file__).resolve().parents[1]/'scripts/ta98_triangle_intersections.py'
        if not validator.exists():self.skipTest('optional local narrow-phase validator unavailable')
        report=_full_mesh_checks(p,loft_faces(9,12),validator)
        self.assertTrue(report['whole_mesh_intersections_checked'])
        self.assertEqual(len(report['interface_pair_checks']),45)
        self.assertTrue(report['interface_intersection_pass'],report)

    def test_single_outer_surface_can_be_transported_without_claiming_wall(self):
        p=straight_sections()[:1]
        q,_,report=sweep_sections(p,p[0].mean(axis=1),dtype=np.float64)
        self.assertTrue(report['section_transport_pass'],report)
        np.testing.assert_allclose(q,p,atol=1e-14)


class StationConstraintTests(unittest.TestCase):
    def test_continuous_extrema_match_independent_depth_determinants(self):
        fields=straight_sections(stations=7,controls=8).reshape(9,-1,3)
        station=np.repeat(np.arange(7),8);faces=loft_faces(7,8)
        constraint=StationPrismConstraints(fields,faces,station,np.zeros(7,bool))
        rng=np.random.default_rng(203);move=rng.normal(size=(7,3))*.003
        maxima,_=constraint.evaluate(move);moved=fields+move[station]
        for k in range(8):
            a,b=moved[k,faces],moved[k+1,faces];depth=b-a
            brute=np.full(len(faces),-np.inf)
            for t in np.linspace(0,1,2001):
                q=a+t*depth
                determinant=np.einsum('fi,fvi->fv',np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),depth)
                brute=np.maximum(brute,determinant.max(1))
            np.testing.assert_allclose(maxima[k],brute,atol=1e-15,rtol=1e-7)
        independent=finite_wall_cells(moved,faces)
        np.testing.assert_allclose(maxima.max(1),independent['maximum_jacobian_by_interval'],rtol=1e-12,atol=1e-20)

    def test_gradient_matches_station_difference_finite_differences(self):
        fields=straight_sections(stations=6,controls=7).reshape(9,-1,3)
        station=np.repeat(np.arange(6),7);faces=loft_faces(6,7)
        constraint=StationPrismConstraints(fields,faces,station,np.zeros(6,bool))
        rng=np.random.default_rng(104);move=rng.normal(size=(6,3))*.004
        _,gradient=constraint.evaluate(move);face=17;low,high=constraint.low[face],constraint.high[face];eps=1e-7
        for d in range(3):
            perturb=np.zeros_like(move);perturb[high,d]=eps
            plus,_=constraint.evaluate(move+perturb);minus,_=constraint.evaluate(move-perturb)
            np.testing.assert_allclose((plus[:,face]-minus[:,face])/(2*eps),gradient[:,face,d],atol=1e-17,rtol=1e-7)
            perturb[:]=0;perturb[low,d]=eps
            plus,_=constraint.evaluate(move+perturb);minus,_=constraint.evaluate(move-perturb)
            np.testing.assert_allclose((plus[:,face]-minus[:,face])/(2*eps),-gradient[:,face,d],atol=1e-17,rtol=1e-7)

    def test_source_faults_and_malformed_constraints_are_explicit(self):
        fields=straight_sections(stations=4,controls=6).reshape(9,-1,3)
        station=np.repeat(np.arange(4),6);faces=loft_faces(4,6);protected=np.array([True,False,False,True])
        faulty=fields.copy();faulty[[3,4]]=faulty[[4,3]]
        with self.assertRaises(TubeGeometryError):StationPrismConstraints(faulty,faces,station,protected)
        diagnostic=StationPrismConstraints(faulty,faces,station,protected,allow_source_faults=True)
        self.assertGreater(diagnostic.sourceInvalidCounts.sum(),0)
        np.testing.assert_array_equal(diagnostic.limit[~diagnostic.sourceHealthy],diagnostic.sourceMaximum[~diagnostic.sourceHealthy])
        healthy=StationPrismConstraints(fields,faces,station,protected)
        np.testing.assert_allclose(healthy.limit,healthy.sourceMaximum*.1)
        for args in ((fields,faces.astype(float),station,protected),(fields,faces,station,np.zeros(3,bool)),(fields,faces,np.repeat([0,1,2,4],6),np.zeros(5,bool))):
            with self.assertRaises(TubeGeometryError):StationPrismConstraints(*args)
        badfaces=faces.copy();badfaces[0]=[0,1,18]
        with self.assertRaises(TubeGeometryError):StationPrismConstraints(fields,badfaces,station,protected)
        for move in (np.zeros((3,3)),np.full((4,3),np.nan)):
            with self.assertRaises(TubeGeometryError):healthy.evaluate(move)
        result=healthy.project(np.tile([.001,0,0],(4,1)),sweeps=0)
        np.testing.assert_array_equal(result[protected],0)


class BiarcTests(unittest.TestCase):
    def test_known_quarter_circle_and_rigid_3d_equivariance(self):
        radius=.04;args=([0,0,0],[1,0,0],[radius,radius,0],[0,1,0])
        p,t,k,r=biarc_path(*args,minimum_bend_radius_m=.039)
        np.testing.assert_allclose([a['radius_m'] for a in r['arcs']],radius,atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(p-[0,radius,0],axis=1),radius,atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(k,axis=1),1/radius,atol=1e-12)
        self.assertAlmostEqual(r['exact_length_m'],radius*np.pi/2,13)
        rot=np.array([[0,0,1],[1,0,0],[0,1,0]]);offset=np.array([.1,-.2,.3])
        q,qt,qk,_=biarc_path(np.array(args[0])@rot.T+offset,np.array(args[1])@rot.T,np.array(args[2])@rot.T+offset,np.array(args[3])@rot.T,minimum_bend_radius_m=.039)
        np.testing.assert_allclose(q,p@rot.T+offset,atol=1e-14)
        np.testing.assert_allclose(qt,t@rot.T,atol=1e-13)
        np.testing.assert_allclose(qk,k@rot.T,atol=1e-11)

    def test_spatial_join_straight_and_parallel_semicircles(self):
        p,t,k,r=biarc_path([0,0,0],[1,.2,.3],[.04,.03,.02],[.2,1,.1])
        self.assertLess(r['tangent_join_error'],1e-12)
        np.testing.assert_array_equal(p[[0,-1]],[[0,0,0],[.04,.03,.02]])
        np.testing.assert_allclose(np.linalg.norm(t,axis=1),1,atol=1e-13)
        self.assertLess(np.linalg.norm(np.diff(p,axis=0),axis=1).max(),.00200001)
        p,t,k,r=biarc_path([0,0,0],[1,0,0],[.04,0,0],[1,0,0])
        np.testing.assert_array_equal(k,0)
        self.assertAlmostEqual(r['exact_length_m'],.04)
        p,t,k,r=biarc_path([0,0,0],[1,0,0],[0,.04,0],[1,0,0])
        np.testing.assert_allclose([a['radius_m'] for a in r['arcs']],.01,atol=1e-14)
        self.assertLess(r['tangent_join_error'],1e-12)
        np.testing.assert_allclose(t[[0,-1]],[[1,0,0],[1,0,0]],atol=1e-14)

    def test_impossible_radius_and_backward_ports_do_not_thin_or_flip(self):
        with self.assertRaises(TubeGeometryError):
            biarc_path([0,0,0],[1,0,0],[.04,.04,0],[0,1,0],minimum_bend_radius_m=.041)
        with self.assertRaises(TubeGeometryError):
            biarc_path([0,0,0],[1,0,0],[-.04,0,0],[1,0,0])
        with self.assertRaises(TubeGeometryError):
            biarc_path([0,0,0],[1,0,0],[.04,0,0],[1,0,0],first_tangent_distance_m=.05)
        with self.assertRaises(TubeGeometryError):biarc_path([0,0,0],[0,0,0],[.04,0,0],[1,0,0])


class ContinuityTests(unittest.TestCase):
    def setUp(self):
        self.a=straight_sections(stations=3,length=.02).reshape(9,-1,3)
        self.f=loft_faces(3,24);self.am=np.arange(48,72);self.bm=np.arange(24)
        self.b=self.a.copy()+[0,0,.02]

    def test_join_has_shared_coordinates_and_only_two_open_ends(self):
        # Phase change and reversed input mouth order do not create another tube.
        bank,f,mapping,report=join_mouths(self.a,self.f,self.am,self.b,self.f,np.roll(self.bm[::-1],7),position_tolerance_m=0)
        self.assertEqual(report['shared_controls'],24)
        self.assertFalse(report['central_lumen_capped'])
        np.testing.assert_array_equal(bank[:,:len(self.a[0])],self.a)
        np.testing.assert_array_equal(mapping[self.bm],self.am)
        self.assertEqual(continuity_audit(bank[-1],f)['mouths'],2)
        self.assertTrue(finite_wall_cells(bank,f)['finite_wall_cells_pass'])

    def test_tolerance_applies_only_to_explicit_collar(self):
        self.b[:,self.bm,0]+=2e-7
        bank,f,mapping,report=join_mouths(self.a,self.f,self.am,self.b,self.f,self.bm,position_tolerance_m=3e-7)
        self.assertAlmostEqual(report['maximum_gap_m'],2e-7,places=12)
        np.testing.assert_array_equal(bank[:,mapping[24:]],self.b[:,24:])
        with self.assertRaises(TubeGeometryError):join_mouths(self.a,self.f,self.am,self.b,self.f,self.bm,position_tolerance_m=1e-7)

    def test_calibre_mismatch_is_not_welded_or_scaled(self):
        self.b[:,:,:2]*=1.15
        with self.assertRaises(TubeGeometryError):join_mouths(self.a,self.f,self.am,self.b,self.f,self.bm,position_tolerance_m=.01,area_tolerance_relative=.01)

    def test_matching_collar_does_not_hide_a_course_kink(self):
        self.b[:,24:,0]+=.01
        with self.assertRaises(TubeGeometryError):join_mouths(self.a,self.f,self.am,self.b,self.f,self.bm)

    def test_damaged_rim_is_rejected_against_expected_collars(self):
        report=continuity_audit(self.a[0],self.f[:-1],expected_collars=[self.a[0,:24],self.a[0,48:]])
        self.assertFalse(report['topological_continuity_pass'])
        self.assertTrue(report['expected_collars_checked'])

    def test_disconnected_duplicate_bodies_and_blind_plug_are_rejected(self):
        bank=np.r_[self.a[0],self.a[0]+[.1,0,0]]
        self.assertFalse(continuity_audit(bank,np.r_[self.f,self.f+len(self.a[0])])['topological_continuity_pass'])
        v=np.r_[self.a[0],self.a[0,self.am].mean(0)[None]]
        cap=np.array([[len(v)-1,int(b),int(a)] for a,b in zip(self.am,np.roll(self.am,-1))])
        result=continuity_audit(v,np.r_[self.f,cap])
        self.assertEqual(result['mouths'],1)
        self.assertFalse(result['topological_continuity_pass'])

    def test_caps_are_annular_only_and_closed_wall_is_a_torus(self):
        n=self.a.shape[1];outer=self.f;inner=self.f[:,[0,2,1]]+n
        rim=np.concatenate([annular_rim(r,r+n) for r in boundary_loops(outer)])
        shell=np.r_[outer,inner,rim];vertices=np.r_[self.a[0],self.a[1]]
        self.assertTrue(continuity_audit(vertices,shell,expected_mouths=0,expected_euler=0)['topological_continuity_pass'])
        for tri in rim:
            phases=tri%24
            self.assertLessEqual(len(set(phases)),2) # No triangle traverses the passage.

    def test_cli_failure_retains_report_and_never_installs(self):
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as d:
            path=Path(d);bank=straight_sections(stations=5).reshape(9,-1,3);faces=loft_faces(5,24)
            np.savez(path/'broken.npz',fields=bank,faces=np.delete(faces,50,axis=0))
            result=subprocess.run([sys.executable,'-B','-m','geometry.tube_tool','audit','--source',str(path/'broken.npz'),'--report',str(path/'report.json')],cwd=root,capture_output=True,text=True)
            self.assertEqual(result.returncode,2,result.stdout+result.stderr)
            report=json.loads((path/'report.json').read_text());self.assertFalse(report['construction_checks_pass']);self.assertFalse(report['atlas_modified'])

    def test_cli_does_not_overwrite_input(self):
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'bank.npz';np.savez(path,fields=self.a,faces=self.f);original=path.read_bytes()
            result=subprocess.run([sys.executable,'-B','-m','geometry.tube_tool','audit','--source',str(path),'--report',str(path)],cwd=root,capture_output=True,text=True)
            self.assertEqual(result.returncode,2)
            self.assertEqual(path.read_bytes(),original)

    def test_atlas_audit_extracts_only_lumen_and_ignores_suppressed_duplicates(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder=root/'.local-models';folder.mkdir()
            n=self.a.shape[1];positions=np.r_[self.a[0],self.a[-1]].astype('<f4')
            signs=np.r_[np.ones(n,'i1'),-np.ones(n,'i1')]
            rims=np.concatenate([annular_rim(r,r+n) for r in boundary_loops(self.f)])
            faces=np.r_[self.f,self.f[:,[0,2,1]]+n,rims].astype('<u4')
            blob=positions.tobytes()+signs.tobytes()+faces.tobytes();(folder/'wall.bin').write_bytes(blob)
            part=dict(id='OESOPHWALL:male:native:epithelium',chunk=0,positions=0,normalSigns=positions.nbytes,
                indices=positions.nbytes+signs.nbytes,vertexCount=len(positions),indexCount=faces.size)
            duplicate=dict(part,id='OESOPHWALL:male:connector:epithelium',suppressed=True)
            cat=dict(chunks=[dict(url='/local-models/wall.bin')],parts=[part,duplicate]);p=root/'atlas.json';p.write_text(json.dumps(cat))
            report=audit_catalogue(p,root)
            self.assertEqual(len(report['lumen_sources']),1)
            self.assertEqual(report['lumen_sources'][0]['lumen_faces'],len(self.f))
            self.assertTrue(report['topological_continuity_pass'])
            self.assertFalse(report['oesophagus_to_anus_pass'])
            duplicate['suppressed']=False;p.write_text(json.dumps(cat));report=audit_catalogue(p,root)
            self.assertGreater(report['duplicate_faces'],0)
            self.assertFalse(report['topological_continuity_pass'])

    def test_zero_thickness_wall_is_not_accepted(self):
        a=np.repeat(self.a[:1],2,axis=0)
        self.assertFalse(finite_wall_cells(a,self.f)['finite_wall_cells_pass'])
        p=a.reshape(2,3,24,3)
        _,_,report=sweep_sections(p,p[0].mean(axis=1),dtype=np.float64)
        self.assertFalse(report['section_transport_pass'])
        self.assertIn('source wall has zero or insufficient correspondence thickness',report['failures'])


if __name__=='__main__':unittest.main()
