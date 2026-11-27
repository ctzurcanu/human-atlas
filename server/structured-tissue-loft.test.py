"""Physical volume and orientation gates for a shared section construction."""
import unittest
import numpy as np
from geometry_tools.algorithms.radial_section_loft import prism_tetrahedra,harmonic_sections,bound_section_deformations,transport_section_sequence,radial_carrier_section,carrier_prism_edge_patches
from geometry_tools.algorithms.labelled_tetrahedra import partition_labelled_tetrahedra
from geometry_tools.implicit_cells import Arrangement
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

class StructuredLoft(unittest.TestCase):
    def carrier(self):
        n=16;angle=np.arange(n)*2*np.pi/n;circle=np.c_[np.cos(angle),np.sin(angle)]
        fractions=np.array([.04,.1,.25,.45,.65,.82,1.])
        p=np.vstack([[[0.,0.]],*[t*circle for t in fractions],circle*1.05]);f=[]
        for j in range(n):f.append([0,1+j,1+(j+1)%n])
        for k in range(len(fractions)):
            a=1+k*n;b=a+n
            for j in range(n):
                jj=(j+1)%n;f.extend([[a+j,b+j,b+jj],[a+j,b+jj,a+jj]])
        f=np.array(f);owners=np.zeros(len(p),int)
        for i,tri in enumerate(f):owners[tri]=i
        q=p[f[owners]];mat=np.stack([q[:,1]-q[:,0],q[:,2]-q[:,0]],axis=2)
        uv=np.linalg.solve(mat,(p-q[:,0])[...,None])[...,0];bary=np.c_[1-uv.sum(1),uv]
        template=dict(positions=p,faces=f,carrier_positions=p,carrier_faces=f,carrier_owner=owners,
            carrier_barycentric=bary,triangle_carrier=np.arange(len(f)),radial_fractions=fractions)
        grey=(np.linalg.norm(p[f],axis=2).max(1)<=.45000001)
        return template,circle,angle,grey

    def test_regional_carrier_preserves_literal_pial_contours_and_shared_tissue(self):
        template,circle,angle,grey=self.carrier();center=np.array([.1,-.2])
        radius=.004*(1+.12*np.cos(3*angle));inner=center+circle*radius[:,None]
        outer=center+circle*(radius+.00005)[:,None]
        p,proof=radial_carrier_section(template,inner,outer,center,grey,.2)
        np.testing.assert_allclose(p[-32:-16],inner,rtol=0,atol=1e-16)
        np.testing.assert_allclose(p[-16:],outer,rtol=0,atol=1e-16)
        self.assertTrue(proof['all_triangles_positive'])
        self.assertAlmostEqual(proof['actual_grey_fraction'],.2,places=7)
        self.assertTrue(proof['common_piecewise_affine_map'])

    def test_regional_carrier_refuses_a_pial_wall_with_reversed_thickness(self):
        template,circle,angle,grey=self.carrier()
        with self.assertRaisesRegex(ValueError,'inverted element'):
            radial_carrier_section(template,circle,circle*.99,[0,0],grey,.2)

    def test_embedded_edge_follows_all_three_volume_maps_without_interface_cracks(self):
        lower=np.array([[0.,0.,0.],[1.,0.,0.],[0.,0.,1.]])
        upper=np.array([[.1,1.,.2],[1.1,1.,0.],[.2,1.,1.2]])
        lam=np.array([[.2,.3,.5],[.6,.1,.3]])
        patches=carrier_prism_edge_patches(lam,lower,upper)
        self.assertEqual(len(patches),3)
        for a,b in zip(patches,patches[1:]):
            np.testing.assert_allclose(a[1][2],b[0][2],atol=1e-15)
            np.testing.assert_allclose(a[2][2],b[3][2],atol=1e-15)
        np.testing.assert_allclose(patches[0][0][2],lam[0]@lower,atol=1e-15)
        np.testing.assert_allclose(patches[-1][1][2],lam[0]@upper,atol=1e-15)
        # Every polygon belongs to a single affine tetrahedron and remains
        # planar. A direct low/high chord generally has a different path.
        for poly in patches:
            p=np.array([row[2]for row in poly]);normal=np.cross(p[1]-p[0],p[2]-p[0])
            self.assertAlmostEqual(float(np.dot(normal,p[3]-p[0])),0.,places=14)

    def test_identity_transport_preserves_nonuniform_reference_controls(self):
        p=np.array([[0.,0.],[2.,0.],[2.,1.],[0.,1.],[.31,.27]])
        f=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
        fixed=np.array([0,1,2,3])
        result,proof=harmonic_sections(p,f,fixed,p[None,fixed],transport_displacements=True)
        np.testing.assert_array_equal(result[0],p)
        self.assertTrue(proof['all_triangles_positive'])

    def test_mean_value_transport_reproduces_affine_coordinates(self):
        p=np.array([[0.,0.],[2.,0.],[2.,1.],[0.,1.],[.31,.27]])
        f=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
        affine=np.array([[1.3,.1],[.2,.8]]);target=p@affine+[.1,-.2]
        fixed=np.arange(4)
        result,proof=harmonic_sections(p,f,fixed,target[None,fixed],weights='mean-value')
        np.testing.assert_allclose(result[0],target,atol=1e-13,rtol=0)
        self.assertEqual(proof['reference_weights'],'mean-value')

    def test_section_fit_retains_its_explicit_distance_and_orientation_limits(self):
        p=np.array([[[0.,0.],[1.,0.],[0.,1.]]]);target=p.copy();target[0,2]=[0.,-1.]
        result,proof=bound_section_deformations(p,target,[[0,1,2]],maximum_displacement_m=.1)
        self.assertFalse(proof['passed'])
        self.assertGreaterEqual(proof['minimum_relative_area'],.2-1e-12)
        self.assertGreater(proof['maximum_target_deviation_m'],.1)
        same,proof=bound_section_deformations(p,p,[[0,1,2]],maximum_displacement_m=.1)
        np.testing.assert_array_equal(same,p)
        self.assertTrue(proof['passed'])

    def test_successive_transport_preserves_shared_affine_stations(self):
        p=np.array([[0.,0.],[2.,0.],[2.,1.],[0.,1.],[.31,.27]])
        f=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]]);fixed=np.arange(4)
        expected=np.stack([p*(1+i*.02)+[i*.03,-i*.01]for i in range(7)])
        result,proof=transport_section_sequence(p,f,fixed,expected[:,fixed],3)
        np.testing.assert_allclose(result,expected,atol=1e-12,rtol=0)
        self.assertTrue(proof['all_triangles_positive'])

    def test_prism_replays_area_times_length(self):
        p=np.array([[0.,0.],[1.,0.],[0.,1.]])
        v,t,parent,proof=prism_tetrahedra(np.array([p,p+[.2,.1]]),[0.,2.],[[0,1,2]])
        self.assertTrue(proof['all_positive_reference_jacobians'])
        self.assertAlmostEqual(proof['volume_m3'],1.)
        self.assertTrue(np.all(np.linalg.det(v[t[:,1:]]-v[t[:,:1]])>0))

    def test_positive_end_sections_do_not_authorize_a_folded_loft(self):
        p=np.array([[0.,0.],[1.,0.],[0.,1.]])
        with self.assertRaisesRegex(ValueError,'negative reference Jacobian'):
            prism_tetrahedra(np.array([p,-p]),[0.,1.],[[0,1,2]])

    def test_domain_pruning_requires_an_affine_dominance_proof(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        scores=np.array([[1.,-1.],[0.,0.],[0.,0.],[0.,0.]])
        result=partition_labelled_tetrahedra(p,[[0,1,2,3]],scores,allowed_tetrahedron_labels=[(0,)])
        self.assertEqual(len(result['cells'][0]),4);self.assertEqual(len(result['cells'][1]),0)
        self.assertTrue(result['evidence']['explicit_domain_restriction_proofs'][0]['omitted_fields_affine_dominated'])
        scores[0,1]=2.
        with self.assertRaisesRegex(ValueError,'not strictly volume-dominated'):
            partition_labelled_tetrahedra(p,[[0,1,2,3]],scores,allowed_tetrahedron_labels=[(0,)])

    def test_material_boundary_keeps_the_neighboring_field_trace(self):
        p=np.array([[0.,0.,0.],[0.,1.,0.],[0.,0.,1.],[-1.,0.,0.],[1.,0.,0.]])
        t=np.array([[0,2,1,3],[0,1,2,4]])
        values=np.c_[np.full(len(p),-1.),p[:,1]-.3]
        def classify(value):
            if arrangement.current_tet_index==0:return'grey'
            return'white-upper'if value[1]>0 else'white-lower'
        arrangement=Arrangement(p,values,classify)
        positions,surfaces,cells=arrangement.build(t);volume=0.
        self.assertEqual(set(cells),{'grey','white-upper','white-lower'})
        for uses in cells.values():
            f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
            proof=topology(Mesh(positions,f,'test:metres'));self.assertTrue(proof['closed_oriented_manifold'])
            volume+=proof['oriented_triangle_integral_m3']
        self.assertAlmostEqual(volume,1/3)
        grey_neighbors={u.get('neighbor')for u in cells['grey']if u.get('neighbor')}
        self.assertEqual(grey_neighbors,{'white-upper','white-lower'})

if __name__=='__main__':unittest.main()
