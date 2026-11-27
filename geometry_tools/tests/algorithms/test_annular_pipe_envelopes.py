import unittest
import numpy as np
import manifold3d as m3d
from geometry_tools.algorithms.annular_pipe_envelopes import recover,recover_paired_sheets,fluid_through_mouths,PipeEnvelopeError

def fixture():
    n,k=12,9;theta=np.arange(n)*2*np.pi/n;x=np.linspace(0,.014,k);v=[]
    for radius in [.002,.0014]:
        for p in x:v.extend(np.column_stack([np.full(n,p),1.2+radius*np.cos(theta)+.0003*np.sin(p*130),radius*np.sin(theta)]))
    v=np.asarray(v).astype('<f4').astype(float);f=[]
    for layer in range(2):
        for i in range(k-1):
            for j in range(n):
                a=layer*n*k+i*n+j;b=layer*n*k+i*n+(j+1)%n;c=b+n;d=a+n
                f.extend([[a,b,c],[a,c,d]]if layer==0 else[[a,c,b],[a,d,c]])
    for i in [0,k-1]:
        for j in range(n):
            a=i*n+j;b=i*n+(j+1)%n;ia=a+n*k;ib=b+n*k
            f.extend([[a,ia,ib],[a,ib,b]]if i==0 else[[a,ib,ia],[a,b,ib]])
    return v,np.asarray(f)

def solid(v,f):
    return m3d.Manifold(m3d.Mesh64(np.ascontiguousarray((v-[0,1.2,0])*1000),np.asarray(f,dtype='<u8')))

class EnvelopeTests(unittest.TestCase):
    def test_bent_finite_wall_is_recovered_by_the_exact_shared_disk_difference(self):
        v,f=fixture();r=recover(v,f);source=solid(r.positions,r.original_faces);outer=solid(r.positions,r.outer_faces);fluid=solid(r.positions,r.fluid_faces)
        self.assertEqual(source.status(),m3d.Error.NoError);self.assertEqual(outer.status(),m3d.Error.NoError);self.assertEqual(fluid.status(),m3d.Error.NoError)
        difference=(source-(outer-fluid))+((outer-fluid)-source)
        self.assertLess(abs(difference.volume()),1e-10)
        self.assertLess(abs((source^fluid).volume()),1e-10)
        self.assertEqual(r.report['annularCapFaces'],48);self.assertEqual(r.report['mouths'],2)
    def test_solid_cylinder_is_not_mislabeled_as_a_hollow_wall(self):
        v,f=fixture();r=recover(v,f)
        with self.assertRaises(PipeEnvelopeError):recover(r.positions,r.outer_faces)
    def test_explicit_sheet_lineage_preserves_a_nonplanar_native_annular_end(self):
        v,f=fixture();theta=np.arange(12)*2*np.pi/12;v[-12:,0]+=.00002*np.cos(theta*2)
        with self.assertRaises(PipeEnvelopeError):recover(v,f)
        r=recover_paired_sheets(v,f,(np.arange(len(v))>=len(v)//2).astype(int))
        source=solid(r.positions,r.original_faces);difference=solid(r.positions,r.outer_faces)-solid(r.positions,r.fluid_faces)
        self.assertLess(abs(((source-difference)+(difference-source)).volume()),1e-10)
    def test_unclosed_surface_is_rejected_instead_of_filled_arbitrarily(self):
        v,f=fixture()
        with self.assertRaises(PipeEnvelopeError):recover(v,f[:-1])
    def test_through_mouth_cutters_leave_a_planar_native_wall_unchanged(self):
        v,f=fixture();r=recover(v,f);p,q,rows=fluid_through_mouths(r)
        source=solid(r.positions,r.original_faces);wall=solid(r.positions,r.outer_faces)-solid(p,q)
        self.assertLess(abs(((source-wall)+(wall-source)).volume()),1e-10)
        self.assertTrue(all(row['extensionMeters']>0 for row in rows))

if __name__=='__main__':unittest.main()
