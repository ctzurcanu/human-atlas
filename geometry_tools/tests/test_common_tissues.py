import unittest,tempfile
from pathlib import Path
import numpy as np
import manifold3d as md
from geometry_tools.boolean_cells import BooleanCells
from geometry_tools.harmonic_wall import layers,tetrahedralize
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_tissue_axis_split import split_shared_cells_at_axis


def arrays(m):
    q=m.to_mesh64();return np.asarray(q.vert_properties)[:,:3],np.asarray(q.tri_verts,int)


class SharedTissueTests(unittest.TestCase):
    def test_regional_annulus_is_preserved_without_a_lumen_disk(self):
        n=8;angle=np.arange(n)*2*np.pi/n
        v=np.array([[r*np.cos(a),r*np.sin(a),z]for z in [0.,1.,2.]for r in [1.,.8]for a in angle]);faces=[]
        for j in range(2):
            for k in range(n):
                nxt=(k+1)%n
                for offset,reverse in [(0,False),(n,True)]:
                    a=j*2*n+offset+k;b=j*2*n+offset+nxt;c=b+2*n;d=a+2*n
                    f=np.array([[a,b,c],[a,c,d]])
                    faces.extend(f[:,[0,2,1]]if reverse else f)
        rim=[]
        for j in [0,1,2]:
            fs=[]
            for k in range(n):
                a=j*2*n+k;b=j*2*n+(k+1)%n;fs.extend([[a,b,b+n],[a,b+n,a+n]])
            fs=np.array(fs)
            if j==0:faces.extend(fs[:,[0,2,1]])
            elif j==2:faces.extend(fs)
            else:rim=fs
        faces=np.array(faces);p,t,_,_,e=tetrahedralize(v,faces,rim)
        self.assertTrue(e['boundary_faces_literal'])
        radial=np.linalg.norm(v[:,:2],axis=1);outer=faces[np.all(abs(radial[faces]-1)<1e-12,axis=1)];inner=faces[np.all(abs(radial[faces]-.8)<1e-12,axis=1)]
        r=layers(v,faces,outer,inner,['a','b','c'],tetrahedra=(p,t),method='boundary-harmonic')
        for uses in r['cells'].values():
            fs=np.concatenate([r['surfaces'][u['surface']]if u['side']==1 else r['surfaces'][u['surface']][:,[0,2,1]]for u in uses])
            self.assertTrue(audit(Mesh(r['positions'],fs,'fixture'),require_closed=True,expected_genus=[1])['passed'])

    def test_source_face_overlay_and_durable_replay(self):
        b=BooleanCells()
        def s(m,n):return b.source(*arrays(m),n)
        k=s(md.Manifold.cube([2,2,2]),'outer');i=s(md.Manifold.cube([1.6,1.6,1.6]).translate([.2,.2,.2]),'inner');p=s(md.Manifold.cube([1,3,3]).translate([.8,-.5,-.5]),'partition')
        solids={'shell':k-i,'left':i-p,'right':i^p}
        with tempfile.TemporaryDirectory()as td:
            path=Path(td)/'cells.npz';b.save(solids,path);replay,inputs=BooleanCells.load(path);r=replay.assemble(inputs)
        volumes={}
        for name,uses in r['cells'].items():
            f=np.concatenate([r['surfaces'][q['surface']]if q['side']==1 else r['surfaces'][q['surface']][:,[0,2,1]]for q in uses]);q=audit(Mesh(r['positions'],f,'fixture'),require_closed=True)
            self.assertTrue(q['passed'],q['failures']);volumes[name]=q['topology']['oriented_triangle_integral_m3']
        self.assertAlmostEqual(sum(volumes.values()),8.)
        self.assertAlmostEqual(volumes['shell'],8.-1.6**3)

    def test_harmonic_annular_coats_preserve_the_lumen(self):
        wall=md.Manifold.cylinder(2.,1.,1.,16)-md.Manifold.cylinder(2.,.8,.8,16);v,f=arrays(wall)
        radial=np.linalg.norm(v[:,:2],axis=1);outer=f[np.all(abs(radial[f]-1)<1e-12,axis=1)];inner=f[np.all(abs(radial[f]-.8)<1e-12,axis=1)]
        r=layers(v,f,outer,inner,['outer','middle','inner']);total=0.
        for name,uses in r['cells'].items():
            fs=np.concatenate([r['surfaces'][q['surface']]if q['side']==1 else r['surfaces'][q['surface']][:,[0,2,1]]for q in uses]);q=audit(Mesh(r['positions'],fs,'fixture'),require_closed=True,expected_genus=[1])
            self.assertTrue(q['passed'],q['failures']);total+=q['topology']['oriented_triangle_integral_m3']
        self.assertAlmostEqual(total,wall.volume(),places=12)
        self.assertTrue(r['evidence']['boundary_faces_literal'])
        # The adjoining cell stays whole while receiving the same subdivided
        # interface when another cell is split by an oblique regional plane.
        q=split_shared_cells_at_axis(r['positions'],r['surfaces'],r['cells'],{'middle':['a','b']},normal=[.1,.2,1.],coordinate=1.)
        for hit in q['edge_intersections']:
            a,b=hit['edge'];t=hit['fraction']
            np.testing.assert_allclose(q['positions'][hit['control']],r['positions'][a]*(1-t)+r['positions'][b]*t,atol=1e-14,rtol=0)
        for sid,parent in q['source_surface_ancestry'].items():
            self.assertEqual(len(parent),len(q['surfaces'][sid]))
            self.assertTrue(np.all((parent>=0)&(parent<len(r['surfaces'][sid.rsplit(':',1)[0]]))))
        for name,uses in q['cells'].items():
            fs=np.concatenate([q['surfaces'][u['surface']]if u['side']==1 else q['surfaces'][u['surface']][:,[0,2,1]]for u in uses])
            self.assertTrue(audit(Mesh(q['positions'],fs,'fixture'),require_closed=True)['passed'])

    def test_unequal_coat_fractions_preserve_volume_and_reject_bad_order(self):
        wall=md.Manifold.cylinder(.02,.01,.01,8)-md.Manifold.cylinder(.02,.008,.008,8);v,f=arrays(wall)
        radial=np.linalg.norm(v[:,:2],axis=1);outer=f[np.all(abs(radial[f]-.01)<1e-12,axis=1)];inner=f[np.all(abs(radial[f]-.008)<1e-12,axis=1)]
        names=['covering','muscle','mucosa'];fractions=[.08,.92]
        r=layers(v,f,outer,inner,names,method='boundary-harmonic',interface_fractions=fractions)
        total=0.
        for uses in r['cells'].values():
            fs=np.concatenate([r['surfaces'][u['surface']]if u['side']==1 else r['surfaces'][u['surface']][:,[0,2,1]]for u in uses])
            check=audit(Mesh(r['positions'],fs,'fixture'),require_closed=True,expected_components=1,expected_genus=[1])
            self.assertTrue(check['passed']);total+=check['topology']['oriented_triangle_integral_m3']
        self.assertAlmostEqual(total,wall.volume(),places=15)
        self.assertEqual(r['evidence']['interface_fractions'],fractions)
        for bad in [[.8,.2],[0.,.8],[.2,1.],[.5],[-.1,.6],[.2,float('nan')]]:
            with self.assertRaises(ValueError):layers(v,f,outer,inner,names,interface_fractions=bad)

    def test_inner_coat_covers_the_entire_lining_when_remainder_contains_end_annuli(self):
        wall=md.Manifold.cylinder(.02,.01,.01,8)-md.Manifold.cylinder(.02,.008,.008,8);v,f=arrays(wall)
        radius=np.linalg.norm(v[:,:2],axis=1)
        is_inner=np.all(abs(radius[f]-.008)<1e-12,axis=1)
        inner=f[is_inner];remainder=f[~is_inner]
        lips=np.intersect1d(np.unique(inner),np.unique(remainder))
        r=layers(v,f,remainder,inner,['outer','middle','mucosa'],method='boundary-harmonic',rim_inner_vertices=lips)
        inner_ids=set(map(int,np.flatnonzero(is_inner)));covered=set()
        for name,uses in r['cells'].items():
            for u in uses:
                ids=set(map(int,r['source_faces'][u['surface']]))&inner_ids
                if ids:self.assertEqual(name,'mucosa');covered|=ids
        self.assertEqual(covered,inner_ids)
        for name,uses in r['cells'].items():
            fs=np.concatenate([r['surfaces'][u['surface']]if u['side']==1 else r['surfaces'][u['surface']][:,[0,2,1]]for u in uses])
            self.assertTrue(audit(Mesh(r['positions'],fs,'fixture'),require_closed=True,expected_genus=[1])['passed'])
        with self.assertRaises(ValueError):layers(v,f,remainder,inner,['a','b'],method='boundary-harmonic',rim_inner_vertices=lips,free_source_vertices=lips)
        with self.assertRaises(ValueError):layers(v,f,remainder,inner,['a','b'],method='distance-ratio',rim_inner_vertices=lips)


if __name__=='__main__':unittest.main()
