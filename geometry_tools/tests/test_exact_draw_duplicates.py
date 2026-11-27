"""Guard the appearance and winding boundaries of literal draw cleanup."""
import json
from pathlib import Path
import unittest
import numpy as np

from geometry_tools.algorithms.exact_draw_duplicates import exact_draw_duplicates
from geometry_tools.model import Mesh,AtlasReader


def fixture(*,reverse=False,field=None,near=False,copies=2):
    tri=np.array([[0.,0.,0.],[.001,0.,0.],[0.,.001,0.]])
    vertices=np.tile(tri,(copies,1));faces=np.arange(3*copies).reshape(-1,3)
    faces[1]=faces[1,[1,2,0]] if not reverse else faces[1,[0,2,1]]
    attrs={'normals':np.tile([0,0,32767],(3*copies,1)).astype('<i2'),
           'normalSigns':np.ones(3*copies,'i1'),
           'uvs':np.tile([[0.,0.],[1.,0.],[0.,1.]],(copies,1)).astype('<f4'),
           'colors':np.full((3*copies,3),128,'u1')}
    if field is not None:attrs[field][3,0]+=1
    if near:vertices[3,2]=1e-12
    return Mesh(vertices,faces,'test:metres',attrs)


class ExactDrawDuplicateRules(unittest.TestCase):
    def test_cyclic_winding_and_all_source_fields_are_identical(self):
        m=fixture();keep,p=exact_draw_duplicates(m)
        np.testing.assert_array_equal(keep,[0]);self.assertEqual(len(p['duplicateCopiesRemoved']),1)
        r=p['duplicateCopiesRemoved'][0]
        np.testing.assert_array_equal(m.vertices[m.faces[r['removedFace'],r['sourceCornerOrder']]],m.vertices[m.faces[0]])

    def test_opposed_sheet_is_preserved(self):
        self.assertEqual(len(exact_draw_duplicates(fixture(reverse=True))[0]),2)

    def test_different_normals_are_preserved(self):
        m=fixture(field='normals');self.assertEqual(len(exact_draw_duplicates(m)[0]),2)

    def test_different_uv_is_preserved(self):
        self.assertEqual(len(exact_draw_duplicates(fixture(field='uvs'))[0]),2)

    def test_different_color_is_preserved(self):
        self.assertEqual(len(exact_draw_duplicates(fixture(field='colors'))[0]),2)

    def test_nearby_geometry_is_never_welded(self):
        self.assertEqual(len(exact_draw_duplicates(fixture(near=True))[0]),2)

    def test_every_removed_face_maps_to_a_surviving_source_face(self):
        keep,p=exact_draw_duplicates(fixture(copies=3));self.assertEqual(len(keep),1)
        self.assertEqual(len(p['duplicateCopiesRemoved']),2)
        for r in p['duplicateCopiesRemoved']:self.assertIn(r['retainedSourceFace'],keep)


BASE=Path('.local-models/geometry-workbench/biodigital-male-circulation-v1/blood-installed-v1915')


@unittest.skipUnless((BASE/'circulation-atlas.json').exists(),'Local unredistributed packed fixture')
class PackedCiliaryCleanup(unittest.TestCase):
    def test_actual_bank_fields_and_surface_support_remain_literal(self):
        p=json.loads((BASE/'prepare-report.json').read_text());old=json.loads(Path(p['sourceCatalogue']).read_text());new=json.loads((BASE/'circulation-atlas.json').read_text());source=json.loads(Path(p['sourceReport']).read_text());ids={r['partId']for r in source['rows']};reader=AtlasReader()
        self.assertEqual(old['concepts'],new['concepts']);self.assertEqual(old['chunks'],new['chunks'][:-1])
        allowed={'positions','normals','normalSigns','uvs','colors','indices','chunk','bounds','indexCount'}
        for a,b in zip(old['parts'],new['parts']):
            if a['id']not in ids:self.assertEqual(a,b);continue
            self.assertEqual({k:v for k,v in a.items()if k not in allowed},{k:v for k,v in b.items()if k not in allowed})
            self.assertEqual(a['vertexCount'],b['vertexCount']);self.assertEqual(a['indexCount']-b['indexCount'],180)
            banks=[reader._bank(reader.asset(x['chunks'][q['chunk']]['url']))for x,q in [(old,a),(new,b)]]
            for key,dtype,width in [('positions','<f4',3),('normals','<i2',3),('normalSigns','i1',1),('uvs','<f4',2),('colors','u1',3)]:
                if key in a:np.testing.assert_array_equal(reader.view(banks[0],a[key],a['vertexCount'],dtype,width),reader.view(banks[1],b[key],b['vertexCount'],dtype,width))
            meshes=[reader.read_part(x,q,frame='model:local-ta98:assembled-metres',catalogue_path=path)
                    for x,q,path in [(old,a,p['sourceCatalogue']),(new,b,BASE/'circulation-atlas.json')]]
            def oriented_triangle_key(triangle):
                return min(tuple(np.roll(triangle,i,axis=0).ravel())for i in range(3))
            self.assertEqual({oriented_triangle_key(t)for t in meshes[0].triangles},
                             {oriented_triangle_key(t)for t in meshes[1].triangles})


if __name__=='__main__':unittest.main()
