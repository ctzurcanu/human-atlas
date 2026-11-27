"""A shared-only bank retains the exact emitted anatomy and normal signs."""
from pathlib import Path
import tempfile,unittest
import numpy as np
from geometry_tools.tissue import layer_arrangement,normal_words,pack,decode_packed
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.finite_bank_artifacts import pack_shared_only

class SharedArtifacts(unittest.TestCase):
    def test_small_cell_decodes_only_its_used_vertices_without_changing_attributes(self):
        v=np.array([[0.,0.,0.],[.01,0,0],[0,.01,0],[0,0,.01],
                    [.1,0,0],[.11,0,0],[.1,.01,0],[.1,0,.01]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        surfaces={'a':f,'b':f+4};cells={k:[dict(surface=k,side=1)]for k in surfaces}
        packed=pack_finite_cells(v,surfaces,normal_words(v,surfaces),cells,dict(a='a',b='b'))
        decoded=decode_packed(packed,'test:metres',[0,0,0])
        for name,mesh in decoded.items():
            row=packed.cells[name];indices=np.ndarray((row['indexCount']//3,3),'<u4',buffer=packed.data,offset=row['indices'])
            raw_v=np.ndarray((row['vertexCount'],3),'<f4',buffer=packed.data,offset=row['positions'])
            raw_n=np.ndarray((row['vertexCount'],3),'<i2',buffer=packed.data,offset=row['normals'])
            self.assertEqual(len(mesh.vertices),len(np.unique(indices)))
            self.assertLess(len(mesh.vertices),row['vertexCount'])
            np.testing.assert_array_equal(mesh.vertices[mesh.faces],raw_v[indices])
            np.testing.assert_array_equal(mesh.attributes['normals'][mesh.faces],raw_n[indices])
        # A chart's normal words are independent of unrelated unused anatomy.
        isolated=normal_words(v[:4],{'a':f})['a']
        np.testing.assert_array_equal(normal_words(v,surfaces)['a'],isolated)

    def test_decodable_bank_is_identical_without_per_cell_array_copies(self):
        square=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])*.01
        fields=np.array([square+[0,0,h]for h in [.003,.002,.001,0.]])
        v,s,c,_=layer_arrangement(fields,np.array([[0,1,2],[0,2,3]]),['outer','middle','inner'])
        ids={n:n for n in c};words=normal_words(v,s)
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp)
            original=pack(v,s,c,ids,'test:metres',root/'original',corner_words=words)
            shared=pack_shared_only(v,s,c,ids,'test:metres',root/'shared',corner_words=words)
            self.assertTrue(original['passed']and shared['passed'])
            self.assertEqual((root/'original/cells.bin').read_bytes(),(root/'shared/cells.bin').read_bytes())
            self.assertEqual(original['parts'],shared['parts'])
            self.assertEqual(original['checks']['cell_volumes_m3'],shared['checks']['cell_volumes_m3'])
            self.assertEqual(sorted(p.name for p in (root/'shared').glob('*.npz')),['arrangement.npz'])
            self.assertEqual(len(list((root/'original').glob('*.npz'))),4)

if __name__=='__main__':unittest.main()
