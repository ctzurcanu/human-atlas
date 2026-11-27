"""Source ownership survives pose changes and near-coincident planar ears."""
import unittest
import tempfile
from pathlib import Path
import json
import numpy as np

from geometry_tools.boolean_cells import triangulate_overlay,BooleanCells
from geometry_tools.urinary_wall_layers import outlet_membership


class UrinaryOwnershipTests(unittest.TestCase):
    def test_literal_adjacent_cells_keep_triangles_and_share_the_archived_source_bank(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
        a=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
        b=np.array([[0,1,2],[0,4,1],[1,4,2],[2,4,0]])
        graph=BooleanCells();solids={'upper':graph.literal_cell(p,a,'upper'),'lower':graph.literal_cell(p,b,'lower')}
        for name,faces in [('upper',a),('lower',b)]:
            mesh=solids[name].to_mesh64()
            np.testing.assert_array_equal(mesh.vert_properties[mesh.tri_verts],p[faces])
        with tempfile.TemporaryDirectory()as td:
            path=Path(td)/'cells.npz';graph.save(solids,path)
            with np.load(path)as z:
                meta=json.loads(str(z['manifest_json']))
                self.assertEqual(len({q['v_key']for q in meta['primitives']}),1)
            restored,cells=BooleanCells.load(path)
            primitives=list(restored.primitives.values());self.assertIs(primitives[0][0],primitives[1][0])
            assembled=restored.assemble(cells)
            self.assertEqual(set(assembled['cells']),{'upper','lower'})

    def test_outlet_source_ownership_is_pose_independent(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
        t=np.array([[0,1,2,3]])
        # Two different source owners at the same height cannot be separated
        # by the old vertical-coordinate rule.
        u,e=outlet_membership(p,f,t,[tuple(p[0])],[tuple(p[1])])
        self.assertEqual(u[0],0.);self.assertEqual(u[1],1.)
        self.assertTrue(np.all((u>=0)&(u<=1)))
        rotation=np.array([[0.,0,1],[0,1,0],[-1,0,0]])
        q=p@rotation.T+[4.,-2,3]
        v,_=outlet_membership(q,f,t,[tuple(q[0])],[tuple(q[1])])
        np.testing.assert_allclose(u,v,atol=1e-14)
        with self.assertRaises(ValueError):
            outlet_membership(p,f,t,[tuple(p[0])],[tuple(p[0])])

    def test_near_coincident_ear_cannot_fill_the_contact_notch(self):
        # A real local source-plane contour: the sixth and seventh controls
        # differ by 1.5 femtometres in physical coordinates. Zero-uncertainty
        # triangulation produced a reversed 3.65 mm2 face across the notch.
        polygon=np.array([
            [9755.918820746243,834.0652193799615],[6188.578508548439,13024.249194979668],
            [4874.29523563385,17515.369607046247],[4139.427169062197,20026.533708408475],
            [3806.4400763958693,21164.40494184196],[3758.2417823821306,21329.106377631426],
            [3758.24178237468,21329.106377631426],[8464.133679024875,5248.305845841765],
            [-129.92715580761433,1125.3380414098501],[-2434.937323048711,21089.722013883293],
            [-2312.185734167695,20026.533708408475],[-2022.2564885020256,17515.369607046247],
            [-1503.729183755815,13024.249194979668],[-96.29792802035809,834.0652193799615],
            [0.,0.],[10000.,0.]])
        factor=4977625.405892116;expected=1.2650964247064879e-6
        flat,t=triangulate_overlay([polygon],expected,factor,2e-12)
        q=flat[t]/factor;a=q[:,1]-q[:,0];b=q[:,2]-q[:,0]
        areas=(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2
        self.assertLess(abs(areas).sum()-expected,1e-16)
        self.assertGreater(areas.min(),-1e-16)
        np.testing.assert_array_equal(flat,polygon)


if __name__=='__main__':unittest.main()
