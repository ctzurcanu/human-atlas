import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_root_junction import union_source_junction
from geometry_tools.algorithms.graft_collar_union import union_graft_collars
from geometry_tools.tests.test_source_root_collar import fixture


def fourway(paired=False,field_alias=False):
    m,cells,roots=fixture(4e-7 if paired else 1e-8)
    v=m.vertices.tolist();faces=m.faces.copy()
    if paired:
        # Match the tested source pattern: four literal shared points, with
        # only two planar chart coordinates differing (not four warped ones).
        for i in range(26,30):
            v[i][1]-=4e-7
            if i<28:v[i][0]+=4e-7
    if field_alias:
        alias=len(v);v.append(m.vertices[26].tolist());faces[cells[2]]=np.where(faces[cells[2]]==26,alias,faces[cells[2]])
    root=roots[2]if paired else roots[0]
    shift=np.array([-.0006,-.0006 if paired else .0006,0])
    target=list(range(len(v),len(v)+6));v.extend((np.array(v)[root]+shift).tolist())
    add=[]
    for i in range(6):
        j=(i+1)%6;add.extend([[root[i],root[j],target[j]],[root[i],target[j],target[i]]])
    add=np.array(add)
    if not paired:add=add[:,[0,2,1]]
    cells.append(list(range(len(faces),len(faces)+12)));roots.append(root)
    v=np.array(v);colors=np.tile(np.arange(len(v),dtype=np.uint8)[:,None],(1,3))
    return Mesh(v,np.r_[faces,add],'fixture-metres',dict(colors=colors)),cells,roots


class SourceRootJunctionTest(unittest.TestCase):
    def test_fourway_remaps_actual_corner_alias_and_retains_its_fields(self):
        m,cells,roots=fourway(field_alias=True)
        out,p=union_source_junction(m,cells,roots,roots[0],branch_count=4,maximum_correction_m=2e-7,maximum_root_area_relative_change=.005)
        np.testing.assert_array_equal(out.vertices[:len(m.vertices)],m.vertices)
        np.testing.assert_array_equal(out.attributes['colors'][:len(m.vertices)],m.attributes['colors'])
        # Raw26 is the declared representative; raw30 is used by the cell.
        edit=next(x for x in p['rootSeamEdits']if x['oldRawId']==30)
        np.testing.assert_array_equal(out.attributes['colors'][edit['newAliasRawId']],m.attributes['colors'][30])
        self.assertEqual(p['carrierAudit']['mouths'],4)
        self.assertTrue(p['carrierAudit']['topological_continuity_pass'])

    def test_two_literal_charts_form_one_four_port_graft_without_source_motion(self):
        m,cells,roots=fourway(paired=True)
        out,p=union_graft_collars(m,cells,[roots[0],roots[2]],junction_count=1)
        np.testing.assert_array_equal(out.vertices[:len(m.vertices)],m.vertices)
        self.assertEqual(p['sourceRootChartCount'],2)
        self.assertEqual(p['sourceJunctionCount'],1)
        self.assertEqual(p['carrierAudit']['mouths'],4)
        self.assertEqual(p['carrierAudit']['component_euler'],[-2])
        self.assertEqual(p['rootClosuresRetained'],0)
        self.assertEqual(p['emittedLumenDisks'],0)

    def test_branch_count_is_explicit_and_cannot_change_by_inference(self):
        m,cells,roots=fourway()
        with self.assertRaisesRegex(ValueError,'Explicit three/four-branch'):
            union_source_junction(m,cells,roots,roots[0],branch_count=3,maximum_correction_m=2e-7,maximum_root_area_relative_change=.005)


if __name__=='__main__':unittest.main()
