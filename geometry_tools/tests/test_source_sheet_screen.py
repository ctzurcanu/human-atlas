"""Avoid confusing a hollow source envelope with an axial discontinuity."""
import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_sheet_screen import screen_source_components


def tube(z0,z1,radius):
    theta=np.arange(12)*2*np.pi/12
    ring=np.c_[radius*np.cos(theta),radius*np.sin(theta),np.zeros(12)]
    v=np.r_[ring+[0,0,z0],ring+[0,0,z1]]
    faces=np.array([[i,(i+1)%12,(i+1)%12+12] for i in range(12)]+[[i,(i+1)%12+12,i+12] for i in range(12)])
    normals=np.tile(np.c_[np.cos(theta),np.sin(theta),np.zeros(12)],(2,1))
    return v,faces,normals

def pair(*,nested=True,opposed=True):
    a=tube(0,.02,.004)
    b=tube(0,.02,.00398) if nested else tube(.020001,.04,.004)
    v=np.r_[a[0],b[0]];f=np.r_[a[1],b[1][:,::-1]+24]
    normals=np.r_[a[2],b[2]*(-1 if opposed else 1)]
    return Mesh(v,f,'test:metres',dict(normals=np.rint(normals*32767).astype('<i2'),normalSigns=np.ones(48,'i1')))

class SourceSheetScreen(unittest.TestCase):
    def test_complete_nested_opposed_sheets_are_not_a_mouth_to_weld(self):
        r=screen_source_components(pair());self.assertTrue(r['pairedOpposedSourceSheets']);self.assertFalse(r['automaticWeldAllowed'])

    def test_two_axial_fragments_with_close_ends_are_not_coextensive_sheets(self):
        r=screen_source_components(pair(nested=False));self.assertFalse(r['pairedOpposedSourceSheets']);self.assertFalse(r['automaticWeldAllowed'])

    def test_near_identical_same_facing_copies_are_not_classified_as_opposed_sheets(self):
        r=screen_source_components(pair(opposed=False));self.assertFalse(r['pairedOpposedSourceSheets'])

if __name__=='__main__':unittest.main()
