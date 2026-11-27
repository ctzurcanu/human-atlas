import unittest
import numpy as np
from geometry_tools.algorithms.nonsingular_source_annuli import recover_nonsingular_source_annuli
from geometry_tools.algorithms.source_tube_topology import recover_regular_source_annuli


def singular_grid():
    angle=np.arange(6)*np.pi/3
    ring=np.c_[np.cos(angle),np.zeros(6),np.sin(angle)]*.01
    vertices=np.concatenate([ring+[0,y*.03,0]for y in range(7)]).tolist()
    faces=[]
    def strip(first,last):
        for i in range(6):
            j=(i+1)%6;faces.extend([[first[i],first[j],last[j]],[first[i],last[j],last[i]]])
    for y in range(6):strip(np.arange(6)+y*6,np.arange(6)+(y+1)*6)
    previous=np.arange(6)+18
    for z in range(1,4):
        current=np.arange(len(vertices),len(vertices)+6)
        vertices.extend((ring+[0,.09,z*.03]).tolist());strip(previous,current);previous=current
    return np.array(vertices),np.array(faces)


class NonsingularAnnuliTest(unittest.TestCase):
    def test_singularities_remain_core_and_never_become_certified_annuli(self):
        v,f=singular_grid();before=v.copy();r=recover_nonsingular_source_annuli(v,f,np.arange(len(f)).reshape(-1,2))
        self.assertEqual(len(r['singular_source_edges']),6)
        self.assertGreater(len(r['strips']),0)
        forbidden=set(r['singular_source_vertices']);physical=np.array(r['raw_to_physical'])
        owned=[]
        for s in r['strips']:
            self.assertFalse(set(physical[f[s['source_face_ids']]].ravel())&forbidden)
            owned.extend(s['source_face_ids'])
        owned+=r['unassigned_source_face_ids']
        self.assertEqual(sorted(owned),list(range(len(f))))
        self.assertEqual(len(set(owned)),len(owned));np.testing.assert_array_equal(v,before)
        self.assertFalse(r['core_anatomical_identity_established'])

    def test_strict_whole_carrier_classifier_still_rejects_nonmanifold_input(self):
        v,f=singular_grid()
        with self.assertRaisesRegex(ValueError,'Nonmanifold'):
            recover_regular_source_annuli(v,f,np.arange(len(f)).reshape(-1,2))

    def test_unpaired_surface_remains_explicitly_owned_core(self):
        v,f=singular_grid();pairs=np.arange(len(f)).reshape(-1,2)
        v=np.r_[v,[[1,0,0],[1,1,0],[1,0,1]]];f=np.r_[f,[[len(v)-3,len(v)-2,len(v)-1]]]
        r=recover_nonsingular_source_annuli(v,f,pairs)
        self.assertIn(len(f)-1,r['unassigned_source_face_ids'])


if __name__=='__main__':unittest.main()
