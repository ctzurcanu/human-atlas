import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.current_source_regions import (
    gather_current_regions,recover_current_sections,scatter_region_pose,scatter_whole_constraints)


def fixture():
    angle=np.arange(6)*2*np.pi/6
    p=np.array([[np.cos(a),np.sin(a),z]for z in range(5)for a in angle],float)
    f=[]
    for z in range(4):
        for i in range(6):
            a=z*6+i;b=z*6+(i+1)%6;c=b+6;d=a+6;f.extend([[a,b,c],[a,c,d]])
    f=np.array(f);records={};meshes={}
    for pid,faceids in [('proximal',np.arange(24)),('distal',np.arange(24,48))]:
        used=np.unique(f[faceids])[::-1];lookup=np.full(len(p),-1,int);lookup[used]=np.arange(len(used))
        attrs=dict(normals=np.tile(np.array([20,30,40],np.int16),(len(used),1)),normalSigns=np.ones(len(used),np.int8),uvs=np.c_[used,used*2].astype(np.float32))
        meshes[pid]=Mesh(p[used],lookup[f[faceids]],'fixture:m',attrs,dict(source_vertex_ids=used.tolist()))
        records[pid]=dict(sourceDerivedRegion=dict(sourcePartId='tube-source',sourceId='tube',sourceFaceIds=faceids.tolist(),boundaryContractSha256='same-partition'))
    return records,meshes


class CurrentRegionTests(unittest.TestCase):
    def gather(self,records,meshes):
        return gather_current_regions(records,meshes,family_id='tube-source',expected_owner_ids=['proximal','distal'])

    def test_gather_scatter_exact_current_values_and_ownership(self):
        records,meshes=fixture();assembly=self.gather(records,meshes)
        whole=assembly.whole;attrs={k:v.copy()for k,v in whole.attributes.items()};attrs['normals'][:,0]+=5
        posed=Mesh(whole.vertices+np.array([.1,.2,.3]),whole.faces,whole.frame,attrs,whole.source)
        parts=scatter_region_pose(assembly,posed)
        for pid,old in meshes.items():
            np.testing.assert_array_equal(parts[pid].vertices,old.vertices+[.1,.2,.3])
            np.testing.assert_array_equal(parts[pid].faces,old.faces)
            np.testing.assert_array_equal(parts[pid].attributes['normals'],old.attributes['normals']+[5,0,0])
            np.testing.assert_array_equal(parts[pid].attributes['normalSigns'],old.attributes['normalSigns'])
            np.testing.assert_array_equal(parts[pid].attributes['uvs'],old.attributes['uvs'])
            self.assertEqual(parts[pid].source,old.source)
        shared=set(assembly.local_to_whole['proximal'])&set(assembly.local_to_whole['distal'])
        scatter=scatter_whole_constraints(assembly,np.array(sorted(shared)),posed.vertices[sorted(shared)])
        for pid in meshes:self.assertEqual(len(scatter[pid]['rawIds']),6)

    def test_virtual_chart_cut_cancels_before_real_section_recovery(self):
        records,meshes=fixture();assembly=self.gather(records,meshes);guide=recover_current_sections(assembly)
        self.assertEqual(len(guide['boundary_raw']),2)
        self.assertTrue(len(guide['profiles'])>=3)
        self.assertTrue(all(len(ids)==6 for ids in guide['boundary_raw']))
        # The shared chart rim is an interior true tube profile, not a mouth.
        shared=set(assembly.local_to_whole['proximal'])&set(assembly.local_to_whole['distal'])
        self.assertTrue(any(set(ids)==shared for ids in guide['profiles']))
        self.assertFalse(any(set(ids)==shared for ids in guide['boundary_raw']))

    def test_rejects_mixed_ancestry_shared_drift_and_duplicate_face_ownership(self):
        records,meshes=fixture();records['distal']['sourceDerivedRegion']['sourcePartId']='another-source'
        with self.assertRaises(ValueError):self.gather(records,meshes)
        records,meshes=fixture();old=meshes['distal'];p=old.vertices.copy();p[-1]+=[0,0,.01]
        meshes['distal']=Mesh(p,old.faces,old.frame,old.attributes,old.source)
        with self.assertRaises(ValueError):self.gather(records,meshes)
        records,meshes=fixture();old=meshes['distal'];attrs={k:v.copy()for k,v in old.attributes.items()};attrs['normals'][-1,0]+=1
        meshes['distal']=Mesh(old.vertices,old.faces,old.frame,attrs,old.source)
        with self.assertRaises(ValueError):self.gather(records,meshes)
        records,meshes=fixture();records['distal']['sourceDerivedRegion']['sourceFaceIds'][0]=0
        with self.assertRaises(ValueError):self.gather(records,meshes)


if __name__=='__main__':unittest.main()
