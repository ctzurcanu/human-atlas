"""Replay the rejected HG producer and distinguish stored-fragment defects."""
from __future__ import annotations
import sys,importlib.util,time
from pathlib import Path
import numpy as np
from fractions import Fraction as Q

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts,self_contact_summary
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.exact_source_intersections import edge_plane,contains

s=importlib.util.spec_from_file_location('hg',Path(__file__).with_name('repair-tongue-hyoglossus-origin.py'));hg=importlib.util.module_from_spec(s);s.loader.exec_module(hg)

def exact_intersection(first,second):
    rational=lambda p:tuple(Q.from_float(float(x))for x in p)
    a=tuple(rational(p)for p in first);b=tuple(rational(p)for p in second);hits=set()
    for source,target in [(a,b),(b,a)]:
        hits.update(p for p in source if contains(p,target))
        for i in range(3):
            q=edge_plane(source[i],source[(i+1)%3],target)
            if q is not None and contains(q,target):hits.add(q)
    shared=set(a)&set(b);array=np.array([[float(x)for x in p]for p in hits]).reshape(-1,3)
    length=max((float(np.linalg.norm(a-b))for a in array for b in array),default=0.)
    return dict(exactRationalIntersectionPoints=array.tolist(),sharedActualSourceCorners=len(shared),
        allHitsAreSharedCorners=bool(hits and hits<=shared),intersectionSpanM=length)

def main(output):
    began=time.monotonic();output.mkdir(exist_ok=False);rows=[]
    bones=[Mesh.load(hg.WORK/f'source-v1/meshes/native-{i:02}.npz')for i in [13,14,15,16,17]];v=[];f=[];offset=0
    for bone in bones:v.append(bone.vertices);f.append(bone.faces+offset);offset+=len(bone.vertices)
    hyoid=Mesh(np.vstack(v),np.vstack(f),bones[0].frame);make_solid(hyoid.vertices,hyoid.faces);tree=BVH(hyoid)
    for side in ['left','right']:
        source=hg.WORK/f'mucosal-insertions-v2-certified-rim/HG-{side}.npz';muscle=Mesh.load(source)
        footprint=Mesh.load(hg.WORK/f'native-attachment-candidate-v1/hyoglossus-{side}-bone-origin-surface.npz');footprint,_=hg.adjacent_native_footprint(footprint,hyoid,tree)
        sheet,_=hg.broad_sheet(footprint,muscle)
        native,oid=make_solid(muscle.vertices,muscle.faces,properties=muscle.attributes['normals'][muscle.faces]);addition,_=make_solid(sheet.vertices,sheet.faces)
        replay=native.trim_by_plane([0.,1.,0.],1.6188)+addition;raw=replay.to_mesh64();tri=np.asarray(raw.vert_properties,float)[np.asarray(raw.tri_verts,int),:3]
        float64=Mesh(tri.reshape(-1,3),np.arange(tri.size//3).reshape(-1,3),muscle.frame)
        old_path=hg.WORK/f'hyoglossus-origin-v4-replaced-lower/HG-{side}.npz';old=Mesh.load(old_path)
        encoded=(float64.vertices-[0.,1.625,0.]).astype('<f4').astype(float)+[0.,1.625,0.]
        same=np.array_equal(encoded,old.vertices)and np.array_equal(float64.faces,old.faces)
        if not same:raise ValueError('Replay does not reproduce the exact rejected stored triangle bank')
        before,_=self_contact_summary(float64,max_pairs=100000,max_tests=1000000)
        all_pairs,summary=pair_contacts(old,old,0.,selfcheck=True,max_pairs=100000,max_tests=1000000)
        mask=all_pairs['relation']=='crossing';pairs=all_pairs['pairs'][mask];proof=[]
        for first,second in pairs:
            actual=old.triangles[[first,second]];raw_pair=float64.triangles[[first,second]]
            p=dict(faces=[int(first),int(second)],actualStoredCorners=actual.tolist(),float64ProducerCorners=raw_pair.tolist(),
                storedExactIntersection=exact_intersection(*actual),float64ExactIntersection=exact_intersection(*raw_pair))
            proof.append(p)
        float64.save(output/f'HG-{side}-producer-float64.npz');atomic_npz(output/f'{side}-stored-self-witnesses.npz',**{k:v[mask]for k,v in all_pairs.items()})
        row=dict(side=side,rejectedCandidate=str(old_path),rejectedCandidateSha256=digest_file(old_path),
            exactStoredReplay=True,float64ProducerSelf=before,actualStoredSelf=summary,
            maximumEncodingErrorM=float(np.linalg.norm(encoded-float64.vertices,axis=1).max()),pairs=proof)
        rows.append(row);atomic_json(output/f'{side}-precision.json',row)
    atomic_json(output/'summary.json',dict(schema='tongue-hg-union-fragment-precision/v1',rows=rows,geometryModified=False,
        candidateStillRejected=True,elapsedSeconds=time.monotonic()-began))
    print([dict(side=r['side'],float64Crossings=r['float64ProducerSelf']['crossing_face_pairs'],storedCrossings=r['actualStoredSelf']['crossing_face_pairs'],maxEncodingErrorM=r['maximumEncodingErrorM'],exactSharedOnly=sum(p['storedExactIntersection']['allHitsAreSharedCorners']for p in r['pairs']))for r in rows],flush=True)

if __name__=='__main__':main(hg.WORK/'hyoglossus-origin-v4-precision-audit')
