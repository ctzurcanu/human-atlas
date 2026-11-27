"""Bounded native GG repair stages; no catalogue or skeleton mutation.

Stage A uses existing global-guarded corner-plane constraints to untangle
only the seven witnessed native basal controls. Stage B is separately saved.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.vertex_clearance_lp import repair_vertex_clearance

BASE=ROOT/'.local-models/geometry-workbench/tongue-completion'
SOURCE=BASE/'source-v1/meshes'
CATALOGUE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='cbed9740aca9cd20ca12bf43a051eb7b36f60f90dc1919b64818dde8848f963f'
START=time.monotonic();DEADLINE=START+50


def check():
    if time.monotonic()>=DEADLINE:raise TimeoutError('GG stage exceeded its 50 second native deadline')
    if hashlib.sha256(CATALOGUE.read_bytes()).hexdigest()!=EXPECTED:raise ValueError('Authoritative male catalogue changed; do not use stale source exports')


def normal_words(v,f):
    n=np.zeros_like(v);t=v[f];cross=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    for i in range(3):np.add.at(n,f[:,i],cross)
    length=np.linalg.norm(n,axis=1);n[length>0]/=length[length>0,None]
    return np.rint(np.clip(n,-1,1)*32767).astype('<i2')


def main():
    check();out=BASE/'genioglossus-v2-stage-a';out.mkdir(exist_ok=False)
    rows=[]
    for side,index in [('left',1),('right',2)]:
        check();source=Mesh.load(SOURCE/f'native-{index:02d}.npz')
        known=np.array([[3009,3021],[3014,3021]],int)
        free=np.unique(source.faces[np.unique(known)])
        fixed=np.setdiff1d(np.arange(len(source.vertices)),free)
        t=source.triangles;fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);fn/=np.linalg.norm(fn,axis=1)[:,None]
        words=np.repeat(fn[:,None,:],3,axis=1)
        guard=LocalCrossingGuard(source.vertices,source.faces,known,
            maximum_displacement_m=.00035,maximum_pair_tests=180000)
        result=repair_vertex_clearance(source.vertices,source.vertices,source.faces,words,
            origin=[0,0,0],global_guard=guard,maximum_displacement_m=.00035,
            locked=fixed,wall_seconds=12.,maximum_lp_solves=24,maximum_pair_constraints=8)
        changed=np.flatnonzero(np.any(result['positions']!=source.vertices,axis=1))
        attrs={k:v.copy()for k,v in source.attributes.items()}
        touched=np.flatnonzero(np.any(np.isin(source.faces,changed),axis=1))
        shading=np.unique(source.faces[touched])if len(touched)else np.array([],int)
        attrs['normals'][shading]=normal_words(result['positions'],source.faces)[shading]
        candidate=Mesh(result['positions'],source.faces,source.frame,attrs,
            {'record':source.source['record'],'catalogue_sha256':EXPECTED,
             'source_geometry_sha256':source.geometry_sha,'repair':'Bounded native basal self-crossing repair',
             'anatomical_acceptance':False})
        path=out/f'genioglossus-{side}.npz';candidate.save(path)
        diag=diagnostics(candidate,self_intersections=True)
        row=dict(side=side,mesh=str(path),source_geometry_sha256=source.geometry_sha,
            mesh_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            free_source_controls=free.tolist(),changed_source_controls=changed.tolist(),
            displacements_m=np.linalg.norm(candidate.vertices[changed]-source.vertices[changed],axis=1).tolist(),
            upper_and_outside_patch_controls_literal=bool(np.array_equal(candidate.vertices[fixed],source.vertices[fixed])),
            source_faces_literal=bool(np.array_equal(candidate.faces,source.faces)),
            solver=result['evidence'],topology=topology(candidate),open3d=diag,
            numerical_self_repair_passed=bool(result['evidence']['passed']and diag['self_intersection_pairs']==0),
            not_yet_fixed=['Mandibular origin gap','Native GG hyoid-body interface crossings'],
            installed=False,anatomical_acceptance=False)
        rows.append(row)
    report=dict(schema='native-gg-local-untangling/v2-stage-a',source_catalogue_sha256=EXPECTED,
        elapsed_seconds=time.monotonic()-START,native_deadline_seconds=50,skeleton_modified=False,
        candidates=rows,installed=False,anatomical_acceptance=False)
    (out/'candidate.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
