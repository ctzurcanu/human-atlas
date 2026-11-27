"""Measure flagged GG contact seams without changing candidate geometry."""
from pathlib import Path
import importlib.util
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import winding_contains

BASE=ROOT/'.local-models/geometry-workbench/tongue-completion';OUT=BASE/'genioglossus-v2-stage-b'
SOURCE=BASE/'source-v1/meshes';START=time.monotonic()


def main():
    spec=importlib.util.spec_from_file_location('local_gg_attachment',ROOT/'server/attach-tongue-genioglossus.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    bones=[Mesh.load(SOURCE/f'native-{i:02d}.npz')for i in [5,13,14,15,16,17]]
    bv,bf,_,_,_,_=module.local_mandible_classifier(bones[0]);jaw=Mesh(bv,bf,bones[0].frame)
    vv=[];ff=[];count=0
    for bone in bones[1:]:vv.append(bone.vertices);ff.append(bone.faces+count);count+=len(bone.vertices)
    hyoid=Mesh(np.vstack(vv),np.vstack(ff),bones[0].frame)
    report=[];tolerance=1e-10
    for side in ['left','right']:
        mesh=Mesh.load(OUT/f'genioglossus-{side}.npz');points=np.unique(mesh.triangles.reshape(-1,3),axis=0);tests=[]
        for bone in bones:
            name=bone.source['record']['id'];z=np.load(OUT/f'{side}-{name.replace(":","_")}-contacts.npz')
            selected=z['relation']=='crossing';pairs=z['pairs'][selected]
            if not len(pairs):continue
            a=mesh.triangles[pairs[:,0]];b=bone.triangles[pairs[:,1]]
            nb=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0]);nb/=np.linalg.norm(nb,axis=1)[:,None]
            na=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);na/=np.linalg.norm(na,axis=1)[:,None]
            ha=np.einsum('tki,ti->tk',a-b[:,0,None,:],nb);hb=np.einsum('tki,ti->tk',b-a[:,0,None,:],na)
            excursion=np.minimum(np.maximum(-ha.min(1),0),np.maximum(ha.max(1),0))
            other=np.minimum(np.maximum(-hb.min(1),0),np.maximum(hb.max(1),0))
            separation=np.minimum(excursion,other)
            true=np.flatnonzero(separation>tolerance)
            tests.append(dict(bone_part_id=name,strict_crossing_face_pairs=len(pairs),
                one_plane_minimum_two_sided_excursion_max_m=float(separation.max()),
                pairs_straddling_both_planes_beyond_declared_tolerance=len(true),
                pair_examples_beyond_tolerance=pairs[true[:10]].tolist(),
                declared_plane_tolerance_m=tolerance,
                interpretation='Separating-plane excursions only; no geometry or acceptance gate changed.'))
        containment=[]
        for name,bone in [('actual assembled hyoid',hyoid),('local witnessed mandible classifier',jaw)]:
            result=winding_contains(bone,points,boundary_tolerance_m=tolerance,distance_backend='open3d')
            inside=np.flatnonzero(result['inside']);boundary=np.flatnonzero(result['boundary_ambiguous'])
            containment.append(dict(compartment=name,unique_candidate_vertices=len(points),
                contained_vertices_outside_boundary_tolerance=len(inside),ambiguous_interface_vertices=len(boundary),
                contained_examples=points[inside[:5]].tolist(),boundary_tolerance_m=tolerance))
        report.append(dict(side=side,geometry_sha256=mesh.geometry_sha,flagged_contact_plane_tests=tests,containment=containment))
    result=dict(schema='native-gg-shared-contact-numerics/v1',elapsed_seconds=time.monotonic()-START,
        candidate_modified=False,original_conservative_gate_not_overwritten=True,anatomical_acceptance=False,results=report)
    (OUT/'interface-numerics.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
