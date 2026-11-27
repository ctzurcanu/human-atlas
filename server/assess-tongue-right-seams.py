"""Map actual review pixels to visible surfaces, without changing geometry."""
from pathlib import Path
import json
import sys
import time
from collections import defaultdict
import numpy as np
import open3d as o3d

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file


def rays(mesh,target_yz):
    tri=mesh.triangles;yz=tri[:,:,[1,2]]
    matrix=np.stack([yz[:,1]-yz[:,0],yz[:,2]-yz[:,0]],axis=-1)
    valid=np.abs(np.linalg.det(matrix))>1e-18
    ids=np.flatnonzero(valid);inverse=np.linalg.inv(matrix[valid]);base=yz[valid,0]
    result=[]
    for p in target_yz:
        uv=np.einsum('nij,nj->ni',inverse,p-base);bary=np.c_[1-uv.sum(1),uv]
        hit=np.all(bary>=-1e-12,axis=1)
        if not np.any(hit):result.append(None);continue
        face_ids=ids[hit];weights=bary[hit];x=np.einsum('ni,ni->n',weights,tri[face_ids,:,0])
        chosen=int(np.argmax(x));face=int(face_ids[chosen]);point=weights[chosen]@tri[face]
        normal=np.cross(tri[face,1]-tri[face,0],tri[face,2]-tri[face,0]);normal/=np.linalg.norm(normal)
        artist=None
        if 'normals'in mesh.attributes:
            field=mesh.attributes['normals'].astype(float);field/=np.linalg.norm(field,axis=1)[:,None]
            artist=weights[chosen]@field[mesh.faces[face]];artist/=np.linalg.norm(artist)
        result.append(dict(x=float(x[chosen]),face=face,point=point.tolist(),barycentric=weights[chosen].tolist(),
            geometric_normal=normal.tolist(),presentation_normal=artist.tolist()if artist is not None else None))
    return result


def main():
    began=time.monotonic();work=ROOT/'.local-models/geometry-workbench/tongue-completion'
    candidate=work/'macro-envelopes-v1-artist-normals-final';record=json.loads((candidate/'construction.json').read_text())
    queue=ROOT/'.local-models/geometry-workbench/review/065041ab456bef60079e52f8bf11a784f33f50097ac852b6943582a65f1bdb2d'
    task=json.loads((queue/'request.json').read_text());render=json.loads((queue/'result.json').read_text())
    assert render['state']=='rendered'and task['views'][0]['direction']==[1,0,0]and not task['views'][0].get('section')
    meshes=[(r['key'],Mesh.load(r['mesh']))for r in record['meshes']]
    points=np.concatenate([m.vertices for _,m in meshes]);center=(points.min(0)+points.max(0))/2
    scale=float(np.linalg.norm(np.ptp(points,axis=0)))*1.15;resolution=task['resolution']
    selected=[(256,273),(329,340),(343,365),(323,389),(355,381),(302,340),(357,324),(315,409)]
    grid=[(x,y)for y in range(305,416,10)for x in range(290,371,10)]
    pixels=selected+grid
    # Exact camera basis for the review's right view: screen right=-Z, up=+Y.
    targets=np.array([[center[1]+(.5-(y+.5)/resolution)*scale,
                       center[2]-((x+.5)/resolution-.5)*scale]for x,y in pixels])
    sources={'native-inner-envelope':Mesh.load(work/'intrinsic-territories-v2/native-inner-envelope.npz'),
             'native-outer-tongue':Mesh.load(work/'source-v1/meshes/native-00.npz')}
    hits={name:rays(mesh,targets)for name,mesh in meshes+list(sources.items())}
    image=np.asarray(o3d.io.read_image(str(queue/'view-0.png')))[:,:,:3]
    witnesses=[];deficits=[]
    for index,pixel in enumerate(pixels):
        available=[(name,hits[name][index])for name,_ in meshes if hits[name][index]is not None]
        front=max(available,key=lambda r:r[1]['x'])if available else None
        inner=hits['native-inner-envelope'][index];outer=hits['native-outer-tongue'][index]
        deficit=None
        if inner is not None and front is not None:
            deficit=inner['x']-front[1]['x'];deficits.append(deficit)
        if index<len(selected):
            witnesses.append(dict(pixel=list(pixel),rgb=image[pixel[1],pixel[0]].tolist(),
                visible_part=front[0]if front else None,visible_hit=front[1]if front else None,
                native_inner_hit=inner,native_outer_hit=outer,inner_to_assembled_front_deficit_m=deficit))
    tv=dict(meshes)['transverse-vertical-shared'];tri=tv.triangles;n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    incident=defaultdict(list)
    for face,t in enumerate(tv.faces):
        for k in range(3):incident[tuple(sorted((int(t[k]),int(t[(k+1)%3]))))].append(face)
    angles=[]
    for edge,faces in incident.items():
        if len(faces)==2:
            angle=float(np.degrees(np.arccos(np.clip(n[faces[0]]@n[faces[1]],-1,1))))
            if angle>45:
                angles.append(dict(edge=list(edge),faces=faces,angle_degrees=angle,center=tv.vertices[list(edge)].mean(0).tolist()))
    report=dict(schema='tongue-right-render-geometry-visibility-assessment/v1',
        construction_sha256=digest_file(candidate/'construction.json'),render_result_sha256=digest_file(queue/'result.json'),
        render_path=str(queue/'view-0.png'),no_section=True,no_hidden_face_configuration=True,
        renderer_cavity_enabled=True,renderer_cavity_scope='Workbench BOTH highlights screen/world-space creases',
        meshes=[dict(key=name,geometry_sha256=m.geometry_sha,triangles=len(m.faces))for name,m in meshes],
        exact_camera=dict(center=center.tolist(),ortho_scale=scale,direction=[1,0,0],screen_right=[0,0,-1],screen_up=[0,1,0]),
        witness_pixels=witnesses,grid_samples=len(grid),
        maximum_inner_to_assembled_front_deficit_m=max(deficits,default=None),
        p95_inner_to_assembled_front_deficit_m=float(np.percentile(deficits,95))if deficits else None,
        tv_edges_above_45_degrees=len(angles),tv_crease_examples=sorted(angles,key=lambda x:-x['angle_degrees'])[:20],
        geometry_modified=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(candidate/'right-seam-assessment.json',report)
    print(json.dumps(dict(witnesses=[dict(pixel=r['pixel'],visible=r['visible_part'],
        deficit_m=r['inner_to_assembled_front_deficit_m'])for r in witnesses],
        maximum_inner_front_deficit_m=report['maximum_inner_to_assembled_front_deficit_m'],
        tv_sharp_edges=len(angles),elapsed_seconds=report['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
