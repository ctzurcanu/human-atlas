"""Repair only the declared collinear conus-tip defect and audit native roots."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit,topology
from geometry_tools.spatial import BVH

ap=argparse.ArgumentParser();ap.add_argument('--plan',required=True);ap.add_argument('--output',required=True)
args=ap.parse_args();started=time.monotonic();plan_path=Path(args.plan).resolve();plan=json.loads(plan_path.read_text())
source=plan['source_catalogues']['female'];path=Path(source['path']);assert digest_file(path)==source['sha256']
a=json.loads(path.read_text());reader=AtlasReader();part=next(p for p in a['parts']if p['id']=='LOCAL:female:spinal_cord:864')
m=reader.read_part(a,part,frame='model:local-female-ta98:assembled-metres',catalogue_path=path)
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True);assert not any(out.iterdir()),'A new or empty output stage is required'
v=m.vertices.copy();f=m.faces.copy();ids=np.flatnonzero(v[:,1]==v[:,1].min())
tip=np.array([0.,v[ids[0],1],v[ids[0],2]])
assert len(np.unique(v[ids],axis=0))==5 and np.max(np.ptp(v[ids,1:],axis=0))<1e-15
movement=np.linalg.norm(v[ids]-tip,axis=1);assert movement.max()<=.000321
v[ids]=tip
area=np.linalg.norm(np.cross(v[f][:,1]-v[f][:,0],v[f][:,2]-v[f][:,0]),axis=1)
keep=area>0;f=f[keep]
attrs={k:values.copy()for k,values in m.attributes.items()}
# Recompute normals only where the repaired terminal fan meets its first ring.
unique,inv=np.unique(v,axis=0,return_inverse=True);uf=inv[f]
tip_id=int(inv[ids[0]]);fan=np.flatnonzero(np.any(uf==tip_id,axis=1));affected=np.unique(uf[fan])
q=unique[uf];fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);sums=np.zeros_like(unique)
for k in range(3):np.add.at(sums,uf[:,k],fn)
sums/=np.maximum(np.linalg.norm(sums,axis=1)[:,None],1e-100)
if 'normals'in attrs:
    changed=np.isin(inv,affected);attrs['normals'][changed]=np.rint(sums[inv[changed]]*32767).astype('<i2')
repaired=Mesh(v,f,m.frame,attrs)
proof=audit(repaired,expected_components=1,expected_genus=[0],require_closed=True,max_tests=500000)
assert proof['passed'],proof['failures']
atomic_npz(out/'surface.npz',positions=v,faces=f,retained_source_faces=np.flatnonzero(keep),source_positions=m.vertices,
    source_faces=m.faces,**attrs)
atomic_json(out/'surface.json',dict(schema='native-conus-tip-repair/v1',source_catalogue_sha256=source['sha256'],
    plan_sha256=digest_file(plan_path),source_part=part,source_chunk=a['chunks'][part['chunk']],
    source_geometry_sha256=m.geometry_sha,candidate_geometry_sha256=repaired.geometry_sha,
    source_topology=topology(m),candidate_audit=proof,removed_faces=np.flatnonzero(~keep).tolist(),
    maximum_displacement_m=float(movement.max()),nonterminal_positions_literal=bool(np.array_equal(v[np.setdiff1d(np.arange(len(v)),ids)],m.vertices[np.setdiff1d(np.arange(len(v)),ids)])),
    unchanged_source_catalogue=digest_file(path)==source['sha256'],new_tissues=0,installed=False,anatomical_acceptance=False))
# Each named nerve's closest actual source vertex establishes a reproducible
# geometric attachment clue, never a declaration that a whole nerve is a root.
bvh=BVH(repaired);attachments=[]
for p in a['parts']:
    if p.get('suppressed')or p['id'].startswith('HRA:'):continue
    if not ('cervical_spinal_nerve_'in p['id']or'intercostal_nerve'in p['id']or'cauda_equina_base'in p['id']):continue
    nerve=reader.read_part(a,p,frame=m.frame,catalogue_path=path)
    vertices=np.unique(nerve.vertices[np.unique(nerve.faces)],axis=0)
    near=vertices[np.abs(vertices[:,0])<=.013]
    if not len(near):near=vertices[[np.argmin(np.abs(vertices[:,0]))]]
    # Finite native sources have at most a few thousand near-cord controls.
    result,_=bvh.closest_points(near);j=int(np.argmin(result['distance_m']))
    face=int(result['face_ids'][j]);point=np.asarray(result['barycentric'][j])@repaired.triangles[face]
    attachments.append(dict(part_id=p['id'],name=p['name'],distance_m=float(result['distance_m'][j]),
        nerve_point_m=near[j].tolist(),cord_point_m=point.tolist(),cord_face=face,
        cord_barycentric=result['barycentric'][j].tolist(),source_geometry_sha256=nerve.geometry_sha,
        anatomical_interface_verified=False))
atomic_json(out/'native-root-attachments.json',dict(source_catalogue_sha256=source['sha256'],
    cord_geometry_sha256=repaired.geometry_sha,attachments=attachments,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(triangles_before=len(m.faces),triangles_after=len(f),removed_zero_area_or_collapsed_faces=int(sum(~keep)),
    maximum_tip_movement_mm=float(movement.max()*1000),self_crossings=proof['self_intersections']['defect_pairs'],
    root_sources=len(attachments),elapsed_seconds=time.monotonic()-started)),flush=True)
