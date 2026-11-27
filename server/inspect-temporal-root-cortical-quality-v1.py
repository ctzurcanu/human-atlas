"""Fixed cortical source quality and sampled signed clearance only; no fit."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit,winding_contains
from geometry_tools.open3d_backend import SurfaceIndex
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';I=P/'temporal-source-incidence-v2';C=P/'three-course-finite-educational-v1';CTX=C/'current-context-v1';FRAME='model:local-ta98:assembled-metres';IDS=['ZA:Insula (Subcentral gyrus and ant. and post. sulci*).l:Insula','ZA:Transverse temporal gyri.l:Temporal lobe']

def main(index):
 start=time.monotonic();out=I/('cortical-quality-'+str(index)+'-v1');out.mkdir(exist_ok=False);a=json.loads((CTX/'source-atlas.json').read_bytes());part=next(p for p in a['parts']if p['id']==IDS[index]);reader=AtlasReader();mesh=reader.read_part(a,part,frame=FRAME);mesh.save(out/'fixed-cortical-source.npz');quality=audit(mesh,require_closed=True);atomic_json(out/'actual-source-audit.json',quality);t=mesh.triangles;fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n=mesh.attributes['normals'].astype(float)
 if 'normalSigns'in mesh.attributes:n*=mesh.attributes['normalSigns'][:,None]
 corner=n[mesh.faces];denom=np.linalg.norm(corner,axis=2)*np.linalg.norm(fn,axis=1)[:,None];cos=np.divide(np.einsum('tki,ti->tk',corner,fn),denom,out=np.zeros_like(denom),where=denom>0);normalbasis=dict(nonpositive_corners=int(np.sum(cos<=0)),minimum_corner_geometric_normal_cosine=float(cos.min()),oriented_volume_m3=quality['topology']['oriented_triangle_integral_m3'],source_attribute_basis='Actual Int16 source normal words times source normalSigns, compared to literal triangle orientation.')
 z=np.load(F/'three-course-source-charts.npz');tz=np.load(I/'complete-source-charts.npz');groups=[('M3-centers24-41',z['M3_vertices'][z['M3_rings'][24:42]].mean(1)),('M3-outer24-41',z['M3_vertices'][z['M3_rings'][24:42]].reshape(-1,3)),('main-centers167-180',z['branch_vertices'][z['main_rings'][167:181]].mean(1)),('main-outer167-180',z['branch_vertices'][z['main_rings'][167:181]].reshape(-1,3)),('temporal-centers0-29',tz['temporal_vertices'][tz['component0_rings'][:30]].mean(1)),('temporal-outer0-29',tz['temporal_vertices'][tz['component0_rings'][:30]].reshape(-1,3))];rows=[];inside_valid=quality['passed']and normalbasis['oriented_volume_m3']>0 and len(quality['topology']['components'])==1
 atomic_json(out/'native-nearest-query-pending.json',dict(source=str(out/'fixed-cortical-source.npz'),query_groups=[dict(name=name,points=len(points))for name,points in groups]));surface=SurfaceIndex(mesh)
 for name,points in groups:
  nearest,evidence=surface.closest_points(points,certify=True);winding=winding_contains(mesh,points,boundary_tolerance_m=1e-10)if inside_valid else None;inside=np.array(winding['inside'])if winding else np.zeros(len(points),bool);amb=np.array(winding['boundary_ambiguous'])if winding else np.zeros(len(points),bool);row=dict(query_group=name,sampled_points=len(points),unsigned_distance_min_m=float(nearest['distance_m'].min()),unsigned_distance_median_m=float(np.median(nearest['distance_m'])),inside_sample_count=int(inside.sum())if inside_valid else None,ambiguous_sample_count=int(amb.sum())if inside_valid else None,signed_depth_available=bool(inside_valid),scope='Distances at explicit source stations/outer controls, not a continuous penetration-depth bound.')
  if inside.any():
   witness=int(np.flatnonzero(inside)[np.argmax(nearest['distance_m'][inside])]);fid=int(nearest['face_ids'][witness]);normal=fn[fid]/np.linalg.norm(fn[fid]);row.update(maximum_sampled_inside_depth_m=float(nearest['distance_m'][witness]),inside_witness_m=points[witness].tolist(),nearest_source_surface_point_m=nearest['points'][witness].tolist(),source_face_id=fid,nearest_outward_geometric_normal=normal.tolist())
  rows.append(row);atomic_npz(out/(name+'-samples.npz'),points=points,nearest_points=nearest['points'],source_face_ids=nearest['face_ids'],distance_m=nearest['distance_m'],inside=inside,ambiguous=amb,frame=np.array(FRAME))
 result=dict(schema='fixed-cortical-source-quality-and-local-station-clearance/v1',part_id=IDS[index],geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,source_faces=len(mesh.faces),topology_and_self_audit_passed=quality['passed'],normal_basis=normalbasis,signed_depth_usable=bool(inside_valid),samples=rows,no_source_geometry_modified=True,no_fit_run=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start);atomic_json(out/'quality-depth-proof.json',result);print(json.dumps(result))

if __name__=='__main__':main(int(sys.argv[1]))
