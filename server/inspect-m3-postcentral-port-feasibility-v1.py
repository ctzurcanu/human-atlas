"""Local source-chart/port feasibility only; no walls or geometry are invented."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.open3d_backend import SurfaceIndex

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3'
OUT=B/'m3-postcentral-port-feasibility-v1';FRAME='model:local-ta98:assembled-metres'

def chart(mesh):
 v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];groups=recover_station_rings(v,f);result=[]
 for index,g in enumerate(groups):
  rings,phase=align_station_correspondence(v,f[g['face_ids']],g['rings']);rings=np.array(rings);centers=v[rings].mean(1);stations=g['station_of_vertex'];mouths=[]
  for station in [0,len(rings)-1]:
   ring=v[rings[station]];center=centers[station];_,_,basis=np.linalg.svd(ring-center,full_matrices=False);normal=basis[-1];direction=centers[min(station+1,len(centers)-1)]-centers[max(0,station-1)]
   if normal@direction<0:normal=-normal
   mouths.append(dict(station=station,center=center.tolist(),normal=normal.tolist(),controls=len(ring),area_m2=projected_area(ring),radii_m=np.linalg.norm(ring-center,axis=1).tolist(),maximum_planarity_error_m=float(abs((ring-center)@normal).max())))
  result.append(dict(index=index,source_component=g['source_component'],face_ids=g['face_ids'],rings=rings,centers=centers,station_of_vertex=stations,mouths=mouths,phase=phase))
 return v,f,result

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);source=B/'source-atlas.json';a=json.loads(source.read_bytes());ps={p['id']:p for p in a['parts']};reader=AtlasReader();ids=['ZA:Middle cerebral artery (M3 segment).l','ZA:Postcentral arterial branch.l'];meshes=[reader.read_part(a,ps[pid],frame=FRAME)for pid in ids];charts=[chart(m)for m in meshes];parent,branch=meshes;pv,pf,pg=charts[0];bv,bf,bg=charts[1];assert len(pg)==1 and len(bg)==2
 arrays,summary=pair_contacts(branch,parent,0.);assert summary['complete_within_pair_budget'];cross=arrays['pairs'][np.isin(arrays['relation'],['crossing','coplanar-overlap'])];assert len(cross)==48
 label=np.full(len(bf),-1,int)
 for group in bg:label[group['face_ids']]=group['index']
 roots={int(k):int(np.sum(label[cross[:,0]]==k))for k in np.unique(label[cross[:,0]])};assert len(roots)==1
 root_component=next(iter(roots));root=bg[root_component];root_faces=np.unique(cross[:,0]);parent_faces=np.unique(cross[:,1]);parent_stations=pg[0]['station_of_vertex'][pf[parent_faces]];branch_stations=root['station_of_vertex'][bf[root_faces]];branch_stations=branch_stations[branch_stations>=0];parent_range=[int(parent_stations.min()),int(parent_stations.max())];branch_range=[int(branch_stations.min()),int(branch_stations.max())]
 ports=[]
 for group in bg:
  for mouth in group['mouths']:
   near,_=SurfaceIndex(parent).closest_points(np.array(mouth['center'])[None],certify=True);ports.append(dict(component=group['index'],**mouth,parent_actual_surface_distance_m=float(near['distance_m'][0])))
 other=bg[1-root_component];first=Mesh(bv,bf[root['face_ids']],FRAME);second=Mesh(bv,bf[other['face_ids']],FRAME);twig_arrays,twig_summary=pair_contacts(first,second,0.);assert twig_summary['complete_within_pair_budget']
 # A bounded proposed local patch consists of exact station intervals around
 # the observed cross-pairs. No unrelated source face belongs to this patch.
 p0=max(0,parent_range[0]-2);p1=min(len(pg[0]['rings'])-1,parent_range[1]+2);b0=max(0,branch_range[0]-2);b1=min(len(root['rings'])-1,branch_range[1]+2)
 p_mask=np.all((pg[0]['station_of_vertex'][pf]>=p0)&(pg[0]['station_of_vertex'][pf]<=p1),axis=1);b_mask=np.zeros(len(bf),bool);owned=bf[root['face_ids']];inside=np.all((root['station_of_vertex'][owned]>=b0)&(root['station_of_vertex'][owned]<=b1),axis=1);b_mask[root['face_ids'][inside]]=True
 packet=dict(parent_vertices=pv,parent_faces=pf,parent_station_rings=pg[0]['rings'],parent_local_patch_face_ids=np.flatnonzero(p_mask),branch_vertices=bv,branch_faces=bf,branch_root_component_faces=root['face_ids'],branch_root_station_rings=root['rings'],branch_local_patch_face_ids=np.flatnonzero(b_mask),branch_other_component_faces=other['face_ids'],actual_crossing_face_pairs=cross,frame=np.array(FRAME));atomic_npz(OUT/'source-port-charts.npz',**packet)
 components=[]
 for role,(v,f,groups)in zip(['M3','postcentral'],charts):
  components.append(dict(role=role,components=[dict(index=g['index'],faces=len(g['face_ids']),station_count=len(g['rings']),ring_size=len(g['rings'][0]),mouths=g['mouths'])for g in groups]))
 proof=dict(schema='m3-postcentral-source-port-feasibility/v1',source_catalogue_sha256=digest_file(source),frame=FRAME,installed=False,geometry_modified=False,source_parts=[dict(id=pid,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,bank=m.source['bank_path'],bank_sha256=m.source['bank_sha256'])for pid,m in zip(ids,meshes)],source_components=components,root_source_component=root_component,root_contact_pairs_by_component=roots,root_contacts=summary,other_postcentral_component_contacts=twig_summary,parent_contact_station_range=parent_range,branch_contact_station_range=branch_range,explicit_local_patch=dict(parent_station_interval=[p0,p1],branch_station_interval=[b0,b1],parent_patch_faces=int(p_mask.sum()),branch_patch_faces=int(b_mask.sum()),all_source_faces_outside_patch_fixed=True,other_branch_component_fixed=True),source_chart_archive=str(OUT/'source-port-charts.npz'),source_chart_archive_sha256=digest_file(OUT/'source-port-charts.npz'),outer_inner_correspondence_exists=False,source_finite_wall=False,finite_wall_construction_permitted_now=False,missing_inputs=['No native inner wall/chart or measured arterial thickness exists in these outside casts.','Root overlap identifies the adjoining component but does not establish a shared annular tissue interface or open blood volume.','A local engineering inner profile and thickness contract would be newly constructed, not recovered source anatomy.'],algorithm_plan=[dict(module='tube_station_topology',use='Recovered literal serial outer rings and source face IDs; preserve phase and all protected outer faces.'),dict(module='geometric_root_collar',use='Not applicable to current single outer sheet: requires paired outer/inner annular outlet profiles.'),dict(module='branch_junction_envelopes.close_branch',use='After explicit paired inner profile exists, make root-connected outer/fluid classifiers with internal anchor strictly inside parent fluid.'),dict(module='communicating_vascular_network.boolean_union',use='Union declared outer classifiers and separately union blood classifiers; subtract blood union so no membrane remains. Remove every classifier-only closure from installed tissue.'),dict(module='tube_rim_bridge.annular_mouths',use='Do not use for a lateral root aperture; no paired annular wall source exists.')],remaining_review='Source root component is numerically identified; no reconstruction/Boolean/offset was run. Actual reference/sections, inner thickness and source-face-preserving Boolean patch provenance must precede any candidate.',elapsed_seconds=time.monotonic()-start)
 atomic_json(OUT/'feasibility.json',proof);print(json.dumps({k:proof[k]for k in ['root_source_component','root_contact_pairs_by_component','parent_contact_station_range','branch_contact_station_range','explicit_local_patch','outer_inner_correspondence_exists','finite_wall_construction_permitted_now','elapsed_seconds']}))

if __name__=='__main__':main()
