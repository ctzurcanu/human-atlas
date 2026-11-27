"""Actual installed lumen sections and one analytic body-support constraint; no fitting."""
from pathlib import Path
import json,time,hashlib,sys,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.distance import point_triangle
from geometry_tools.audit import topology,winding_contains
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1';O=B/'installed-body-sinus-landmarks-v1'
def sections(tri,point,normal):
 normal=np.asarray(normal);normal/=np.linalg.norm(normal);u=np.array([1.,0,0]);u-=normal*(u@normal)
 if np.linalg.norm(u)<.1:u=np.array([0.,0,1]);u-=normal*(u@normal)
 u/=np.linalg.norm(u);v=np.cross(normal,u);da=(tri-point)@normal;aa=tri;bb=np.roll(tri,-1,axis=1);db=np.roll(da,-1,axis=1);cross=da*db<0;ids=np.flatnonzero(cross.sum(1)==2);segments=[]
 for i in ids:
  use=np.flatnonzero(cross[i]);p=aa[i,use]+(bb[i,use]-aa[i,use])*(da[i,use]/(da[i,use]-db[i,use]))[:,None];segments.append(p)
 if not segments:return []
 pts=np.array(segments);xyz,inv=np.unique(np.round(pts.reshape(-1,3),11),axis=0,return_inverse=True);edges=set(tuple(sorted(map(int,e)))for e in inv.reshape(-1,2));adj={}
 for a,b in edges:adj.setdefault(a,set()).add(b);adj.setdefault(b,set()).add(a)
 if any(len(n)!=2 for n in adj.values()):return [dict(closed=False,degree_counts={str(k):sum(len(n)==k for n in adj.values())for k in set(map(len,adj.values()))})]
 todo=set(adj);out=[]
 while todo:
  seed=min(todo);order=[];prev=None;cur=seed
  while True:
   order.append(cur);todo.discard(cur);nxt=next(x for x in adj[cur]if x!=prev)
   if nxt==seed:break
   prev,cur=cur,nxt
  p=xyz[order];q=np.c_[(p-point)@u,(p-point)@v];nxt=np.roll(q,-1,axis=0);c=q[:,0]*nxt[:,1]-nxt[:,0]*q[:,1];area=c.sum()/2;cent=((q+nxt)*c[:,None]).sum(0)/(6*area);world=point+u*cent[0]+v*cent[1];out.append(dict(closed=True,points_m=p.tolist(),centre_m=world.tolist(),area_m2=float(abs(area)),equivalent_diameter_m=float(2*np.sqrt(abs(area)/np.pi))))
 return sorted(out,key=lambda x:x['centre_m'][2])
def closest_loops(a,b):
 pa=np.array(a['points_m']);pb=np.array(b['points_m']);best=(float('inf'),None,None)
 for p,qs,swap in[(p,pb,False)for p in pa]+[(p,pa,True)for p in pb]:
  q2=np.roll(qs,-1,axis=0);d=q2-qs;t=np.clip(np.sum((p-qs)*d,axis=1)/np.sum(d*d,axis=1),0,1);q=qs+t[:,None]*d;j=int(np.linalg.norm(q-p,axis=1).argmin());dist=float(np.linalg.norm(q[j]-p))
  if dist<best[0]:best=(dist,q[j]if swap else p,p if swap else q[j])
 return dict(distance_m=best[0],first_m=best[1].tolist(),second_m=best[2].tolist(),midpoint_m=((best[1]+best[2])/2).tolist())
def main():
 start=time.monotonic();O.mkdir(exist_ok=True);C=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';sha=hashlib.sha256(C.read_bytes()).hexdigest();assert sha=='a9ae0fbaaf62fdf600a966e43641f54a4b20621944784e31ec43f43260428638';cat=json.loads(C.read_text());parts={p['id']:p for p in cat['parts']};reader=AtlasReader();ids=['ZA:Left common carotid artery','TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall','TA98PART:ZA:Internal carotid artery.l:A12.2.06.002','ZA:External carotid artery.l'];meshes={id:reader.read_part(cat,parts[id],frame='model:local-ta98:assembled-metres',catalogue_path=C,catalogue_sha=sha)for id in ids};outer=[];inner=[]
 # Draw-group sign/chart order is NOT an outer/inner tissue classification.
 real=[]
 for id,m in meshes.items():
  for g in parts[id]['tissueDrawGroups']:
   if g.get('neighbor')not in ids:real.extend(m.triangles[g['start']//3:(g['start']+g['count'])//3])
 real=np.array(real);centroids=real.mean(1);outer_source=Mesh.load(B/'localized-no-rejoin-wall-v1/outer-union-classifier-diagnostic.npz');inner_source=Mesh.load(B/'localized-no-rejoin-wall-v1/inner-union-classifier-diagnostic.npz');od,_=BVH(outer_source).closest_points(centroids);ind,_=BVH(inner_source).closest_points(centroids);od=od['distance_m'];ind=ind['distance_m'];known=np.minimum(od,ind)<1e-6;outer=real[known&(od<ind)];inner=real[known&(ind<od)];surface_classification=dict(method='Actual installed triangles classified by closest source-bound outer/inner classifier surface; draw-group chart signs are not anatomical layer labels.',outer_source_geometry_sha256=outer_source.geometry_sha,inner_source_geometry_sha256=inner_source.geometry_sha,actual_outer_faces=len(outer),actual_inner_faces=len(inner),unclassified_actual_faces=int((~known).sum()),maximum_accepted_centroid_source_distance_m=float(np.minimum(od,ind)[known].max()));horizontal=[]
 for y in np.arange(1.60541,1.61442,.00025):
  horizontal.append(dict(Y_m=float(y),outer=sections(outer,np.array([0.,y,0]),np.array([0.,1,0])),inner=sections(inner,np.array([0.,y,0]),np.array([0.,1,0]))))
 transition={}
 for name in['outer','inner']:
  first=next(i for i,r in enumerate(horizontal)if len(r[name])==2 and all(x['closed']for x in r[name]));r=horizontal[first];transition[name]=dict(previous_plane_Y_m=horizontal[first-1]['Y_m'],first_separated_plane_Y_m=r['Y_m'],nearest_branch_contours=closest_loops(*r[name]))
 guide=np.load(B/'ica-localized-no-rejoin-v1/world-profiles.npz');arc=guide['source_arc_m'];centres=guide['target_centres'];normal=np.gradient(centres,arc,axis=0);luminal=[]
 for mm in[20,22,24,26,28,30,35,40,45,50,60]:
  s=mm/1000;p=np.array([np.interp(s,arc,centres[:,j])for j in range(3)]);n=np.array([np.interp(s,arc,normal[:,j])for j in range(3)]);loops=sections(inner,p,n);closed=[x for x in loops if x.get('closed')]
  if not closed:luminal.append(dict(source_arc_mm=mm,closed_loops=0));continue
  chosen=min(closed,key=lambda x:np.linalg.norm(np.array(x['centre_m'])-p));luminal.append(dict(source_arc_mm=mm,plane_origin_m=p.tolist(),plane_normal=(n/np.linalg.norm(n)).tolist(),closed_loops=len(closed),selected_inner_loop=chosen,selection_basis='Actual closed inner-surface section nearest the held ICA guide; merged proximal sections are not a distinct ICA lumen.'))
 firstY=transition['outer']['first_separated_plane_Y_m'];origin=np.array(transition['outer']['nearest_branch_contours']['midpoint_m']);before=[r for r in horizontal if r['Y_m']<firstY and len(r['outer'])==1 and r['outer'][0].get('closed')];cranial=np.array(before[-1]['outer'][0]['centre_m'])-np.array(before[-3]['outer'][0]['centre_m']);cranial/=np.linalg.norm(cranial);posterior=np.array([0.,0,-1]);posterior-=cranial*(posterior@cranial);posterior/=np.linalg.norm(posterior);medial=np.array([-1.,0,0]);medial-=cranial*(medial@cranial)+posterior*(medial@posterior);medial/=np.linalg.norm(medial);basis=np.stack([medial,cranial,posterior],axis=1);a=np.array([.0011,.003,.00105]);direction=(medial+posterior)/np.sqrt(2);cy=origin+.0015*cranial;pts=outer.reshape(-1,3);roi=pts[np.abs((pts-cy)@cranial)<.0036];support=float(np.max((roi-cy)@direction));shape_support=float(np.sqrt(np.sum((a*(basis.T@direction))**2)));margin=.0001;query=cy+(support+shape_support+margin)*direction;constraints=[]
 for id in['arterial-outer-surfaces','ZA:Glossopharyngeal nerve (IX).l','ZA:Internal jugular vein.l']:
  if id=='arterial-outer-surfaces':tri=outer;gh=None
  else:m=reader.read_part(cat,parts[id],frame='model:local-ta98:assembled-metres',catalogue_path=C,catalogue_sha=sha);tri=m.triangles;gh=m.geometry_sha
  transformed=(tri-query)@basis/a;d,_,_=point_triangle(np.zeros((len(transformed),3)),transformed);normdist=float(np.sqrt(d.min()));constraints.append(dict(id=id,geometry_sha256=gh,minimum_ellipsoid_normalized_surface_radius=normdist,conservative_physical_surface_gap_lower_bound_m=float((normdist-1)*a.min()),no_surface_contact_bound=normdist>1))
 report=dict(schema='installed-left-carotid-body-and-sinus-landmark-evidence/v1',complete=True,source_catalogue_sha256=sha,actual_bank_sha256=next(iter(meshes.values())).source['bank_sha256'],actual_part_geometry_sha256={id:m.geometry_sha for id,m in meshes.items()},surface_classification=surface_classification,horizontal_sections=horizontal,actual_flow_divider_brackets=transition,ICA_axis_normal_inner_sections=luminal,body_landmark_frame=dict(origin_at_outer_crotch_m=origin.tolist(),basis_columns_medial_cranial_posterior=basis.tolist(),derivation='Flow-divider nearest actual outer contours at first separated plane. Cranial direction from preceding common outer sections; posterior and medial are projected source anatomical axes.'),body_size_constraint=dict(illustrative_diameters_m=(2*a).tolist(),analytic_ellipsoid_volume_mm3=float(4*np.pi*np.prod(a)/3*1e9),primary_2025_volume_range_mm3=[10,25]),one_analytic_clearance_constraint=dict(query_centre_m=query.tolist(),not_installed_or_fitted=True,posteromedial_unit_direction=direction.tolist(),outer_support_from_reference_centre_m=support,shape_support_m=shape_support,margin_m=margin,rule='For this6mm longitudinal interval, body centre projection must exceed actual arterial outer support by ellipsoid support plus0.1mm. This is a conservative one-shot geometric constraint, not an anatomically accepted placement.',surface_checks=constraints),sinus_interpretation='Measured proximal inner-lumen calibre may delimit a provisional immaterial reference region. A generic taper, union bulge or selected medial wall patch does not establish the histological baroreceptor field or a discrete sinus bulb. No extra anatomical solid is justified by these measurements alone.',geometry_mutated=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start);(O/'evidence.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(complete=True,flow_brackets=transition,body_query=query.tolist(),surface_checks=constraints,ICA_sections=[dict(arc_mm=r['source_arc_mm'],loops=r['closed_loops'],diameter_mm=r.get('selected_inner_loop',{}).get('equivalent_diameter_m',0)*1000)for r in luminal],elapsed=report['elapsed_seconds'])))
if __name__=='__main__':main()
