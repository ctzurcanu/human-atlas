"""Complete intersection segments and stable parent exit; no fit/union/install."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict,Counter
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import audit,winding_contains
from geometry_tools.distance import exact_intersection_segment
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';I=P/'temporal-source-incidence-v2';C=P/'three-course-finite-educational-v1';OUT=I/'common-parent-root-zone-v1';FRAME='model:local-ta98:assembled-metres';BAND=1e-8

def coplanar_interval(first,last,triangle):
 p=[[Q(float(x))for x in q]for q in [first,last,*triangle]];a,b,x,y,z=p;sub=lambda u,v:[r-s for r,s in zip(u,v)];u=sub(y,x);w=sub(z,x);n=[u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]];drop=max(range(3),key=lambda i:abs(n[i]));keep=[i for i in range(3)if i!=drop];a,b,x,y,z=[[q[i]for i in keep]for q in p];cross=lambda u,v,w:(v[0]-u[0])*(w[1]-u[1])-(v[1]-u[1])*(w[0]-u[0]);sign=1 if cross(x,y,z)>0 else-1;lo,hi=Q(0),Q(1)
 for u,v in [(x,y),(y,z),(z,x)]:
  f0=sign*cross(u,v,a);f1=sign*cross(u,v,b);d=f1-f0
  if not d:
   if f0<0:return None
  elif d>0:lo=max(lo,-f0/d)
  else:hi=min(hi,-f0/d)
 return(lo,hi)if lo<=hi else None

def cuts(first,last,mesh,known_boundary=False):
 if known_boundary:return[Q(0),Q(1)],[(Q(0),Q(1))],0
 tri=mesh.triangles;ids=np.flatnonzero(np.all(tri.max(1)>=np.minimum(first,last)-1e-12,axis=1)&np.all(tri.min(1)<=np.maximum(first,last)+1e-12,axis=1));values={Q(0),Q(1)};coplanar=[];tests=0
 for fid in ids:
  tests+=1;r=exact_segment_triangle_intersection(first,last,tri[fid])
  if not r['intersects']:continue
  if r['coplanar']:
   interval=coplanar_interval(first,last,tri[fid])
   if interval is not None:values.update(interval);coplanar.append(interval)
  else:values.add(Q(int(r['parameter_numerator']),int(r['parameter_denominator'])))
 return sorted(values),coplanar,tests

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);g=json.loads((C/'generation-complete.json').read_text());paired=np.load(C/'M3-derived-paired-chart.npz');reference={};reference_audits={}
 for kind,key in [('outer','outer_positions'),('fluid','inner_positions')]:
  path=next(r['checkpoint_input']for r in g['native_inputs']if r['role']=='M3-'+kind);q=np.load(C/path);v=q['positions_mm']/1000.+np.array(g['origin_m']);v[:len(paired[key])]=paired[key];mesh=Mesh(v,q['faces'],FRAME);reference[kind]=mesh;mesh.save(OUT/('unchanged-M3-'+kind+'-classifier.npz'));result=audit(mesh,require_closed=True,expected_components=1,expected_genus=[0]);atomic_json(OUT/(kind+'-reference-audit.json'),result);assert result['passed'];reference_audits[kind]=result
 z=np.load(F/'three-course-source-charts.npz');tz=np.load(I/'complete-source-charts.npz');tv=tz['temporal_vertices'];tf=tz['temporal_faces'];labels=np.full(len(tv),-1,int)
 for station,ring in enumerate(tz['component0_rings']):labels[ring]=station
 segments=[];errors=[]
 for name,pv,pf,pr,fullids in [('M3',z['M3_vertices'],z['M3_faces'],z['M3_rings'],np.arange(len(z['M3_faces']))),('postcentral',z['branch_vertices'],z['branch_faces'][z['main_source_face_ids']],z['main_rings'],z['main_source_face_ids'])]:
  ls=np.full(len(pv),-1,int)
  for station,ring in enumerate(pr):ls[ring]=station
  pairs=np.load(I/('M3-contact-pairs.npz'if name=='M3'else'main-postcentral-contact-pairs.npz'))['pairs']
  for a,b in pairs:
   r=exact_intersection_segment(tv[tf[a]],pv[pf[b]]);assert r['relation']=='crossing'and len(r['endpoints'])==2;segments.append(dict(family='source-temporal-'+name,points=r['endpoints'],temporal_face=int(a),other_face=int(fullids[b]),temporal_stations=sorted(set(map(int,labels[tf[a]]))),other_stations=sorted(set(map(int,ls[pf[b]]))),outer_boundary=name=='M3'));errors.extend(r['coordinate_error_upper_m'].ravel())
 witnesses=json.loads((I/'inner-crossing-witnesses-v1/crossing-witnesses.json').read_text())
 for row in witnesses['witnesses']:
  if not row['new_peer_face']:continue
  r=exact_intersection_segment(np.array(row['wall_triangle_m']),np.array(row['temporal_triangle_m']));assert r['relation']=='crossing';segments.append(dict(family='new-'+row['role'],points=r['endpoints'],temporal_face=row['temporal_source_face'],other_face=row['source_face'],temporal_stations=row['temporal_stations'],other_stations=row['source_stations'],outer_boundary=False));errors.extend(r['coordinate_error_upper_m'].ravel())
 # Centerline intervals are included in the same complete piecewise analysis.
 centers=tv[tz['component0_rings']].mean(1)
 for station,(first,last)in enumerate(zip(centers[:-1],centers[1:])):segments.append(dict(family='temporal-centerline',points=np.array([first,last]),temporal_stations=[station,station+1],other_stations=[],temporal_face=None,other_face=None,outer_boundary=False))
 samples=[];pieces=[];tests=0
 for index,s in enumerate(segments):
  first,last=s['points'];oc,op,ot=cuts(first,last,reference['outer'],s['outer_boundary']);fc,fp,ft=cuts(first,last,reference['fluid']);tests+=ot+ft;parameters=sorted(set(oc+fc))
  for low,high in zip(parameters[:-1],parameters[1:]):
   if low==high:continue
   mid=(low+high)/2;point=first+float(mid)*(last-first);samples.append(point);pieces.append(dict(segment=index,low=float(low),high=float(high),length_m=float(high-low)*float(np.linalg.norm(last-first)),outer_exact_boundary=any(a<=mid<=b for a,b in op),fluid_exact_boundary=any(a<=mid<=b for a,b in fp)))
 points=np.array(samples);outer=winding_contains(reference['outer'],points,boundary_tolerance_m=BAND);fluid=winding_contains(reference['fluid'],points,boundary_tolerance_m=BAND)
 for i,piece in enumerate(pieces):
  ob=piece['outer_exact_boundary']or outer['boundary_ambiguous'][i];fb=piece['fluid_exact_boundary']or fluid['boundary_ambiguous'][i];oi=outer['inside'][i];fi=fluid['inside'][i];piece['classification']='parent-outer-boundary'if ob else'outside-parent-outer'if not oi else'parent-fluid-boundary'if fb else'common-parent-fluid'if fi else'common-parent-wall';piece['fluid_inside']=fi;piece['outer_inside']=oi
 grouped=defaultdict(list)
 for piece in pieces:grouped[segments[piece['segment']]['family']].append(piece)
 summaries=[]
 for family,rows in grouped.items():
  counts=Counter(r['classification']for r in rows);lengths={key:sum(r['length_m']for r in rows if r['classification']==key)for key in counts};outside=[r for r in rows if r['classification']=='outside-parent-outer'];summaries.append(dict(family=family,segments=len({r['segment']for r in rows}),piece_classes=dict(counts),length_by_class_m=lengths,outside_temporal_stations=sorted(set(k for r in outside for k in segments[r['segment']]['temporal_stations'])),outside_other_stations=sorted(set(k for r in outside for k in segments[r['segment']]['other_stations'])),outside_temporal_faces=sorted(set(segments[r['segment']]['temporal_face']for r in outside if segments[r['segment']]['temporal_face']is not None))))
 # Contact-curve components use exact rounded rational endpoint identity;
 # no proximity graph is used to invent a common ostial connection.
 contact=[i for i,s in enumerate(segments)if s['family'].startswith('source-')];endpoint={};adj=defaultdict(set)
 for index in contact:
  for p in segments[index]['points']:
   key=tuple(p)
   for other in endpoint.get(key,[]):adj[index].add(other);adj[other].add(index)
   endpoint.setdefault(key,[]).append(index)
 seen=set();patches=[]
 for index in contact:
  if index in seen:continue
  todo=[index];members=[]
  while todo:
   k=todo.pop()
   if k in seen:continue
   seen.add(k);members.append(k);todo.extend(adj[k]-seen)
  owned=[r for r in pieces if r['segment']in set(members)];patches.append(dict(segments=members,families=sorted(set(segments[k]['family']for k in members)),classes=dict(Counter(r['classification']for r in owned)),outside_parent_length_m=sum(r['length_m']for r in owned if r['classification']=='outside-parent-outer'),temporal_stations=sorted(set(j for k in members for j in segments[k]['temporal_stations']))))
 center=grouped['temporal-centerline'];notoutside=[p for p in center if p['classification']!='outside-parent-outer'];last=max(notoutside,key=lambda p:segments[p['segment']]['temporal_stations'][0]+p['high']);laststation=segments[last['segment']]['temporal_stations'][0];stableexit=laststation+last['high'];outerexit=centers[laststation]+last['high']*(centers[laststation+1]-centers[laststation]);atomic_npz(OUT/'actual-segments-and-pieces.npz',segment_endpoints=np.array([s['points']for s in segments]),piece_midpoints=points,piece_segment_ids=np.array([p['segment']for p in pieces]),frame=np.array(FRAME))
 result=dict(schema='actual-temporal-common-parent-root-zone/v1',no_geometry_modified=True,no_fit_or_union_run=True,reference_classifiers='Unchanged complete native M3 outer and authorized derived M3 fluid, independently closed genus0/self-crossing-free.',classification_boundary_band_m=BAND,boundary_band_purpose='Numerical classification band covering <=3.73nm packed inner-coordinate quantization; not anatomical thickness or clinical clearance.',actual_intersection_segment_coordinate_error_upper_m=max(errors,default=0.),exact_segment_parent_triangle_tests=tests,piecewise_method='Exact binary-rational segment/triangle crossings split each complete contact segment at every parent boundary; each resulting open interval is classified by full winding, with coplanar boundary intervals explicit.',families=summaries,contact_curve_components=patches,temporal_centerline_stable_outer_exit_station_parameter=stableexit,temporal_centerline_stable_outer_exit_point_m=outerexit.tolist(),all_later_centerline_pieces_outside_parent=True,segments=[dict(**{k:v for k,v in s.items()if k!='points'},endpoints_m=s['points'].tolist())for s in segments],pieces=pieces,common_parent_zone_acceptance=False,decision='Read actual outside-parent intervals and graph; no endpoint-only or face-centroid-only reroute conclusion.',elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'containment-proof.json',result);print(json.dumps(dict(families=summaries,contact_curve_components=[{k:p[k]for k in ['families','classes','outside_parent_length_m','temporal_stations']}for p in patches],stable_outer_exit_station_parameter=stableexit,exact_segment_parent_triangle_tests=tests,elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
