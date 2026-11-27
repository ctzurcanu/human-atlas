"""Independent v4 source-bound exterior and fixed-context gate; no anatomy mutation."""
from pathlib import Path
from collections import Counter
import hashlib,json,time,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.distance import point_triangle
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1';C=B/'candidate-v4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,default=lambda v:v.tolist() if hasattr(v,'tolist') else v.item())+'\n')
def trikey(t):return tuple(sorted(map(tuple,t)))
def main():
 start=time.monotonic();plan=json.loads((B/'plan.json').read_text());proof=json.loads((C/'proof.json').read_text());contract=json.loads((B/'local-completion-contract-v4.json').read_text());cp=ROOT/plan['source_catalogue'];assert sha(cp)==plan['source_master_sha256']==contract['source_master_sha256']==proof['source_master_sha256'];assert sha(C/'cohort.npz')==proof['candidate_sha256']
 m=Mesh.load(C/'cohort.npz');an={k:v for k,v in np.load(C/'ancestry.npz').items()};course={k:v for k,v in np.load(B/'source-course-arrays.npz').items()};rows={int(x['source']):x for x in plan['sources']};src={}
 for sid in plan['selected_whole_source_indices']:
  row=rows[sid];assert sha(row['path'])==row['sha256'];src[sid]=Mesh.load(row['path'])
 expected=np.concatenate([src[s].vertices for s in plan['selected_whole_source_indices']]);assert np.array_equal(an['source_position_prefix'],expected)and np.array_equal(m.vertices[:len(expected)],expected)
 for f in ['normals','uvs']:assert np.array_equal(m.attributes[f][:len(expected)],np.concatenate([src[s].attributes[f]for s in plan['selected_whole_source_indices']]))
 # Independent retained-port set: native graph external ports, less the two now-declared c1 roots.
 retained=[p for p in plan['fixed_terminal_ports']if not p['root_matches']and not(p['source']in[0,1]and p['component']==1 and p['end']==0)];assert len(retained)==22
 expected_edges=set();portproof=[]
 for p in retained:
  rr=course[p['ring_array']][p['ring_station']];sm=src[p['source']];v=sm.vertices[rr];expected_edges|={frozenset([tuple(v[k]),tuple(v[(k+1)%len(v)])])for k in range(len(v))}
  matches=np.array([np.flatnonzero(np.all(m.vertices==x,axis=1))[0]for x in v]);assert np.array_equal(m.attributes['normals'][matches],sm.attributes['normals'][rr]);assert np.array_equal(m.attributes['uvs'][matches],sm.attributes['uvs'][rr]);portproof.append(dict(source=p['source'],component=p['component'],end=p['end'],position_normal_UV_literal=True))
 vv,inv=np.unique(m.vertices,axis=0,return_inverse=True);ff=inv[m.faces];e=np.concatenate([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]);ue,cnt=np.unique(np.sort(e,axis=1),axis=0,return_counts=True);be=ue[cnt==1];actual={frozenset(map(tuple,vv[x]))for x in be};missing=expected_edges-actual;extra=actual-expected_edges
 # Every new patch must be an explicitly declared real cone, never a flat classifier disk.
 cones={}
 for a in contract['anchors']:
  sid=a['source'];ci=a['component'];rings=course[f's{sid}_c{ci}_rings'];station=4 if sid in[0,1]and ci==1 and a['end']==0 else(0 if a['end']==0 else -1);v=src[sid].vertices[rings[station]];anchor=np.array(a['anchor_m'])
  for k in range(len(v)):cones[(sid,trikey(np.array([v[k],v[(k+1)%len(v)],anchor])))]=True
 orig=an['source_face_ids']>=0;new=~orig;assert set(np.unique(an['face_kinds']))=={'original','junction_faces'};assert np.all(an['source_face_ids'][new]==-2);assert np.all(an['face_kinds'][new]=='junction_faces');maxerr=0.;baryerr=0.;cone_matches=0
 for j in range(len(m.faces)):
  sid=int(an['source_indices'][j]);sf=int(an['source_face_ids'][j]);pt=an['parent_triangles'][j];bc=an['corner_barycentric'][j]
  if sf>=0:assert np.array_equal(pt,src[sid].vertices[src[sid].faces[sf]])
  else:assert (sid,trikey(pt))in cones;cone_matches+=1
  maxerr=max(maxerr,float(np.max(np.abs(bc@pt-m.vertices[m.faces[j]]))));baryerr=max(baryerr,float(max(0,-bc.min(),bc.max()-1,np.max(np.abs(bc.sum(1)-1)))))
 assert maxerr<1e-12 and baryerr<1e-10
 for trim in contract['trimmed_source_faces']:assert not np.any((an['source_indices']==trim['source'])&np.isin(an['source_face_ids'],trim['trimmed_source_face_ids']))
 t=topology(m);tree=BVH(m);pairs,summary=pair_contacts(m,m,0.,first_bvh=tree,second_bvh=tree,selfcheck=True,max_pairs=600000,max_tests=2000000);assert summary['complete_within_pair_budget'];strict=pairs['pairs'][np.isin(pairs['relation'],['crossing','coplanar-overlap'])]
 report=dict(schema='maxillary-v4-independent-exterior-gate/v1',complete=True,source_master_sha256=sha(cp),candidate_sha256=sha(C/'cohort.npz'),contract_sha256=sha(B/'local-completion-contract-v4.json'),triangles=len(m.faces),topology=t,retained_external_profiles=len(retained),actual_boundary_edges=len(actual),missing_external_edges=len(missing),unexpected_external_edges=len(extra),original_position_normal_UV_prefix_literal=True,external_port_proof=portproof,real_cone_output_faces=cone_matches,maximum_parent_reconstruction_error_m=maxerr,maximum_barycentric_outside=baryerr,strict_crossings=len(strict),self_summary=summary,elapsed_seconds=time.monotonic()-start,installable=False)
 report['passed']=not(missing or extra or len(strict))and len(t['components'])==2 and sum(x['euler']for x in t['components'])==-18 and not any(t[k]for k in['nonmanifold_edges','inconsistent_winding_edges','duplicate_faces','zero_area_faces'])
 write(C/'independent-exterior-audit.json',report);print(json.dumps({k:report[k]for k in['complete','passed','triangles','retained_external_profiles','actual_boundary_edges','missing_external_edges','unexpected_external_edges','real_cone_output_faces','strict_crossings','elapsed_seconds']}))
if __name__=='__main__':main()
