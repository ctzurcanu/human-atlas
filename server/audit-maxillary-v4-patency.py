"""Diagnostic closure and source-course patency, never an atlas cap."""
from pathlib import Path
import hashlib,json,time,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import topology,winding_contains
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1';C=B/'candidate-v4'
def main():
 start=time.monotonic();p=json.loads((B/'plan.json').read_text());c=json.loads((B/'local-completion-contract-v4.json').read_text());m=Mesh.load(C/'cohort.npz');a={k:v for k,v in np.load(B/'source-course-arrays.npz').items()};sources={r['source']:Mesh.load(r['path'])for r in p['sources']if r['source']in p['selected_whole_source_indices']};retained=[x for x in p['fixed_terminal_ports']if not x['root_matches']and not(x['source']in[0,1]and x['component']==1 and x['end']==0)]
 vertices=m.vertices.tolist();faces=m.faces.tolist();ed={}
 for f in m.faces:
  for u,v in zip(f,np.roll(f,-1)):ed[(tuple(m.vertices[u]),tuple(m.vertices[v]))]=(int(u),int(v))
 for port in retained:
  rr=a[port['ring_array']][port['ring_station']];ring=sources[port['source']].vertices[rr];centre=len(vertices);vertices.append(ring.mean(0).tolist())
  for u,v in zip(ring,np.roll(ring,-1,axis=0)):
   pair=ed.get((tuple(u),tuple(v)))
   if pair is None:pair=ed[(tuple(v),tuple(u))]
   faces.append([pair[1],pair[0],centre])
 closed=Mesh(np.array(vertices),np.array(faces),m.frame,{},dict(diagnostic_only=True));t=topology(closed);assert t['closed_oriented_manifold'];assert len(t['components'])==2 and all(x['closed_orientable_genus']==0 for x in t['components'])
 probes=[];owners=[]
 for sid in p['selected_whole_source_indices']:
  row=next(r for r in p['sources']if r['source']==sid)
  for cr in row['courses']:
   rings=a[cr['prefix']+'_rings'];v=sources[sid].vertices[rings];lo=4 if sid in[0,1]and cr['component']==1 else 1
   for j in range(lo,len(v)-1):probes.append(v[j].mean(0));owners.append(dict(kind='native-profile-centre',source=sid,component=cr['component'],station=j))
 for j,an in enumerate(c['anchors']):
  sid=an['source'];ci=an['component'];rr=a[f's{sid}_c{ci}_rings'];st=4 if sid in[0,1]and ci==1 and an['end']==0 else(0 if an['end']==0 else -1);centre=sources[sid].vertices[rr[st]].mean(0);anchor=np.array(an['anchor_m'])
  for f in [.1,.25,.5,.75,.9]:probes.append(centre*(1-f)+anchor*f);owners.append(dict(kind='declared-junction-axis',anchor=j,fraction=f))
 w=winding_contains(closed,np.array(probes),boundary_tolerance_m=1e-10);bad=[dict(**owners[j],point_m=np.array(probes[j]).tolist(),winding=w['winding_number'][j],ambiguous=w['boundary_ambiguous'][j])for j in range(len(probes))if not w['inside'][j]or w['boundary_ambiguous'][j]]
 report=dict(schema='maxillary-v4-independent-source-axis-patency/v1',complete=True,passed=not bad,candidate_sha256=hashlib.sha256((C/'cohort.npz').read_bytes()).hexdigest(),diagnostic_closure_only_external_mouths=22,closed_topology=t,probes=len(probes),failures=bad,elapsed_seconds=time.monotonic()-start,anatomical_lumen_calibre_not_claimed_from_centreline_alone=True,atlas_geometry_modified=False)
 (C/'independent-patency.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:report[k]for k in['complete','passed','probes','failures','elapsed_seconds']}))
if __name__=='__main__':main()
