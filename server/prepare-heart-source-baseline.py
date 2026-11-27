"""Source-bound heart inventory and actual native chamber/root evidence."""
from pathlib import Path
import argparse,json,sqlite3,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);args=ap.parse_args();sex=args.sex;b=ROOT/'.local-models/geometry-workbench/heart';out=b/('source-v1-'+sex);out.mkdir(parents=True,exist_ok=False);source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';a=json.loads(source.read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-'+('female-'if sex=='female'else'')+'ta98:assembled-metres';start=time.monotonic();authority=ROOT/'.local-models/ta98-authority.sqlite'
with sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True)as c:
 c.row_factory=sqlite3.Row;terms=[dict(row)for row in c.execute("SELECT * FROM ta98 WHERE id LIKE 'A12.1.%' ORDER BY id")]
metadata=json.loads((ROOT/'app/data/ta98-metadata.json').read_text())['byConceptMatch'];concepts={c['id']:c for c in a['concepts']};hidden=set(a.get('defaultHidden',[]))
def classify(t):
 name=t['name_en'].lower()
 if t['id']=='A12.1.00.001':return'whole-organ'
 if 'valve of foramen'in name:return'tissue-region'
 if any(x in name for x in ['opening','orifice','foramen']):return'opening'
 if any(x in name for x in ['endocardium','myocardium','layer','serosa']):return'wall-layer'
 if any(x in name for x in ['surface','base of','sinus','vestibule','fossa']):return'surface-region'
 if any(x in name for x in ['apex','notch','node','tubercle']):return'point-landmark'if 'node'not in name else'tissue-region'
 if t.get('immaterial')or any(x in name for x in ['border','sulcus','limbus','commissure','vortex']):return'boundary-curve'
 return'tissue-region'
inventory=[]
for t in terms:
 matches=[]
 for p in a['parts']:
  cid=p.get('conceptId');mapped=p.get('ta98Term')or concepts.get(cid,{}).get('ta98Term')or metadata.get(cid,{}).get('term')
  if mapped==t['id']:matches.append(dict(part_id=p['id'],name=p['name'],hidden=p['id']in hidden or bool(p.get('defaultHidden')),suppressed=bool(p.get('suppressed')),kind=concepts.get(cid,{}).get('ta98Kind')or metadata.get(cid,{}).get('kind')))
 inventory.append(dict(authority=t,representation_class=classify(t),source_meshes=matches,missing_distinct_representation=not bool(matches),distinct_representation_verified=False,applicable_to_model=True))
inv=dict(schema='ta98-geometry-inventory/v1',model='local-female-ta98'if sex=='female'else'local-ta98',root_term='A12.1.00.001',authority_path=str(authority),authority_sha256=digest_file(authority),source_catalogue_sha256=digest_file(source),terms=inventory,term_count=len(inventory),includes_pericardial_cross_references=True,anatomical_acceptance=False)
atomic_json(out/'inventory.json',inv)
ids=[p['id']for p in a['parts']if p.get('system')=='cardiac'and not p.get('suppressed')and p['id']not in hidden]
if sex=='male':ids+=['ZA:Ascending aorta','ZA:Pulmonary trunk','ZA:Superior vena cava','ZA:Inferior vena cava (thoracic part)']
else:ids+=['LOCAL:female:aorta:1513']
ids=list(dict.fromkeys(pid for pid in ids if pid in parts));rows=[]
for pid in ids:
 p=parts[pid];m=r.read_part(a,p,frame=frame,catalogue_path=source);proof=audit(m,require_closed=False,max_pairs=300000,max_tests=3000000);atomic_json(out/(pid.replace(':','_').replace('/','_')+'-audit.json'),proof)
 unique,indices=np.unique(m.vertices,axis=0,return_inverse=True)
 try:loops=boundary_loops(indices[m.faces]);mouths=[dict(vertices=q.tolist(),centre=unique[q].mean(0).tolist(),bounds=np.stack([unique[q].min(0),unique[q].max(0)]).tolist())for q in loops]
 except ValueError as exc:mouths=[];proof['boundary_loop_error']=str(exc)
 row=dict(part_id=pid,geometry_sha256=m.geometry_sha,bank_sha256=m.source['bank_sha256'],triangles=len(m.faces),bounds=m.bounds.tolist(),passed=proof['passed'],failures=proof['failures'],topology=proof['topology'],self_crossings=proof['self_intersections']['defect_pairs'],mouths=mouths);rows.append(row);atomic_json(out/'progress.json',dict(parts=rows,elapsed_seconds=time.monotonic()-start),overwrite=True);print(sex,pid,'self',row['self_crossings'],'components',len(proof['topology']['components']),'open',proof['topology']['open_edges'],flush=True)
atomic_json(out/'source-checks.json',dict(model=inv['model'],source_catalogue_sha256=inv['source_catalogue_sha256'],parts=rows,term_count=len(terms),elapsed_seconds=time.monotonic()-start,source_geometry_modified=False,anatomical_acceptance=False))
