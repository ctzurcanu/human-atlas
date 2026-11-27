"""Independent occupied-tissue comparison, including fully nested overlap."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--stage',required=True);ap.add_argument('--catalogue');ap.add_argument('--render-pose',choices=['part','legacy-material'],default='part');args=ap.parse_args();start=time.monotonic()
path=ROOT/(args.catalogue or f'.local-models/ta98-review/bowel-wall-local-atlas/{args.sex}/atlas.json');a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};m=a['finiteCardiacIntegration'];assert m['sex']==args.sex;reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres'
out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def solid(ids):
 records=[]
 for i in ids:
  p=parts[i]
  if args.render_pose=='legacy-material':
   p=dict(p);legacy=a.get('materials',{}).get(p.get('material'),{}).get('geometryInstance')
   if legacy:p['geometryInstance']=legacy
   else:p.pop('geometryInstance',None)
  records.append(p)
 ms=[reader.read_part(a,p,frame=frame,catalogue_path=path)for p in records];offsets=np.r_[0,np.cumsum([len(m.vertices)for m in ms])]
 v,iv=np.unique(np.concatenate([m.vertices for m in ms]),axis=0,return_inverse=True);f=np.concatenate([iv[m.faces+offsets[i]]for i,m in enumerate(ms)])
 _,ix,inverse,count=np.unique(np.sort(f,axis=1),axis=0,return_index=True,return_inverse=True,return_counts=True);assert count.max()<=2
 parity=np.sum(f[:,[0,0,1]]>f[:,[1,2,2]],axis=1)%2
 orientation=np.bincount(inverse,weights=np.where(parity,-1,1))
 assert np.all(orientation[count==2]==0),'Same-facing duplicate faces cannot cancel'
 f=f[ix[count==1]]
 # Literal shared tissue faces cancel; no tolerance weld or hole closure.
 return make_solid(v,f)[0]
heart_ids=sum([m.get(k,[])for k in ['wallPartIds','valvePartIds','papillaryPartIds','chordalPartIds']],[])
finite_ids=[i for i in heart_ids if parts[i].get('pairedTissue')and parts[i].get('sectionSolid')]
heart=solid(finite_ids);rows=[]
for side in ['left','right']:
 ids=[i for i in a['pulmonarySegmentIntegration']['partIds']if i.startswith(f'TA98PULMONARY:{args.sex}:{side}-')];lung=solid(ids);overlap=heart^lung
 rows.append(dict(side=side,lung_part_ids=ids,intersection_m3=float(overlap.volume()),intersection_ml=float(overlap.volume()*1e6),intersection_triangles=int(overlap.num_tri()),closed_heart_tissue_volume_ml=float(heart.volume()*1e6)))
report=dict(source_catalogue_sha256=digest_file(path),render_pose_source=args.render_pose,heart_part_ids=finite_ids,unrepresented_or_nonfinite_valve_parts=[i for i in heart_ids if i not in finite_ids],roots_checked_separately=True,sides=rows,no_heart_lung_tissue_overlap=all(r['intersection_triangles']==0 for r in rows),method='closed oriented occupied tissues; exact-coordinate opposed internal faces cancelled; Float64 Manifold intersections',same_facing_duplicate_cancellation_rejected=True,geometry_modified=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'report.json',report);print(json.dumps({k:v for k,v in report.items()if k in ['sides','no_heart_lung_tissue_overlap','elapsed_seconds']},indent=2),flush=True)
