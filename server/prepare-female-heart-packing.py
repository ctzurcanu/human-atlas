"""A bounded common heart-instance correction against the native sternum.

The independently retained aortic arch and all fixed context remain literal.
This prepares a candidate, never installs or alters geometry banks.
"""
from pathlib import Path
import argparse, copy, json, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.algorithms.projected_clearance import directional_clearance
from geometry_tools.contact import pair_contacts
from importlib import util

ap=argparse.ArgumentParser()
ap.add_argument('--stage',required=True)
ap.add_argument('--relative-scale',type=float,default=.9)
ap.add_argument('--sternal-clearance-m',type=float,default=.002)
args=ap.parse_args()
assert .8<=args.relative_scale<=1.
assert .001<=args.sternal_clearance_m<=.005
start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart'
out=h/args.stage;out.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
original=json.loads(path.read_text());a=copy.deepcopy(original)
spec=util.spec_from_file_location('cardiac_forest',ROOT/'scripts/build-spatial-dependency-forest.py')
module=util.module_from_spec(spec);spec.loader.exec_module(module)
forest=module.build_forest('female',a)
groups=[g for root in forest['forest']for g in root['children']if g['id']=='assembly:heart']
assert len(groups)==1
moving=groups[0]['activePartIds']
assert moving and 'LOCAL:female:aortic_arch:1513'not in moving
forest['sourceSha256']=digest_file(path)
atomic_json(out/'source-forest.json',forest)
parts={p['id']:p for p in a['parts']};reader=AtlasReader()
frame='model:local-female-ta98:assembled-metres'
def read(pid):return reader.read_part(original,next(p for p in original['parts']if p['id']==pid),frame=frame,catalogue_path=path)
heart=read('TA98HEART:female:epicardium');sternum=read('LOCAL:female:sternum:342')
annuli=[read(f'TA98HEART:female:fibrous:{side}-annulus')for side in ['right','left']]
pivot=np.mean([m.vertices.mean(0)for m in annuli],axis=0)
scale=args.relative_scale
scaled=Mesh(pivot+scale*(heart.vertices-pivot),heart.faces,frame)
before=directional_clearance(heart,sternum)
after_scale=directional_clearance(scaled,sternum)
translation=np.array([0.,0.,-max(0.,after_scale['maximum_forward_excess_m']+args.sternal_clearance_m)])
fitted=Mesh(scaled.vertices+translation,scaled.faces,frame)
after=directional_clearance(fitted,sternum)
assert after['maximum_forward_excess_m']<=-args.sternal_clearance_m+1e-10
contacts=[]
for name,m in [('before',heart),('candidate',fitted)]:
 r,proof=pair_contacts(m,sternum,0.,max_pairs=250000,max_tests=2000000)
 assert proof['complete_within_pair_budget']
 atomic_json(out/(name+'-sternum.json'),proof)
 contacts.append(dict(name=name,**proof))
assert contacts[-1]['crossing_face_pairs']==0 and contacts[-1]['coplanar_overlap_face_pairs']==0
for pid in moving:
 p=parts[pid]
 assert p.get('geometryInstance')and p.get('geometryOrigin')is not None,pid
 p['geometryOrigin']=(pivot+scale*(np.array(p['geometryOrigin'])-pivot)+translation).tolist()
 p['geometryInstance']['scale']*=scale
 actual=reader.read_part(a,p,frame=frame,catalogue_path=path)
 source=read(pid)
 assert actual.triangles.shape==source.triangles.shape
 error=float(np.abs(actual.triangles-(pivot+scale*(source.triangles-pivot)+translation)).max())
 assert error<2e-15,(pid,error)
 p['bounds']=actual.bounds.tolist()
 p['cardiacThoracicFit']=dict(schema='cardiac-thoracic-packing/v1',relativeScale=scale,
  pivotM=pivot.tolist(),translationM=translation.tolist(),sourceCatalogueSha256=digest_file(path))
for g in a.get('cardiacFeatureIntegration',{}).get('termGroups',[]):
 g['anchor']['point']=(pivot+scale*(np.array(g['anchor']['point'])-pivot)+translation).tolist()
marker=a['finiteCardiacIntegration']
marker['geometryInstance']=copy.deepcopy(parts[marker['wallPartIds'][0]]['geometryInstance'])
marker['thoracicFit']=dict(schema='cardiac-thoracic-packing/v1',sourceCatalogueSha256=digest_file(path),
 relativeScale=scale,scaleComparedWithMale=marker['geometryInstance']['scale'],pivotM=pivot.tolist(),translationM=translation.tolist(),
 movingPartIds=moving,fixedPartIds=[pid for pid in parts if pid not in moving],sternalProjection=after,
 independentRecipientJoinsVerified=False,clinicalSizeEstablished=False)
for pid,p in parts.items():
 if pid not in moving:assert p==next(q for q in original['parts']if q['id']==pid),pid
atomic_json(out/'atlas.json',a)
fitted.save(out/'epicardium.npz');sternum.save(out/'fixed-sternum.npz')
proof=dict(prepared=True,source_catalogue_sha256=digest_file(path),moving_parts=moving,
 relative_scale=scale,scale_compared_with_male=marker['geometryInstance']['scale'],pivot_m=pivot.tolist(),translation_m=translation.tolist(),
 before=before,after_scale=after_scale,after=after,sternal_contacts=contacts,
 banks_modified=False,all_nonmembers_literal=True,internal_coordinates_and_interfaces_preserved=True,
 wider_packing_review_pending=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'packing.json',proof)
print(json.dumps({k:proof[k]for k in ['moving_parts','relative_scale','scale_compared_with_male','translation_m','elapsed_seconds']}),flush=True)
