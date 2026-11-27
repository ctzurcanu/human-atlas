"""One 60-second externally supervised final construction/export stage."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.review import prepare
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1';OUT=WORK/'artery-exterior-corridor-v37'
def run(name):
 spec=importlib.util.spec_from_file_location(name.replace('-','_'),ROOT/'server'/name);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.main();return mod

def main():
 began=time.monotonic();run('construct-cystic-artery-exterior-corridor-v37.py');g=json.loads((OUT/'gates.json').read_text());
 if not g.get('allRequiredGeometryGatesPassed'):return
 if time.monotonic()-began>32:
  atomic_json(OUT/'final-stage-export-held.json',dict(reason='Construction used the reserved export budget; no emitted candidate or review accepted',elapsedSeconds=time.monotonic()-began,installed=False));return
 run('partition-cystic-artery-local-stage3.py');sys.argv=['pack','artery'];run('pack-cystic-artery-local-stage3.py');stage=OUT/'bank/artery';records=json.loads((stage/'records.json').read_text());manifest=json.loads((stage/'manifest.json').read_text());assert manifest['emittedGeometryChecksPassed'];parent=Mesh.load(stage/'parent-emitted.npz');branch=Mesh.load(stage/'cystic-emitted.npz');path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(path)=='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();current_parent=reader.read_part(a,parts['TA98HEPATIC:male:artery'],frame=parent.frame,catalogue_path=path);parent_literal=parent.geometry_sha==current_parent.geometry_sha;parent_words_literal=False
 if np.array_equal(parent.triangles,current_parent.triangles):
  words=lambda m:m.attributes['normals'][m.faces].astype(np.int64)*m.attributes['normalSigns'][m.faces,None];parent_words_literal=np.array_equal(words(parent),words(current_parent))
 es=importlib.util.spec_from_file_location('native_envelope',ROOT/'server/cystic-biliary-actual-envelope.py');em=importlib.util.module_from_spec(es);es.loader.exec_module(em);outer,op=em.vessel_envelope(branch,OUT/'owned/artery/partition.npz',OUT,label='candidate-artery',bank_sha256=manifest['bankSha256']);duct_outer=Mesh.load(OUT/'actual-biliary-outer-computational-envelope.npz');origin=np.array(parts['TA98CYSTICBLOOD:male:artery']['geometryOrigin']);bs,_=make_solid((outer.vertices-origin)*1000,outer.faces);ds,_=make_solid((duct_outer.vertices-origin)*1000,duct_outer.faces);occupation=float((bs^ds).volume());checks=[]
 contexts=[(pid,reader.read_part(a,parts[pid],frame=branch.frame,catalogue_path=path))for pid in ['TA98CYSTICBLOOD:male:vein','TA98HEPATIC:male:portal','TA98HEPATIC:male:hepatic_vein','TA98HEPATIC:male:biliary','FJ3080','TA98PART:FJ3079:A05.8.01.061']]+[('GB-innermost',Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz'))]
 for pid,m in contexts:
  _,c=pair_contacts(branch,m,0.,max_pairs=500000,max_tests=5000000);checks.append(dict(id=pid,geometrySha256=m.geometry_sha,**c))
 passed=parent_literal and parent_words_literal and occupation==0 and all(c['complete_within_pair_budget']and not c['crossing_face_pairs']and not c['coplanar_overlap_face_pairs']for c in checks);proof=dict(schema='local-artery-stage3-complete-emitted/v1',sourceCatalogueSha256=digest_file(path),bankSha256=digest_file(stage/'cells-and-regions.bin'),actualParentGeometryLiteral=parent_literal,actualParentCornerNormalWordsLiteral=parent_words_literal,completeEmittedOuterEnvelopeDuctOccupation_mm3=occupation,outerEnvelope=op,contextChecks=checks,wholeEmittedGatesPassed=passed,installed=False,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'emitted-final-proof.json',proof)
 if passed:
  vein=reader.read_part(a,parts['TA98CYSTICBLOOD:male:vein'],frame=branch.frame,catalogue_path=path);vein.save(OUT/'fixed-installed-vein.npz');duct=reader.read_part(a,parts['TA98HEPATIC:male:biliary'],frame=branch.frame,catalogue_path=path);duct.save(OUT/'fixed-installed-hepatic-biliary.npz');entries=[dict(mesh=str(stage/'cystic-emitted.npz'),name='Emitted corrected cystic artery',color=[.86,.13,.10,1]),dict(mesh=str(OUT/'fixed-installed-vein.npz'),name='Fixed installed cystic vein',color=[.15,.27,.83,1]),dict(mesh=str(OUT/'fixed-installed-hepatic-biliary.npz'),name='Fixed actual hepatic-biliary wall',color=[.72,.47,.38,.45])];review=prepare(None,dict(frame=branch.frame,meshes=entries,views=[dict(name='actual-local-arterial-corridor-front',direction=[0,0,1]),dict(name='actual-local-arterial-corridor-oblique',direction=[1,.2,1])],resolution=1000,wall_seconds=40));atomic_json(OUT/'emitted-review.json',review)
 print(json.dumps(dict(completeEmittedPass=passed,parentLiteral=parent_literal,parentNormalsLiteral=parent_words_literal,occupiedDuctEnvelope_mm3=occupation,allContextCrossings={c['id']:c['crossing_face_pairs']for c in checks},elapsedSeconds=time.monotonic()-began)))
if __name__=='__main__':main()
