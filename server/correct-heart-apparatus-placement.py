"""Compose the earlier cardiac registration with the documented body resize.

Keep literal source vertex/normal banks. Translate the two anterior AV cusps
into the current body frame, then carry that same correction into the already
installed female similarity. Named posterior correspondences are source-bound.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/args.stage;out.mkdir(exist_ok=False)
male_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';male=json.loads(male_path.read_text());pp={p['id']:p for p in male['parts']};reg=pp['TA98CARDIAC:FJ2418']['registration'];S=np.asarray(reg['scale']);T=np.asarray(reg['translation']);H=male['heightAdjustment']['scale'];t=np.asarray(male['heightAdjustment']['translationMeters']);delta=H*T+t-S*t-T;reader=AtlasReader();proofs=[]
identities={
 'ZA:Inferior leaflet of right atrioventricular valve':('A12.1.02.005','Posterior leaflet of right atrioventricular valve','FJ2433'),
 'ZA:Inferior papillary muscle of right ventricle':('A12.1.02.018','Posterior papillary muscle of right ventricle','FJ2430'),
}
for sex in ['male','female']:
 path=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json';a=json.loads(path.read_text());prior=copy.deepcopy(a);parts={p['id']:p for p in a['parts']};concepts={c['id']:c for c in a['concepts']};marker=a['finiteCardiacIntegration'];changed=[]
 for p in a['parts']:
  source=p.get('cardiacTemplateSourcePartId',p['id'])
  if source in ['TA98HEART:male:tricuspid-anterior-cusp','TA98HEART:male:mitral-anterior-cusp']:
   before=reader.read_part(a,p,frame=('model:local-ta98:assembled-metres'if sex=='male'else'model:local-female-ta98:assembled-metres'),catalogue_path=path);d=delta.copy()
   if sex=='female':
    q=p['geometryInstance'];d=Rotation.from_quat(q['rotation']).apply(d*q['scale'])
   else:p['geometryInstance']=dict(rotation=[0.,0.,0.,1.],scale=1.)
   p['geometryOrigin']=(np.asarray(p['geometryOrigin'])+d).tolist();p['bounds']=(np.asarray(p['bounds'])+d).tolist();p['cardiacRegistrationCorrection']=dict(source_catalogue_sha256=digest_file(path),method='Compose preserved original cardiac registration with the recorded whole-body similarity; transport into recipient similarity once.',translation_m=d.tolist(),source_translation_correction_m=delta.tolist(),vertex_and_normal_words_preserved=True)
   material=p['material'];a['materials'][material].update(geometryOrigin=p['geometryOrigin'],geometryInstance=copy.deepcopy(p['geometryInstance']))
   after=reader.read_part(a,p,frame=before.frame,catalogue_path=path);assert np.array_equal(before.faces,after.faces);assert np.max(abs(after.vertices-before.vertices-d))<3e-16
   proof=audit(after,require_closed=True,max_pairs=1300000,max_tests=14000000);atomic_json(out/f'{sex}-{source.rsplit(":",1)[-1]}-audit.json',proof);assert proof['passed'],proof['failures'];proofs.append(dict(sex=sex,id=p['id'],source_id=source,translation_m=d.tolist(),geometry_sha256=after.geometry_sha,passed=True));changed.append(p['id'])
  if source in identities:
   term,name,source_id=identities[source];p.update(name=name,ta98Term=term);p['cardiacIdentityEvidence']=dict(named_source=source_id,correspondence_report=str(b/'apparatus-correspondence-v121/correspondence.json'),reference='https://www.jacc.org/doi/10.1016/j.jcmg.2018.07.032',qualification='Source-qualified attitudinal inferior/posterior correspondence; complete attachments remain unverified.')
   c=concepts[p['conceptId']];c.update(name=name,ta98Term=term,ta98Kind='exact');c['coverageLimitation']='Named source tissue. Detailed annular/chordal attachments remain unverified.'
   for g in marker['termGroups']:
    if g['partIds']==[p['id']]:g.update(name=name,term=term,kind='exact')
   changed.append(p['id'])
 marker['apparatusPlacementCorrection']=dict(translation_m=delta.tolist(),source_frame='current male body',originating_catalogue_sha256=digest_file(male_path),evidence=str(b/'apparatus-correspondence-v121/correspondence.json'));target=out/sex;target.mkdir();atomic_json(target/'atlas.json',a);atomic_json(target/'changes.json',dict(catalogue_sha256=digest_file(path),changed_parts=changed,geometry_banks_added=0,controls=len(marker['termGroups']),direct_mapped_terms=len({g['term']for g in marker['termGroups']if g['kind']=='exact'}),anatomical_acceptance=False))
atomic_json(out/'construction.json',dict(passed=True,source_registration=reg,body_similarity=male['heightAdjustment'],translation_correction_m=delta.tolist(),parts=proofs,geometry_banks_added=0,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False));print(json.dumps(dict(translation_correction_mm=(delta*1000).tolist(),passed=True,parts=len(proofs),new_geometry_banks=0)),flush=True)
