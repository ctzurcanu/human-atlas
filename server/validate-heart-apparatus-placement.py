"""Independent replay of AV cusp registration in current male/female frames."""
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--candidate');ap.add_argument('--report',required=True);args=ap.parse_args();reader=AtlasReader();base=ROOT/args.candidate if args.candidate else ROOT/'.local-models/ta98-review/bowel-wall-local-atlas';paths={s:base/s/'atlas.json'for s in ['male','female']};atlases={s:json.loads(p.read_text())for s,p in paths.items()};pp={s:{p['id']:p for p in a['parts']}for s,a in atlases.items()};male=atlases['male'];reg=pp['male']['TA98CARDIAC:FJ2418']['registration'];S=np.array(reg['scale']);T=np.array(reg['translation']);height=male['heightAdjustment'];H=height['scale'];t=np.array(height['translationMeters']);composed=T*H+t-S*t;rows=[];male_meshes={}
for kind,sid in [('tricuspid','FJ2421'),('mitral','FJ2420')]:
 pid=f'TA98HEART:male:{kind}-anterior-cusp';p=pp['male'][pid];source=reader.read_part(male,pp['male'][sid],frame='model:local-ta98:assembled-metres',catalogue_path=paths['male']);actual=reader.read_part(male,p,frame=source.frame,catalogue_path=paths['male']);expected=source.triangles*S+composed;assert expected.shape==actual.triangles.shape
 # The bank was rounded once near its anatomical origin. The current pose
 # must preserve exactly that rounding envelope, not a millimetre tolerance.
 local=expected-np.asarray(p['geometryOrigin']);ulp=.5*np.abs(np.spacing(local.astype('float32')).astype(float));error=np.abs(actual.triangles-expected);bound=ulp+32*np.finfo(float).eps*np.maximum(1,np.abs(expected));assert np.all(error<=bound),(pid,error.max(),bound.max())
 stale_expected=source.triangles*S+T;stale_error=float(np.linalg.norm(stale_expected-actual.triangles,axis=2).max());assert stale_error>.004
 rows.append(dict(sex='male',part=pid,source=sid,maximum_coordinate_error_m=float(error.max()),maximum_analytical_rounding_bound_m=float(bound.max()),stale_pre_height_translation_error_m=stale_error,triangles=len(actual.faces),passed=True));male_meshes[pid]=actual
for pid,m in male_meshes.items():
 p=next(p for p in atlases['female']['parts']if p.get('cardiacTemplateSourcePartId')==pid);actual=reader.read_part(atlases['female'],p,frame='model:local-female-ta98:assembled-metres',catalogue_path=paths['female']);fit=atlases['female']['finiteCardiacIntegration']['transfer'];r=np.array(fit['rotation']);scale=fit['scale'];translation=np.array(fit['translation_m']);expected=m.triangles@r.T*scale+translation;error=np.abs(actual.triangles-expected);bound=48*np.finfo(float).eps*np.maximum(1,np.abs(expected));assert np.all(error<=bound),(pid,error.max())
 rows.append(dict(sex='female',part=p['id'],male_template=pid,maximum_coordinate_error_m=float(error.max()),maximum_roundoff_bound_m=float(bound.max()),triangles=len(actual.faces),passed=True))
for sex,a in atlases.items():
 for term,source in [('A12.1.02.005','FJ2433'),('A12.1.02.018','FJ2430')]:
  g=next(g for g in a['finiteCardiacIntegration']['termGroups']if g['kind']=='exact'and g['term']==term);p=pp[sex][g['partIds'][0]];assert p['cardiacIdentityEvidence']['named_source']==source;assert p['ta98Term']==term
 source=json.loads((ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json').read_text());original={p['id']:p for p in source['parts']};owned=set(g['id']for g in a['parts']if g.get('cardiacRegistrationCorrection')or g.get('cardiacIdentityEvidence'));assert a['chunks']==source['chunks'];assert len(a['parts'])==len(source['parts'])
 assert all(p==original[p['id']]for p in a['parts']if p['id']not in owned),'An unrelated source/anchor changed'
atomic_json(Path(args.report),dict(passed=True,catalogue_sha256={s:digest_file(p)for s,p in paths.items()},parts=rows,part_count_unchanged=True,banks_unchanged=True,unrelated_parts_literal=True));print(json.dumps(rows,indent=2))
