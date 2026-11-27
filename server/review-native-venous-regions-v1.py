"""Queue exact emitted native source regions, without anatomical acceptance."""
from pathlib import Path
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-native-venous-regions-v1'
OUT=BASE/'literal-forearm-pack-v1';CANDIDATE=OUT/'proposed-atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'


def main():
    assert digest_file(CANDIDATE)==EXPECTED
    atlas=json.loads(CANDIDATE.read_bytes());plan=json.loads((BASE/'source-plan.json').read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    assert json.loads((OUT/'controls-proof.json').read_bytes())['passed']
    target=OUT/'native-review';assert not target.exists();target.mkdir();requests=[]
    for side in ['left','right']:
        rows=[];groups=[g for g in atlas['nativeVenousRegionIntegration']['groups'] if g['side']==side]
        for group in groups:
            for role in ['remainderPartId','regionPartId']:
                p=parts[group[role]];m=reader.read_part(atlas,p,frame='model:local-ta98:assembled-metres',catalogue_path=CANDIDATE,catalogue_sha=EXPECTED)
                path=target/(hashlib.sha256(p['id'].encode()).hexdigest()[:16]+'.npz');m.save(path)
                region=role=='regionPartId';cephalic=group['regionTerm']=='A12.3.08.021'
                color=([.12,.55,.9,1] if cephalic else [.50,.25,.85,1]) if region else [.32,.42,.56,1]
                rows.append(dict(mesh=str(path),name=(p['name'] if region else group['wholeName']+' · remaining literal source course'),part_id=p['id'],color=color))
        radius='ZA:Radius.'+side[0]+':Bone'
        rows.append(dict(mesh=str(ROOT/plan['sources'][radius]['mesh']),name='Literal fixed native radius · context',part_id=radius,color=[.83,.79,.69,.45]))
        point=next(r for r in plan['regions'] if r['term']=='A12.3.08.021' and r['side']==side)['anatomical_landmarks']
        center=[(a+b)/2 for a,b in zip(point['proximal_m'],point['distal_m'])]
        spec=dict(model='local-ta98',frame='model:local-ta98:assembled-metres',meshes=rows,resolution=800,wall_seconds=60,
            views=[dict(name=side+' actual emitted source partitions · anterior',direction=[0,0,1]),
                   dict(name=side+' diagnostic longitudinal uncapped source surfaces',direction=[1 if side=='left' else -1,0,0],section=dict(normal=[1,0,0],point=center,clear_outer=True))])
        atomic_json(target/(side+'-spec.json'),spec);requests.append(dict(side=side,**prepare(None,spec)))
    result=dict(candidate_sha256=EXPECTED,requests=requests,installed=False,anatomical_acceptance=False,view_review_pending=True,
                limits=['Source outer surfaces only; diagnostic cuts intentionally uncapped.','Named boundaries provisional native radius landmarks.','No repaired blood lumen, valve or upstream/downstream communication is claimed.'])
    atomic_json(target/'queue-receipt.json',result);print(json.dumps(result))


if __name__=='__main__':main()
