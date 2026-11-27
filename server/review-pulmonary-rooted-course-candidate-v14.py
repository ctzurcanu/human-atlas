"""Queue two literal whole-family initialization views; parent owns Blender UI."""
import json
import time
from pathlib import Path

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.review import prepare

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();candidate=json.loads((BASE/'rooted-course-initialization-v13/candidate-evidence.json').read_text())
    family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue'])
    assert digest_file(path)==candidate['source_catalogue_sha256'];atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']}
    out=BASE/'rooted-course-render-v14';out.mkdir();reader=AtlasReader();guides=[]
    for i,pid in enumerate(family['fixed_current_airway_ids']+family['fixed_current_territory_ids']):
        mesh=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');dest=out/f'fixed-guide-{i:02d}.npz';mesh.save(dest)
        airway=pid in family['fixed_current_airway_ids'];color=[.27,.72,.86,.65]if airway else[.62,.67,.72,.055]
        guides.append(dict(mesh=str(dest.resolve()),part_id=pid,name='FIXED CURRENT HRA '+parts[pid]['name'],color=color,
                           source_geometry_sha256=mesh.geometry_sha))
    requests=[]
    for system in ['arterial','venous']:
        rows=[]
        for row in candidate['candidate_meshes']:
            pid=row['id'];vein='vein'in pid.lower();show=(system=='venous')==vein
            if not show:continue
            color=[.25,.43,.88,1.]if vein else[.90,.30,.23,1.]
            if pid in ['ZA:Pulmonary trunk','ZA:Bifurcation of pulmonary trunk']:color=[.82,.56,.18,1.]
            if vein and pid in [f'ZA:{side} {level} pulmonary vein'for side in ['Left','Right']for level in ['inferior','superior']]:color=[.38,.28,.68,1.]
            rows.append(dict(mesh=row['mesh'],part_id=pid,name='CANDIDATE INITIALIZATION '+parts[pid]['name'],color=color,
                             source_geometry_sha256=row['geometry_sha256']))
        spec=dict(model='local-ta98',frame='model:male:assembled-metres',meshes=rows+guides,resolution=1100,wall_seconds=60,
          views=[dict(name=system+' source-course INITIALIZATION versus fixed current HRA; anterior',direction=[0,0,1]),
                 dict(name=system+' source-course INITIALIZATION versus fixed current HRA; posterior',direction=[0,0,-1]),
                 dict(name=system+' source-course INITIALIZATION; left oblique',direction=[1,.25,.3])])
        (out/(system+'-spec.json')).write_text(json.dumps(spec,indent=2)+'\n')
        request=prepare(None,spec)
        requests.append(dict(system=system,palette=rows+guides,**request))
    evidence=dict(schema='pulmonary-rooted-source-course-render-queue/v14',
      candidate_evidence_sha256=digest_file(BASE/'rooted-course-initialization-v13/candidate-evidence.json'),
      source_catalogue_sha256=candidate['source_catalogue_sha256'],requests=requests,
      render_review_pending=True,candidate_is_not_anatomical_acceptance=True,
      no_wall_or_patency_claim=True,geometry_changed=False,catalogue_changed=False,
      elapsed_seconds=time.monotonic()-began)
    (out/'queue-receipt.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),requests=[{k:v for k,v in r.items()if k!='palette'}for r in requests],
      receipt_sha256=digest_file(out/'queue-receipt.json'),elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__':main()
