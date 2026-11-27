"""Before/after render requests with identical nonanatomical framing extent."""
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue'])
    assert digest_file(path)==family['source_catalogue_sha256'];atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']}
    after=json.loads((BASE/'rooted-course-render-v14/queue-receipt.json').read_text());out=BASE/'rooted-course-comparison-v15';out.mkdir();reader=AtlasReader();requests=[]
    for scene in after['requests']:
        system=scene['system'];native=[];bounds=[]
        for i,row in enumerate(scene['palette']):
            m=Mesh.load(row['mesh']);bounds.append(m.bounds)
            if row['name'].startswith('CANDIDATE INITIALIZATION '):
                pid=row['part_id'];before=reader.read_part(atlas,parts[pid],frame=m.frame);file=out/(system+f'-source-{i:02d}.npz');before.save(file)
                native.append(dict(row,mesh=str(file.resolve()),name='NATIVE SOURCE BEFORE '+parts[pid]['name'],source_geometry_sha256=before.geometry_sha));bounds.append(before.bounds)
            else:native.append(row)
        extent=np.array([np.array(bounds)[:,0].min(0),np.array(bounds)[:,1].max(0)])
        vertices=np.array([[extent[x,0],extent[y,1],extent[z,2]]for x in [0,1]for y in [0,1]for z in [0,1]])
        faces=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
        file=out/(system+'-nonanatomical-fixed-framing-box.npz');Mesh(vertices,faces,'model:male:assembled-metres',source=dict(nonanatomical_framing_only=True)).save(file)
        framing=dict(mesh=str(file.resolve()),part_id='REVIEW:nonanatomical-common-framing',name='TRANSPARENT NONANATOMICAL FIXED FRAMING EXTENT',color=[0.,0.,0.,0.])
        for phase,palette in [('before',native),('after',scene['palette'])]:
            spec=dict(model='local-ta98',frame='model:male:assembled-metres',meshes=palette+[framing],resolution=1100,wall_seconds=60,
              views=[dict(name=system+' '+phase+' versus fixed HRA; anterior',direction=[0,0,1]),
                     dict(name=system+' '+phase+' versus fixed HRA; posterior',direction=[0,0,-1]),
                     dict(name=system+' '+phase+' versus fixed HRA; left oblique',direction=[1,.25,.3])])
            (out/(system+'-'+phase+'-spec.json')).write_text(json.dumps(spec,indent=2)+'\n')
            req=prepare(None,spec);requests.append(dict(system=system,phase=phase,palette=palette,common_framing_extent_m=extent.tolist(),**req))
    evidence=dict(schema='pulmonary-source-course-before-after-same-framing-render/v15',
      source_catalogue_sha256=family['source_catalogue_sha256'],candidate_evidence_sha256=digest_file(BASE/'rooted-course-initialization-v13/candidate-evidence.json'),
      requests=requests,common_framing_is_transparent_nonanatomical_geometry=True,
      candidate_rejected_for_installation=True,geometry_changed=False,catalogue_changed=False,
      elapsed_seconds=time.monotonic()-began)
    (out/'queue-receipt.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),requests=[{k:v for k,v in r.items()if k not in ['palette','common_framing_extent_m']}for r in requests],
      receipt_sha256=digest_file(out/'queue-receipt.json'),elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__':main()
