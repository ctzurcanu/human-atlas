"""Before/after diagnostic of rejected v28, preserving all current identities."""
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT=BASE/'coherent-pulmonary-root-before-after-v30'
FRAME='model:male:assembled-metres'


def main():
    began=time.monotonic();OUT.mkdir()
    candidate=json.loads((BASE/'coherent-pulmonary-root-translation-v28/candidate-evidence.json').read_text())
    disposition=json.loads((BASE/'coherent-pulmonary-root-cheap-context-v29/aggregate-disposition.json').read_text());assert disposition['installation_rejected']
    family=json.loads((BASE/'family-contract-v2.json').read_text());assert digest_file(Path(family['source_catalogue']))==candidate['source_catalogue_sha256']
    atlas=json.loads(Path(family['source_catalogue']).read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    moving={r['id']for r in candidate['literal_source_replay']if r['moving']};before=[];after=[];context=[];bounds=[]
    for i,row in enumerate(candidate['candidate_meshes']):
        pid=row['id']
        if pid not in moving and pid not in ['ZA:Pulmonary trunk','ZA:Bifurcation of pulmonary trunk']:continue
        source=reader.read_part(atlas,parts[pid],frame=FRAME);path=OUT/('literal-current-part-%02d.npz'%i)
        Mesh(source.vertices,source.faces,FRAME,source.attributes,dict(source.source,part_id=pid,current_source_part_id=pid)).save(path)
        m=Mesh.load(row['mesh']);bounds.extend([source.bounds,m.bounds])
        is_moving=pid in moving;color=[.90,.35,.17,1.]if is_moving else[.64,.18,.16,1.]
        entry=dict(mesh=str(path.resolve()),part_id=pid,name=('SOURCE CURRENT 'if is_moving else'FIXED CURRENT ')+parts[pid]['name'],color=color)
        if is_moving:
            before.append(entry);after.append(dict(entry,mesh=row['mesh'],name='RECORDED-T INITIALIZATION '+parts[pid]['name']))
        else:context.append(entry)
    for i,pid in enumerate(family['fixed_current_airway_ids']+family['fixed_current_territory_ids']):
        source=reader.read_part(atlas,parts[pid],frame=FRAME);path=OUT/('literal-fixed-HRA-%02d.npz'%i)
        Mesh(source.vertices,source.faces,FRAME,source.attributes,dict(source.source,part_id=pid,current_source_part_id=pid)).save(path);bounds.append(source.bounds)
        color=[.12,.76,.91,.65]if pid in family['fixed_current_airway_ids']else[.63,.69,.74,.04]
        context.append(dict(mesh=str(path.resolve()),part_id=pid,name='FIXED CURRENT '+parts[pid]['name'],color=color))
    bb=np.array([np.asarray(bounds)[:,0].min(0),np.asarray(bounds)[:,1].max(0)])
    v=np.array([[bb[x,0],bb[y,1],bb[z,2]]for x in [0,1]for y in [0,1]for z in [0,1]])
    ff=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
    marker=OUT/'transparent-nonanatomical-shared-framing.npz';Mesh(v,ff,FRAME,source=dict(nonanatomical_framing_only=True)).save(marker)
    frame=dict(mesh=str(marker.resolve()),part_id='REVIEW:nonanatomical-shared-framing',name='NONANATOMICAL TRANSPARENT FRAME ONLY',color=[0,0,0,0],nonanatomical_framing_only=True)
    requests=[]
    for phase,rows in [('before',before),('after',after)]:
        spec=dict(model='local-ta98',frame=FRAME,meshes=rows+context+[frame],resolution=1200,wall_seconds=60,
                  views=[dict(name='Coherent PA root '+phase+'; anterior',direction=[0,0,1]),
                         dict(name='Coherent PA root '+phase+'; posterior',direction=[0,0,-1]),
                         dict(name='Coherent PA root '+phase+'; right oblique',direction=[-1,.2,.4])])
        specpath=OUT/(phase+'-spec.json');specpath.write_text(json.dumps(spec,indent=2)+'\n')
        request=prepare(None,spec)
        task=json.loads((Path(request['folder'])/'request.json').read_text())
        assert all(row.get('part_id')for row in task['meshes'])
        assert sum(row.get('nonanatomical_framing_only',False)for row in task['meshes'])==1
        requests.append(dict(phase=phase,spec=str(specpath.resolve()),spec_sha256=digest_file(specpath),**request))
    result=dict(schema='rejected-coherent-PA-root-translation-before-after-review/v30',
                candidate_evidence_sha256=digest_file(BASE/'coherent-pulmonary-root-translation-v28/candidate-evidence.json'),
                rejection_evidence_sha256=digest_file(BASE/'coherent-pulmonary-root-cheap-context-v29/aggregate-disposition.json'),
                source_catalogue_sha256=candidate['source_catalogue_sha256'],common_frame_bounds_m=bb.tolist(),requests=requests,
                all_current_part_identities_preserved_in_record=True,frame_only_marker_explicit=True,
                diagnostic_not_acceptance=True,installation_rejected=True,elapsed_seconds=time.monotonic()-began)
    (OUT/'queue-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),sha256=digest_file(OUT/'queue-receipt.json'),requests=requests,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
