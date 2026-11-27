"""Source-bound before/after S2 arterial construction with fixed recipients."""
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();candidate=json.loads((BASE/'right-S2-terminal-bound-candidate-v22/candidate-evidence.json').read_text())
    family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue']);assert digest_file(path)==candidate['source_catalogue_sha256']
    atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();out=BASE/'right-S2-before-after-review-v23';out.mkdir();before=[];after=[];guides=[];bounds=[]
    for i,row in enumerate(candidate['selected_fragment_meshes']):
        m=Mesh.load(row['mesh']);old=reader.read_part(atlas,parts[row['part_id']],frame=m.frame);ff=np.array(m.source['source_face_ids']);used,inv=np.unique(old.faces[ff],return_inverse=True)
        old=Mesh(old.vertices[used],inv.reshape(-1,3),m.frame,{k:x[used]for k,x in old.attributes.items()},dict(source_part_id=row['part_id'],source_component=row['component'],source_face_ids=ff.tolist()))
        file=out/f'literal-before-fragment-{i}.npz';old.save(file);bounds.extend([old.bounds,m.bounds])
        recurrent='Apical segmental artery'in row['part_id'];color=[.95,.48,.10,1.]if recurrent else[.67,.25,.88,1.]
        name=('A2 recurrent/descending'if recurrent else'A2 ascending')+f'; source component{row["component"]}'
        before.append(dict(mesh=str(file.resolve()),part_id=row['part_id']+f'#component{row["component"]}',name='BEFORE '+name,color=color))
        after.append(dict(mesh=row['mesh'],part_id=row['part_id']+f'#component{row["component"]}',name='AFTER '+name,color=color))
    audit=json.loads((BASE/'evidence.json').read_text());arr=np.load(BASE/'literal-components.npz',allow_pickle=False)
    pid='ZA:Apical segmental artery of right lung';pi=next(i for i,r in enumerate(audit['sources'])if r['id']==pid);m=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres')
    ff=arr[f'p{pi}_c9_source_face_ids'];used,inv=np.unique(m.faces[ff],return_inverse=True);sub=Mesh(m.vertices[used],inv.reshape(-1,3),m.frame,{k:x[used]for k,x in m.attributes.items()},dict(part_id=pid,source_component=9,source_face_ids=ff.tolist()))
    file=out/'fixed-source-upper-trunk-component9.npz';sub.save(file);bounds.append(sub.bounds)
    guides.append(dict(mesh=str(file.resolve()),part_id=pid+'#upper-trunk9',name='FIXED actual superior-trunk parent; source component9',color=[.72,.24,.18,1.]))
    gids=['ZA:Right pulmonary artery']+family['fixed_current_airway_ids']+family['fixed_current_territory_ids']
    for i,pid in enumerate(gids):
        m=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres');file=out/f'fixed-guide-{i}.npz';m.save(file);bounds.append(m.bounds)
        if pid=='ZA:Right pulmonary artery':color=[.72,.24,.18,1.]
        elif pid=='LOCAL:male:pulmonary:bronchus:right_posterior_bronchus':color=[.15,.81,.95,.85]
        elif pid.startswith('LOCAL:male:pulmonary:bronchus:'):color=[.32,.64,.74,.28]
        elif pid=='TA98PULMONARY:male:right-superior:posterior':color=[.75,.85,.77,.14]
        else:color=[.58,.66,.72,.035]
        guides.append(dict(mesh=str(file.resolve()),part_id=pid,name='FIXED CURRENT '+parts[pid]['name'],color=color))
    bb=np.array([np.array(bounds)[:,0].min(0),np.array(bounds)[:,1].max(0)]);v=np.array([[bb[x,0],bb[y,1],bb[z,2]]for x in [0,1]for y in [0,1]for z in [0,1]])
    f=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
    file=out/'transparent-nonanatomical-common-framing.npz';Mesh(v,f,'model:male:assembled-metres',source=dict(nonanatomical_framing_only=True)).save(file)
    guide=dict(mesh=str(file.resolve()),part_id='REVIEW:nonanatomical-common-framing',name='TRANSPARENT NONANATOMICAL SHARED FRAMING EXTENT',color=[0.,0.,0.,0.])
    requests=[]
    for phase,rows in [('before',before),('after',after)]:
        spec=dict(model='local-ta98',frame='model:male:assembled-metres',meshes=rows+guides+[guide],resolution=1100,wall_seconds=60,
          views=[dict(name='Right S2 '+phase+'; anterior',direction=[0,0,1]),dict(name='Right S2 '+phase+'; posterior',direction=[0,0,-1]),dict(name='Right S2 '+phase+'; right lateral oblique',direction=[-1,.2,.35])])
        (out/(phase+'-spec.json')).write_text(json.dumps(spec,indent=2)+'\n');request=prepare(None,spec)
        requests.append(dict(phase=phase,palette=rows+guides,common_framing_bounds_m=bb.tolist(),**request))
    evidence=dict(schema='right-S2-proper-ring-before-after-fixed-recipient-render/v23',source_catalogue_sha256=candidate['source_catalogue_sha256'],candidate_evidence_sha256=digest_file(BASE/'right-S2-terminal-bound-candidate-v22/candidate-evidence.json'),
      requests=requests,render_review_pending=True,finite_wall_and_context_acceptance_pending=True,catalogue_changed=False,elapsed_seconds=time.monotonic()-began)
    (out/'queue-receipt.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),requests=[{k:v for k,v in r.items()if k not in ['palette','common_framing_bounds_m']}for r in requests],
      receipt_sha256=digest_file(out/'queue-receipt.json'),elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__':main()
