"""Compare the two actual right S2 source courses without assigning anatomy."""
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();audit=json.loads((BASE/'evidence.json').read_text());contract=json.loads((BASE/'family-contract-v2.json').read_text())
    source=Path(audit['source_catalogue']);assert digest_file(source)==contract['source_catalogue_sha256']
    atlas=json.loads(source.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();arrays=np.load(BASE/'literal-components.npz',allow_pickle=False)
    out=BASE/'right-posterior-takeoff-review-v6';out.mkdir();rows=[];palette=[]
    for ident,chosen,color,label in [
        ('ZA:Apical segmental artery of right lung',[0,1,2,3],[.95,.43,.12,1.],'Upper-trunk S2-course hypothesis; original source object calls this apical'),
        ('ZA:Posterior segmental artery of right lung',[0,1,2],[.59,.20,.80,1.],'Interlobar S2-course hypothesis; named posterior source'),
        ('ZA:Apical segmental artery of right lung',[9],[.88,.60,.44,.55],'Literal common upper/apical trunk; proposed origin context only'),
    ]:
        pi,rec=next((i,r)for i,r in enumerate(audit['sources'])if r['id']==ident);m=reader.read_part(atlas,parts[ident],frame='model:male:assembled-metres')
        for ci in chosen:
            ids=arrays[f'p{pi}_c{ci}_source_face_ids'];used,inverse=np.unique(m.faces[ids],return_inverse=True)
            sub=Mesh(m.vertices[used],inverse.reshape(-1,3),m.frame,{k:v[used]for k,v in m.attributes.items()},dict(part_id=ident,source_face_ids=ids.tolist(),source_geometry_sha256=m.geometry_sha))
            file=out/f'p{pi}-c{ci}.npz';sub.save(file);name=label+f'; component {ci}'
            rows.append(dict(mesh=str(file.resolve()),name=name,part_id=ident+f'#component-{ci}',color=color));palette.append(dict(name=name,part_id=ident,component=ci,color=color,source_face_ids=ids.tolist(),mesh_sha256=digest_file(file)))
    for gi,(ident,color,label)in enumerate([
        ('ZA:Posterior segmental bronchus of right lung (BII)',[.95,.88,.65,.75],'Original ZA BII source guide; fixed reference'),
        ('ZA:Apical segmental bronchus of right lung (BI)',[.75,.85,.66,.35],'Original ZA BI source guide; fixed reference'),
        ('LOCAL:male:pulmonary:bronchus:right_posterior_bronchus',[.12,.66,.90,.85],'Current installed HRA BII; fixed different-source recipient'),
        ('LOCAL:male:pulmonary:bronchus:right_apical_bronchus',[.20,.80,.73,.30],'Current installed HRA BI; fixed different-source context'),
        ('ZA:Right pulmonary artery',[.60,.38,.30,.65],'Native main RPA/interlobar root; fixed root context'),
        ('ZA:Superior lobar artery of right lung',[.80,.28,.22,.80],'Native superior-lobar origin context'),
        ('TA98PULMONARY:male:right-superior:posterior',[.38,.50,.65,.08],'Current installed provisional S2 territory; fixed'),
    ]):
        p=parts[ident];m=reader.read_part(atlas,p,frame='model:male:assembled-metres');file=out/f'guide-{gi}.npz';m.save(file)
        rows.append(dict(mesh=str(file.resolve()),name=label,part_id=ident,color=color));palette.append(dict(name=label,part_id=ident,color=color,mesh_sha256=digest_file(file),license=p.get('license'),source_geometry_sha256=m.geometry_sha))
    spec=dict(model='local-ta98',frame='model:male:assembled-metres',meshes=rows,resolution=1000,wall_seconds=60,
              views=[dict(name='Right S2 two-source-takeoff hypotheses; posterior',direction=[0,0,-1]),
                     dict(name='Right S2 upper versus interlobar takeoff; lateral',direction=[-1,.10,0]),
                     dict(name='Right S2 source courses; posterolateral oblique',direction=[-1,.35,-.8])])
    atomic_json(out/'spec.json',spec);request=prepare(None,spec)
    evidence=dict(schema='right-posterior-pulmonary-source-hypothesis-render/v6',source_catalogue_sha256=contract['source_catalogue_sha256'],
                  request=request,palette=palette,anatomical_assignments_unaccepted=True,geometry_changed=False,catalogue_changed=False,
                  primary_reference_urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC4454828/','https://pmc.ncbi.nlm.nih.gov/articles/PMC11461042/'],
                  note='Takeoff topology and original BII-course witnesses motivate two hypotheses. No closest-plane classification or current HRA fit is established.',elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'queue-receipt.json',evidence);print(json.dumps(dict(request=request,output=str(out),receipt_sha256=digest_file(out/'queue-receipt.json'),elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__':main()
