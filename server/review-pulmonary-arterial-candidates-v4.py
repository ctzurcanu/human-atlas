"""Queue literal source branch components and fixed, separately sourced guides."""
import colorsys
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.review import prepare

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began = time.monotonic()
    audit = json.loads((BASE/'evidence.json').read_text())
    contract = json.loads((BASE/'family-contract-v2.json').read_text())
    path = Path(audit['source_catalogue'])
    assert digest_file(path) == contract['source_catalogue_sha256']
    atlas = json.loads(path.read_text())
    parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    target = BASE/'arterial-review-v4'
    target.mkdir()
    requests = []
    source_arrays = np.load(BASE/'literal-components.npz', allow_pickle=False)
    for side, segment, original_guide, installed_guide, territory in [
        ('right', 'Anterior', 'ZA:Anterior segmental bronchus of right lung (BIII)', 'right_anterior_bronchus', 'right-superior:anterior'),
        ('right', 'Posterior', 'ZA:Posterior segmental bronchus of right lung (BII)', 'right_posterior_bronchus', 'right-superior:posterior'),
        ('left', 'Anterior', 'ZA:Anterior segmental bronchus of left lung (BIII)', 'left_anterior_bronchus', 'left-superior:anterior'),
        ('left', 'Posterior', 'ZA:Apicoposterior segmental bronchus of left lung (BI+BII)', 'left_posterior_bronchus', 'left-superior:apicoposterior'),
    ]:
        ident = f'ZA:{segment} segmental artery of {side} lung'
        pi, record = next((i, r)for i,r in enumerate(audit['sources'])if r['id']==ident)
        m = reader.read_part(atlas, parts[ident], frame='model:male:assembled-metres')
        assert m.geometry_sha == record['geometry_sha256']
        rows = []
        palette = []
        for ci, component in enumerate(record['components']):
            ids = source_arrays[f'p{pi}_c{ci}_source_face_ids']
            used, inv = np.unique(m.faces[ids], return_inverse=True)
            sub = Mesh(m.vertices[used], inv.reshape(-1, 3), m.frame,
                       {k: v[used] for k,v in m.attributes.items()},
                       dict(part_id=ident, parent_geometry_sha256=m.geometry_sha,
                            original_face_ids=ids.tolist(), literal_source_subset=True))
            output = target/f'{side}-{segment.lower()}-component-{ci}.npz'
            sub.save(output)
            color = [*colorsys.hsv_to_rgb((ci*.61803398875)%1., .75, .85), 1.]
            name = f'{side} {segment.lower()} arterial source component {ci}; identity unassigned'
            rows.append(dict(mesh=str(output.resolve()), name=name, part_id=ident+f'#source-component-{ci}', color=color))
            palette.append(dict(component=ci, name=name, color=color, source_geometry_sha256=sub.geometry_sha,
                                source_face_ids=ids.tolist()))
        guides = [
            (original_guide, 'ORIGINAL ZA hidden bronchial guide; source reference only', [.93,.88,.69,.55]),
            ('LOCAL:male:pulmonary:bronchus:'+installed_guide, 'FIXED INSTALLED HRA bronchial guide; different source registration', [.30,.75,.92,.80]),
            ('TA98PULMONARY:male:'+territory, 'FIXED INSTALLED HRA territory; provisional anatomical boundary', [.45,.58,.70,.10]),
            ('ZA:'+side.title()+' pulmonary artery', 'FIXED NATIVE main pulmonary artery; source outer surface', [.72,.33,.30,.65]),
        ]
        for gi, (gid, label, color) in enumerate(guides):
            gp = parts[gid]
            gm = reader.read_part(atlas, gp, frame=m.frame)
            output = target/f'{side}-{segment.lower()}-guide-{gi}.npz'
            gm.save(output)
            rows.append(dict(mesh=str(output.resolve()), name=label, part_id=gid, color=color))
            palette.append(dict(part_id=gid, name=label, color=color, geometry_sha256=gm.geometry_sha,
                                default_hidden=gid in atlas['defaultHidden'], suppressed=gp.get('suppressed',False),
                                license=gp.get('license'), provenance=gp.get('provenance')))
        spec = dict(model='local-ta98', frame=m.frame, meshes=rows, resolution=1000, wall_seconds=60,
                    views=[dict(name=f'{side} {segment.lower()} native branch components and separately sourced fixed guides; anterior', direction=[0,0,1]),
                           dict(name=f'{side} {segment.lower()} native branch components; posterior', direction=[0,0,-1]),
                           dict(name=f'{side} {segment.lower()} source takeoff and course; lateral oblique', direction=[1 if side=='left'else -1,.3,.4])])
        specpath = target/f'{side}-{segment.lower()}-spec.json'
        atomic_json(specpath, spec)
        request = prepare(None, spec)
        requests.append(dict(side=side, segment=segment, palette=palette, **request))
    receipt = dict(schema='pulmonary-source-arterial-candidate-render-queue/v4',
                   source_catalogue_sha256=contract['source_catalogue_sha256'],
                   family_contract_sha256=digest_file(BASE/'family-contract-v2.json'),
                   requests=requests, geometry_changed=False, catalogue_changed=False,
                   anatomy_accepted=False, review_pending=True,
                   note='No ascending/descending term assigned; root anatomical/render review must resolve component takeoff/course and current-recipient registration.',
                   elapsed_seconds=time.monotonic()-began)
    atomic_json(target/'queue-receipt.json', receipt)
    print(json.dumps(dict(output=str(target), requests=[{k:v for k,v in r.items()if k!='palette'}for r in requests],
                          receipt_sha256=digest_file(target/'queue-receipt.json'), elapsed_seconds=receipt['elapsed_seconds'])))


if __name__=='__main__':
    main()
