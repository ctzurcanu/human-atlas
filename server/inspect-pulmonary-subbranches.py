"""Recover literal source cylinders for four upper-lobe arterial assemblies.

This is inspection evidence, not an assignment of ascending/descending names.
Source topology, source attributes and artery coordinates remain unchanged.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.tube_station_topology import recover_station_rings

SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '6526654a0e412947121518636d1b8004a934093f8a4cf459afcfa55c80a1b63c'
OUTPUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-upper-pulmonary-subbranches-v1'
ASSEMBLIES = [
    ('right-anterior', 'ZA:Anterior segmental artery of right lung', 'A12.2.01.104', ['A12.2.01.105', 'A12.2.01.106']),
    ('right-posterior', 'ZA:Posterior segmental artery of right lung', 'A12.2.01.107', ['A12.2.01.108', 'A12.2.01.109']),
    ('left-anterior', 'ZA:Anterior segmental artery of left lung', 'A12.2.01.205', ['A12.2.01.206', 'A12.2.01.207']),
    ('left-posterior', 'ZA:Posterior segmental artery of left lung', 'A12.2.01.208', ['A12.2.01.209', 'A12.2.01.210']),
]


def main():
    raw = SOURCE.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    assert sha == EXPECTED, 'Source revision changed; explicitly rebase inspection'
    assert not OUTPUT.exists(), 'Keep prior inspection evidence'
    atlas = json.loads(raw)
    parts = {p['id']: p for p in atlas['parts']}
    inventory = json.loads((ROOT / 'reports/ta98-blood-circulation-inventory-installed-v1.json').read_bytes())
    authority = {r['id']: r for r in inventory['terms']}
    OUTPUT.mkdir()
    reader = AtlasReader(ROOT)
    evidence = []
    for slug, part_id, term, children in ASSEMBLIES:
        assert all(authority[c]['parent_id'] == term for c in children)
        part = parts[part_id]
        mesh = reader.read_part(atlas, part, frame='model:local-ta98:assembled-metres', catalogue_path=SOURCE, catalogue_sha=sha)
        mesh.save(OUTPUT / f'{slug}-source.npz')
        components = recover_station_rings(mesh.vertices, mesh.faces)
        rows = []
        for i, component in enumerate(components):
            rings = np.asarray(component['rings'])
            centres = mesh.vertices[rings].mean(axis=1)
            path = OUTPUT / f'{slug}-component-{i}.npz'
            np.savez_compressed(path, positions=mesh.vertices, faces=mesh.faces[component['face_ids']],
                                source_face_ids=component['face_ids'], rings=rings, centres=centres,
                                frame=np.array(mesh.frame), source_geometry_sha256=np.array(mesh.geometry_sha))
            rows.append(dict(component=i, file=str(path.relative_to(ROOT)), stations=len(rings),
                             ring_vertices=rings.shape[1], source_faces=len(component['face_ids']),
                             ends_m=centres[[0,-1]].tolist(),
                             source_centreline_length_m=float(np.linalg.norm(np.diff(centres,axis=0),axis=1).sum())))
        evidence.append(dict(part_id=part_id, authority_parent=term, planned_children=children,
                             source_geometry_sha256=mesh.geometry_sha, source_bank_sha256=mesh.source['bank_sha256'],
                             source_record=part, components=rows, source_attributes_literal=True,
                             branch_identity_assignment='not yet established from source topology alone'))
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == sha
    result=dict(schema='pulmonary-upper-subbranch-source-inspection/v1', source_catalogue_sha256=sha,
                assemblies=evidence, anatomy_accepted=False, source_unchanged=True,
                planned_fixed_context=['pulmonary artery trunks', 'segmental bronchi', 'upper-lobe segments', 'skeleton'],
                primary_anatomical_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC11461042/',
                reference_limit='Branch-pattern evidence; not a registration or patient-specific dimensional source')
    (OUTPUT/'inspection.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUTPUT), assemblies=[dict(id=r['part_id'],components=len(r['components']),
                     stations=[c['stations']for c in r['components']])for r in evidence], source_unchanged=True)))


if __name__ == '__main__':
    main()
