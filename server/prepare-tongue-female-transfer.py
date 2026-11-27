"""Read-only native anatomy preflight for transferring the completed tongue.

This exports immutable target witnesses, not a resized tongue or an accepted
fit. Mandibular dimension ratios are explicit initialization hypotheses; they
cannot establish attachment, bone clearance or a correct whole-organ contour.
"""
from pathlib import Path
import argparse
import json
import sqlite3
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology

BASE = ROOT/'.local-models/geometry-workbench/tongue-completion'
FEMALE = ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
MALE = ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'


def main(output):
    began=time.monotonic()
    output=output.resolve()
    if not output.is_relative_to(BASE) or output.exists():
        raise ValueError('Use a new local tongue-transfer preflight directory')
    source_audit=BASE/'source-v1/source-audit.json'
    inventory=ROOT/'reports/ta98-tongue-inventory.json'
    source=json.loads(source_audit.read_text())
    expected={name:row['sha256'] for name,row in source['sourceCatalogues'].items()}
    if digest_file(FEMALE)!=expected['female'] or digest_file(MALE)!=expected['male']:
        raise ValueError('Rebind source/target inventory explicitly after a catalogue change')
    female=json.loads(FEMALE.read_text())
    concepts={row['id']:row for row in female['concepts']}
    parts={row['id']:row for row in female['parts']}
    requested=[
        ('LOCAL:female:tongue:1089','native replacement target; open source, not a solid-volume template'),
        ('LOCAL:female:genioglossus:239','native extrinsic course and mandibular interface context'),
        ('LOCAL:female:l_hyoglossus:242','native left extrinsic/hyoid interface context'),
        ('LOCAL:female:r_hyoglossus:255','native right extrinsic/hyoid interface context'),
        ('LOCAL:female:l_styloglossus:243','native left styloid-course guide'),
        ('LOCAL:female:r_styloglossus:256','native right styloid-course guide'),
        ('LOCAL:female:mandible:960','immutable fixed mandibular surface'),
        ('LOCAL:female:hyoid:1002','immutable fixed whole-hyoid surface'),
        ('LOCAL:female:skull_cranium:999','immutable fixed cranial/styloid context fragment'),
        ('LOCAL:female:skull_cranium:1000','immutable fixed cranial/styloid context fragment'),
        ('LOCAL:female:skull_cranium:1001','immutable fixed cranial/styloid context fragment'),
        ('LOCAL:female:mylohyoid:238','fixed oral-floor adjoining tissue'),
        ('LOCAL:female:l_geniohyoid:241','fixed left oral-floor adjoining tissue'),
        ('LOCAL:female:r_geniohyoid:254','fixed right oral-floor adjoining tissue'),
        ('LOCAL:female:uvula:1088','fixed palatal endpoint guide; not a substitute for soft-palate geometry'),
        ('LOCAL:female:pharynx:1087','fixed posterior adjoining pharyngeal tissue'),
        ('LOCAL:female:nasopharynx:1047','fixed palatal/oral roof neighborhood context')]
    if any(part_id not in parts for part_id,_ in requested):
        raise ValueError('An explicitly declared native target part is missing')
    output.mkdir(parents=True)
    reader=AtlasReader()
    rows=[]
    for index,(part_id,role) in enumerate(requested):
        if time.monotonic()-began>30:
            raise TimeoutError('Female transfer preflight exceeded its 30 second internal deadline')
        part=parts[part_id]
        mesh=reader.read_part(female,part,frame='model:local-female-ta98:assembled-metres',
                             catalogue_path=FEMALE,catalogue_sha=expected['female'])
        path=output/f'target-{index:02d}.npz'
        mesh.save(path)
        info=topology(mesh)
        rows.append(dict(part_id=part_id,role=role,name=part['name'],
            ta98_term=concepts.get(part['conceptId'],{}).get('ta98Term'),
            mesh=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,
            source_bank_sha256=mesh.source['bank_sha256'],source_bank_path=mesh.source['bank_path'],
            source_part=part,bounds_m=mesh.bounds.tolist(),topology=info,
            fixed=part['system']=='skeletal' or part_id not in [row[0]for row in requested[:6]],
            solid_occupancy_accepted=False))
    male_row=next(row for row in source['contextParts'] if row['partId']=='ZA:Mandible:Bone')
    male_jaw=Mesh.load(male_row['meshSavePath'])
    female_jaw=Mesh.load(next(row['mesh']for row in rows if row['part_id']=='LOCAL:female:mandible:960'))
    ratios=np.diff(female_jaw.bounds,axis=0)[0]/np.diff(male_jaw.bounds,axis=0)[0]
    center_shift=female_jaw.bounds.mean(0)-male_jaw.bounds.mean(0)
    authority=ROOT/'.local-models/ta98-authority.sqlite'
    with sqlite3.connect(authority) as connection:
        palate_terms=[row[0]for row in connection.execute(
            "SELECT id FROM ta98 WHERE name_en='soft palate' AND name_la='palatum molle'")]
    if set(palate_terms)!={'A05.1.01.104','A05.2.01.003'}:
        raise ValueError('Resolve actual authoritative soft-palate entries before transfer')
    soft_palate=[part['id']for part in female['parts']
                 if concepts.get(part['conceptId'],{}).get('ta98Term') in palate_terms]
    report=dict(schema='tongue-female-transfer-native-preflight/v1',
        source_catalogues={name:dict(path=str(path),sha256=digest_file(path))
                           for name,path in [('male',MALE),('female',FEMALE)]},
        inventory=dict(path=str(inventory),sha256=digest_file(inventory)),
        authority=dict(path=str(authority),sha256=digest_file(authority),soft_palate_terms=palate_terms),
        source_audit=dict(path=str(source_audit),sha256=digest_file(source_audit)),
        target_frame='model:local-female-ta98:assembled-metres',target_parts=rows,
        dimension_hypotheses=dict(mandible_extent_ratios_xyz=ratios.tolist(),
            uniform_width_ratio=float(ratios[0]),source_to_target_box_center_shift_m=center_shift.tolist(),
            selected_fit=False,anatomical_acceptance=False,
            caution='These ratios initialize registration only; do not shrink height from the open female tongue bounding box.'),
        dependency_forest=[
            dict(assembly='Tongue shared tissue bank',members=['mucosa','aponeurosis','superior-longitudinal'],
                 transform_together=True,literal_reciprocal_charts_required=True),
            dict(assembly='Tongue intrinsic territories',members=['inferior-longitudinal-left','inferior-longitudinal-right','transverse-vertical-shared'],
                 transform_together=True,shared_reference_controls_preserved=True),
            dict(assembly='Bilateral extrinsic courses',members=['genioglossus-left','genioglossus-right','hyoglossus-left','hyoglossus-right','styloglossus-left','styloglossus-right','palatoglossus-left','palatoglossus-right'],
                 initial_transform_together=True,subsequent_fit='explicit source/target attachment witnesses; muscle origins and tissue entries move coherently')],
        adjoining_interfaces=[
            dict(name='mandibular genioglossus origin',host='LOCAL:female:mandible:960',state='actual surface witness required'),
            dict(name='body and horns of hyoid',host='LOCAL:female:hyoid:1002',state='actual finite footprint required; no whole-bone motion'),
            dict(name='styloid origins',host='native cranial fragments',state='actual styloid surface witness required; no box-only origin'),
            dict(name='palatal attachment and faucial mucosa',host='native palate/uvular neighborhood',
                 directly_mapped_soft_palate_parts=soft_palate,
                 state='no directly mapped native soft-palate part found; establish the actual palatal surface before fitting palatoglossus'if not soft_palate else'actual palatal witness required'),
            dict(name='oral floor and posterior pharynx',state='fixed context requires full fitted assembly checks')],
        source_anatomy_modified=False,skeleton_modified=False,geometry_fit_created=False,
        installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    for name,path in [('male',MALE),('female',FEMALE)]:
        if digest_file(path)!=expected[name]:raise ValueError('Installed catalogue changed during read-only export')
    atomic_json(output/'preflight.json',report)
    print(json.dumps(dict(output=str(output),target_parts=len(rows),uniform_width_ratio=float(ratios[0]),
                         directly_mapped_soft_palate_parts=len(soft_palate),
                         elapsed_seconds=report['elapsed_seconds'],geometry_fit_created=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
