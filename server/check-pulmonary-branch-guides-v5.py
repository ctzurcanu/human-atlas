"""Measure source branch courses against separately identified airway guides.

These are geometric correspondence witnesses, not a closest-plane anatomical
classifier and not acceptance of vessel/airway intersections or patency.
"""
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic()
    source=json.loads((BASE/'native-curves-v3/evidence.json').read_text())
    contract=json.loads((BASE/'family-contract-v2.json').read_text())
    path=Path('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')
    assert digest_file(path)==contract['source_catalogue_sha256']
    atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    arrays=np.load(BASE/'native-curves-v3/literal-section-centres.npz',allow_pickle=False)
    target=BASE/'upper-branch-guide-witnesses-v5';target.mkdir()
    rows=[];guides=[];saved={}
    sets={
        'right': [('apical','ZA:Apical segmental bronchus of right lung (BI)','right_apical_bronchus'),
                  ('posterior','ZA:Posterior segmental bronchus of right lung (BII)','right_posterior_bronchus'),
                  ('anterior','ZA:Anterior segmental bronchus of right lung (BIII)','right_anterior_bronchus')],
        'left': [('apical','ZA:Apicoposterior segmental bronchus of left lung (BI+BII)','left_apical_bronchus'),
                 ('posterior','ZA:Apicoposterior segmental bronchus of left lung (BI+BII)','left_posterior_bronchus'),
                 ('anterior','ZA:Anterior segmental bronchus of left lung (BIII)','left_anterior_bronchus')],
    }
    for side,definitions in sets.items():
        selected=[c for c in source['curves']if c['part_id'] in [f'ZA:{s} segmental artery of {side} lung'for s in ['Apical','Anterior','Posterior']]]
        queries=[];slices=[]
        for c in selected:
            center=arrays[c['key']+'_centers']
            # Include endpoints and seven internal course witnesses; no invented path.
            stations=np.unique(np.rint(np.linspace(0,len(center)-1,9)).astype(int))
            begin=len(queries);queries.extend(center[stations]);slices.append((c,stations,begin,len(queries)))
        points=np.asarray(queries)
        saved[side+'_literal_course_query_points']=points
        for segment,za,hra in definitions:
            for role,ident in [('original-ZA-hidden-source-reference',za),('fixed-installed-HRA-bronchial-guide','LOCAL:male:pulmonary:bronchus:'+hra)]:
                part=parts[ident];mesh=reader.read_part(atlas,part,frame='model:male:assembled-metres')
                result,proof=BVH(mesh).closest_points(points,backend='bvh')
                key=side+'_'+segment+'_'+('ZA'if role.startswith('original')else'HRA')
                saved[key+'_distance_m']=result['distance_m'];saved[key+'_nearest_points']=result['points'];saved[key+'_face_ids']=result['face_ids'];saved[key+'_barycentric']=result['barycentric']
                guides.append(dict(key=key,id=ident,role=role,segment=segment,geometry_sha256=mesh.geometry_sha,
                                   license=part.get('license'),provenance=part.get('provenance'),queries=len(points),
                                   source_bank_sha256=reader.file_sha(reader.asset(atlas['chunks'][part['chunk']]['url'])),
                                   nearest_triangle_float64_certificate=proof))
                for c,stations,begin,end in slices:
                    d=result['distance_m'][begin:end]
                    rows.append(dict(curve=c['index'],part_id=c['part_id'],source_component=c['component'],guide_key=key,
                                     native_sample_station_indices=stations.tolist(),minimum_center_to_airway_surface_mm=float(d.min()*1000),
                                     median_center_to_airway_surface_mm=float(np.median(d)*1000),
                                     maximum_center_to_airway_surface_mm=float(d.max()*1000),
                                     component_label_not_assigned_by_distance=True))
    np.savez_compressed(target/'literal-nearest-triangle-witnesses.npz',**saved)
    evidence=dict(schema='upper-pulmonary-native-course-airway-guide-witnesses/v5',
                  source_catalogue_sha256=contract['source_catalogue_sha256'],
                  curve_evidence_sha256=digest_file(BASE/'native-curves-v3/evidence.json'),
                  arrays_sha256=digest_file(target/'literal-nearest-triangle-witnesses.npz'),
                  guides=guides,course_queries=rows,
                  component_assignments=[],anatomy_accepted=False,geometry_changed=False,catalogue_changed=False,
                  limitations=['Original-ZA and installed-HRA guides remain distinct source registrations.',
                               'Centre-to-airway-surface distance does not establish vessel wall/lumen clearance.',
                               'Near a segmental bronchus is insufficient to name ascending/descending takeoff or intersegmental venous identity.'],
                  elapsed_seconds=time.monotonic()-began)
    (target/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(target),source_components=len({r['curve']for r in rows}),guide_sources=len(guides),
                          evidence_sha256=digest_file(target/'evidence.json'),elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__':main()
