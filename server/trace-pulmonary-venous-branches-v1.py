"""Read-only literal source venous incidence and named airway witnesses.

Names and distances alone do not classify a vein as intra/intersegmental.
The endpoint graph records coordinates of disconnected casts, never patency.
"""
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    started = time.monotonic()
    source = json.loads((BASE/'native-curves-v3/evidence.json').read_text())
    contract = json.loads((BASE/'family-contract-v2.json').read_text())
    path = Path('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')
    assert digest_file(path) == contract['source_catalogue_sha256']
    atlas = json.loads(path.read_text()); parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader(); arrays = np.load(BASE/'native-curves-v3/literal-section-centres.npz', allow_pickle=False)
    output = BASE/'venous-trace-evidence-v1'; output.mkdir()
    veins = [c for c in source['curves'] if c['system'] == 'venous']
    charts = {}; witnesses = []; guide_records = []; root_nodes = []
    for side in ['right', 'left']:
        selected = [c for c in veins if ('of '+side+' lung') in c['part_id']]
        queries = []; slices = []
        for c in selected:
            centers = arrays[c['key']+'_centers']
            at = np.unique(np.rint(np.linspace(0, len(centers)-1, 9)).astype(int))
            begin = len(queries); queries.extend(centers[at]); slices.append((c, at, begin, len(queries)))
        q = np.asarray(queries); charts[side+'_literal_queries'] = q
        guides = [p for p in atlas['parts'] if p['id'].startswith('ZA:') and 'segmental bronchus of '+side+' lung' in p['id']]
        for gi, guide in enumerate(guides):
            mesh = reader.read_part(atlas, guide, frame='model:male:assembled-metres')
            result, proof = BVH(mesh).closest_points(q, backend='bvh')
            key = side+'_guide_'+str(gi)
            for field in ['distance_m', 'points', 'face_ids', 'barycentric']:
                charts[key+'_'+field] = result[field]
            guide_records.append(dict(key=key, id=guide['id'], geometry_sha256=mesh.geometry_sha,
                                      original_source_guide=True, currently_suppressed=guide.get('suppressed'),
                                      nearest_triangle_certificate=proof))
            for c, at, begin, end in slices:
                d = result['distance_m'][begin:end]
                witnesses.append(dict(curve=c['index'], part_id=c['part_id'], component=c['component'], guide_id=guide['id'],
                                      guide_key=key, stations=at.tolist(), distances_mm=(d*1000).tolist(),
                                      median_mm=float(np.median(d)*1000), minimum_mm=float(d.min()*1000), maximum_mm=float(d.max()*1000),
                                      no_identity_assigned_by_distance=True))
        # Literal endpoint-to-other-station incidence, at submicron tolerance.
        for c in selected:
            center = arrays[c['key']+'_centers']; ei = c['larger_radius_endpoint']; p = center[0 if ei == 0 else -1]
            matches = []
            for other in selected:
                if c['index'] == other['index']: continue
                centers = arrays[other['key']+'_centers']; dist = np.linalg.norm(centers-p, axis=1)
                at = int(dist.argmin())
                if dist[at] < 1e-6:
                    matches.append(dict(curve=other['index'], part_id=other['part_id'], component=other['component'],
                                        station=at, distance_m=float(dist[at]),
                                        endpoint=0 if at==0 else 1 if at==len(centers)-1 else None))
            root_nodes.append(dict(curve=c['index'], part_id=c['part_id'], component=c['component'],
                                   larger_radius_endpoint=ei, root_candidate_m=p.tolist(),
                                   literal_station_matches=matches, physiological_drainage_not_established=True))
    # Four closed source casts do not have a recoverable serial flow centreline.
    root_casts=[]
    for name in ['Right superior','Right inferior','Left superior','Left inferior']:
        ident='ZA:'+name+' pulmonary vein'; mesh=reader.read_part(atlas,parts[ident],frame='model:male:assembled-metres')
        side=name.split()[0].lower(); selected=[c for c in veins if 'of '+side+' lung' in c['part_id']]
        q=np.asarray([arrays[c['key']+'_centers'][0 if c['larger_radius_endpoint']==0 else -1]for c in selected])
        result,proof=BVH(mesh).closest_points(q,backend='bvh');key=name.lower().replace(' ','_')+'_cast'
        for field in ['distance_m','points','face_ids','barycentric']:charts[key+'_'+field]=result[field]
        root_casts.append(dict(id=ident,geometry_sha256=mesh.geometry_sha,closed_outer_cast=True,serial_centerline_recovered=False,
                               query_curve_indices=[c['index']for c in selected],larger_end_to_cast_distance_mm=(result['distance_m']*1000).tolist(),
                               witness_key=key,nearest_triangle_certificate=proof,lumen_and_shared_annulus_not_established=True))
    np.savez_compressed(output/'literal-venous-airway-witnesses.npz',**charts)
    evidence=dict(schema='pulmonary-venous-literal-incidence-airway-witnesses/v1',source_catalogue_sha256=contract['source_catalogue_sha256'],
                  curve_evidence_sha256=digest_file(BASE/'native-curves-v3/evidence.json'),arrays_sha256=digest_file(output/'literal-venous-airway-witnesses.npz'),
                  venous_components=len(veins),source_guides=guide_records,course_witnesses=witnesses,literal_root_candidates=root_nodes,root_casts=root_casts,
                  accepted_intra_intersegmental_assignments=[],geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
                  limitations=['Coincident source centre points are a course-graph witness, not a welded annulus or communicating lumen.',
                               'A bronchial course distance does not define an intersegmental plane, tributary territory or physiological drainage.',
                               'The actual source bronchi and current installed HRA bronchi are separate registrations.',
                               'Source names cannot qualify intrasegmental/intersegmental ownership without course and drainage review.'],
                  elapsed_seconds=time.monotonic()-started)
    (output/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(output),vein_components=len(veins),source_guides=len(guide_records),
                          witness_rows=len(witnesses),evidence_sha256=digest_file(output/'evidence.json'),elapsed_seconds=evidence['elapsed_seconds'])))


if __name__=='__main__': main()
