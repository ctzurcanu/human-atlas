"""Read registered bounds/metadata only; no triangle or contact testing."""
from pathlib import Path
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.db import Database

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
REVISION='31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2'


def main():
    assert digest_file(SOURCE)==EXPECTED
    atlas=json.loads(SOURCE.read_bytes());parts={p['id']:p for p in atlas['parts']}
    previous_path=BASE/'feasibility-v2/evidence.json';previous=json.loads(previous_path.read_bytes())
    tested={(r['branch'],r['part_id']) for r in previous['context']}
    db=Database();db.assert_current_catalogues('local-ta98');branches=[];rows=[]
    tongue_ids=set(atlas['finiteTongueIntegration']['partIds'])
    def classification(p):
        pid=p['id'];url=p.get('provenance',{}).get('url','')
        if pid.startswith('TA98TONGUE:male:') and pid in tongue_ids and url=='https://github.com/Z-Anatomy/Models-of-human-anatomy':
            return 'native-derived-tongue','Explicit current finiteTongueIntegration part, sourcePartId ZA:Tongue; world coordinates include its recorded geometryOrigin once.'
        if pid.startswith('ZA:') and url=='https://github.com/Z-Anatomy/Models-of-human-anatomy':
            return 'native-ZA','Original baked native source or retained native source ID; current registry metadata equals the actual catalogue record.'
        if pid.startswith('TA98PART:ZA:') and url=='https://github.com/Z-Anatomy/Models-of-human-anatomy':
            return 'native-derived-TA98PART','Literal ZA source-derived partition identity and provenance; current actual registered record is used.'
        if pid.startswith('TA98ARTERIAL:male:') and url=='https://github.com/Z-Anatomy/Models-of-human-anatomy':
            return 'native-derived-installed-arterial','Current installed arterial derivative with native ZA provenance; not excluded because it lacks the ZA prefix.'
        if p.get('license')=='unverified-local-use-only':return 'unregistered-local-overlay','Local Zygote source; source alignment/overlap review is pending, not certified native registration.'
        if pid.startswith('FJ') or url=='https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html':return 'independent-BP3D-overlay','Licensed source, but catalogue position alone does not establish native registration.'
        if pid.startswith('O3M:'):return 'independent-O3M-overlay','Separate imported source frame; exact catalogue occupancy retained as a candidate, not certified native registration.'
        return 'registration-scope-unresolved','Registered in the current numerical model frame; anatomical registration not established solely by that frame ID.'
    with db.connect() as c:
        current=c.execute('SELECT a.catalogue_sha256,c.path,c.metadata FROM atlas_geom.current_catalogues a JOIN atlas_geom.catalogues c ON c.sha256=a.catalogue_sha256 WHERE a.model_id=%s',('local-ta98',)).fetchone()
        assert current['catalogue_sha256']==REVISION and current['metadata']['json_sha256']==EXPECTED
        for branch in ['ascending-palatine','tonsillar','glandular-0','glandular-1']:
            path=BASE/(branch+'-fields.npz')
            with np.load(path,allow_pickle=False) as z:
                fields=z['fields'];assert fields.shape[0]==2 and fields.shape[-1]==3
                low=fields.min(axis=(0,1));high=fields.max(axis=(0,1));frame=str(z['frame'])
            assert frame=='model:local-ta98:assembled-metres'
            branches.append(dict(branch=branch,source_fields=str(path.relative_to(ROOT)),source_fields_sha256=digest_file(path),
                swept_bounds_m=[low.tolist(),high.tolist()],method='Exact min/max of both saved paired source fields; no hull, displacement, contact or occupancy test.'))
            hits=c.execute('''SELECT m.id,m.part_id,m.name,m.frame_id,m.geometry_sha256,m.attributes_sha256,m.metadata,m.source,m.hidden,m.suppressed,m.error
              FROM atlas_geom.meshes m WHERE m.catalogue_sha256=%s AND NOT m.suppressed AND m.geometry_sha256 IS NOT NULL AND m.error IS NULL
              AND m.bbox && cube(%s::float8[],%s::float8[]) ORDER BY m.part_id''',(REVISION,low.tolist(),high.tolist())).fetchall()
            for hit in hits:
                p=parts[hit['part_id']];assert not p.get('suppressed') and hit['metadata']['record']==p
                assert hit['frame_id']==frame
                category,evidence=classification(p)
                row=dict(branch=branch,part_id=p['id'],name=p['name'],registry_mesh_id=hit['id'],registry_revision=REVISION,
                    geometry_sha256=hit['geometry_sha256'],attributes_sha256=hit['attributes_sha256'],frame=frame,
                    catalogue_record_sha256=hashlib.sha256(json.dumps(p,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                    catalogue_record_matches_registry=True,source_category=category,source_frame_evidence=evidence,
                    in_previous_context=(branch,p['id']) in tested,hidden=hit['hidden'],suppressed=False,
                    registered_bounds_intersect_swept_bounds=True,geometry_tested_here=False,
                    source_id=p.get('sourceId'),provenance=p.get('provenance'),license=p.get('license'),
                    geometryOrigin=p.get('geometryOrigin'),sourceOffset=p.get('sourceOffset'),geometryInstance=p.get('geometryInstance'),
                    source_bank_sha256=hit['source'].get('bank_sha256'),source_catalogue_path=hit['source'].get('catalogue_path'))
                rows.append(row)
    native=lambda r:r['source_category'].startswith('native-')
    omitted=[r for r in rows if native(r) and not r['in_previous_context']]
    allparts=set(parts);stale=sorted({pid for branch,pid in tested if pid not in allparts or parts[pid].get('suppressed')})
    focus_ids=sorted(pid for pid in allparts if ('styloglossus' in pid.lower() or 'hyoglossus' in pid.lower()) and not parts[pid].get('suppressed'))
    focus=[dict(part_id=pid,branches=[r['branch'] for r in rows if r['part_id']==pid],omitted_branches=[r['branch'] for r in omitted if r['part_id']==pid],
        source_record_sha256=hashlib.sha256(json.dumps(parts[pid],sort_keys=True,separators=(',',':')).encode()).hexdigest()) for pid in focus_ids]
    result=dict(schema='facial-current-derivative-context-inventory/v1',source_catalogue=str(SOURCE.relative_to(ROOT)),source_catalogue_sha256=EXPECTED,
        registry_revision=REVISION,registry_catalogue_metadata_json_sha256=current['metadata']['json_sha256'],
        previous_context=str(previous_path.relative_to(ROOT)),previous_context_sha256=digest_file(previous_path),previous_context_source_sha256=previous['source_catalogue_sha256'],
        branch_bounds=branches,all_active_registered_overlap_pairs=len(rows),unique_overlapping_parts=len({r['part_id'] for r in rows}),
        native_derived_overlap_pairs=sum(native(r) for r in rows),omitted_native_derived_pairs=len(omitted),omitted_native_derived_parts=sorted({r['part_id'] for r in omitted}),
        omitted_candidates=omitted,unregistered_or_unresolved_candidates=[r for r in rows if not native(r)],all_registered_overlap_candidates=rows,
        explicit_current_styloglossus_hyoglossus=focus,stale_or_suppressed_previous_context_ids=stale,
        no_collision_or_occupied_volume_conclusion=True,no_geometry_testing_construction_import_or_install=True,source_unchanged=True)
    target=BASE/'nearby-derivative-context-inventory-v24.json';atomic_json(target,result);assert digest_file(SOURCE)==EXPECTED
    print(json.dumps({k:result[k] for k in ['source_catalogue_sha256','registry_revision','all_active_registered_overlap_pairs','unique_overlapping_parts','native_derived_overlap_pairs','omitted_native_derived_pairs','omitted_native_derived_parts','explicit_current_styloglossus_hyoglossus','stale_or_suppressed_previous_context_ids']}))


if __name__=='__main__':main()
