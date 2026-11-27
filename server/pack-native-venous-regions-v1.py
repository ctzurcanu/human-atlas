"""Partition literal native forearm-vein source faces; no geometry fitting."""
from pathlib import Path
import copy
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_continuity import boundary_loops

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-native-venous-regions-v1'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'


def main():
    began=time.monotonic();assert digest_file(SOURCE)==EXPECTED
    plan=json.loads((BASE/'source-plan.json').read_bytes());assert plan['source_catalogue_sha256']==EXPECTED
    original=json.loads(SOURCE.read_bytes());atlas=copy.deepcopy(original)
    census=json.loads((ROOT/plan['inventory']).read_bytes());authority={r['id']:r for r in census['terms']}
    out=BASE/'literal-forearm-pack-v1';assert not out.exists();out.mkdir()
    regions=[r for r in plan['regions'] if r['term'] in ['A12.3.08.021','A12.3.08.022']];assert len(regions)==4
    # Save authority and interface contract BEFORE the physical partition.
    terms=set(['A12.3.08.014','A12.3.08.015','A12.3.08.018','A12.3.08.021','A12.3.08.022'])
    def under(term,parent):
        seen=set()
        while term and term not in seen:
            if term==parent:return True
            seen.add(term);term=authority.get(term,{}).get('parent_id')
        return False
    terms|={t for t in authority if under(t,'A12.3.08.015') or under(t,'A12.3.08.018')}
    inventory=dict(schema='native-forearm-venous-region-contract/v1',source_catalogue_sha256=EXPECTED,
        authority_inventory=plan['inventory'],authority_inventory_sha256=plan['inventory_sha256'],
        terms=[{k:authority[t].get(k) for k in ['id','name_en','name_la','parent_id','representation_status','active_exact_source_part_ids']} for t in sorted(terms)],
        whole_assemblies=[dict(whole_concept_id=r['source_part_id'],whole_term=r['source_exact_term'],region_term=r['term'],side=r['side'],
            planned_parts=['literal source-face remainder','literal forearm outer-surface source region'],
            interfaces=r['shared_interfaces'],boundary_status=r['boundary_status']) for r in regions],
        entity_type='Named macroscopic venous-tree segments represented by literal native outer-surface patches; NOT completed finite wall volumes.',
        geometry_change='Original faces partitioned once; original position, normal and UV words copied literally. No changed contacts, no new caps, no movement.',
        source_surface_open_ends_retained=True,lumen_patency_unverified=True,
        held_azygos=dict(term='A12.3.07.002',reason='111 inherited actual right-main-bronchus crossing face pairs; not proposed for acceptance or installation.'),
        primary_references=plan['primary_references'],installed=False)
    inventory_file=ROOT/'reports/ta98-male/native-forearm-venous-region-plan-v1.json';atomic_json(inventory_file,inventory)
    reader=AtlasReader();parts={p['id']:p for p in atlas['parts']};concepts={c['id']:c for c in atlas['concepts']}
    bank=bytearray();proofs=[];groups=[];newparts=[];newconcepts=[]
    def append(data):
        while len(bank)%4:bank.append(0)
        offset=len(bank);bank.extend(data);return offset
    chunk=len(atlas['chunks']);bank_path=out/'geometry.bin'
    for region in regions:
        pid=region['source_part_id'];source=plan['sources'][pid];p=parts[pid];assert p==source['record']
        assert p.get('geometryScale',1)==1 and not any(k in p for k in ['geometryOrigin','sourceOffset','geometryInstance','normalSigns'])
        raw_path=reader.asset(original['chunks'][p['chunk']]['url']);assert digest_file(raw_path)==source['bank_sha256']
        raw=raw_path.read_bytes();n=p['vertexCount'];f=np.frombuffer(raw,dtype='<u4',count=p['indexCount'],offset=p['indices']).reshape(-1,3).copy()
        selected=np.asarray(region['source_face_ids'],dtype=np.int64);remainder=np.setdiff1d(np.arange(len(f)),selected)
        assert len(selected)>0 and len(remainder)>0 and not np.intersect1d(selected,remainder).size
        assert np.array_equal(np.sort(np.r_[selected,remainder]),np.arange(len(f)))
        positions=raw[p['positions']:p['positions']+n*12];normals=raw[p['normals']:p['normals']+n*6];uvs=raw[p['uvs']:p['uvs']+n*8]
        offsets=dict(positions=append(positions),normals=append(normals),uvs=append(uvs));world=np.frombuffer(positions,dtype='<f4').reshape(-1,3)
        rid=f'TA98PART:{pid}:{region["term"]}';remainder_concept=f'TA98SOURCE:{pid}:other-course'
        common=dict(chunk=chunk,**offsets,vertexCount=n,nativeVenousRegionRevision='male-native-forearm-source-regions-v1',
                    nativeVenousSourcePartId=pid,nativeVenousSourceGeometrySha256=source['geometry_sha256'],sourceOuterSurfaceOnly=True,
                    sourceLumenPatencyUnverified=True,anatomicalBoundaryStatus=region['boundary_status'])
        # The original concept stays the whole assembly and remains the target
        # of old direct selections/relations. Its retained physical source ID
        # owns only the non-forearm source faces, never a duplicate envelope.
        p.update(common,indices=append(f[remainder].astype('<u4').tobytes()),indexCount=int(3*len(remainder)),
                 bounds=[world[np.unique(f[remainder])].min(0).astype(float).tolist(),world[np.unique(f[remainder])].max(0).astype(float).tolist()],
                 nativeVenousRegionRole='remaining-source-course')
        q=copy.deepcopy(source['record']);q.update(common,id=rid,conceptId=rid,name=region['name'].capitalize()+f' ({region["side"]})',
                 sourceId=source['record']['sourceId']+':forearm',ta98Term=region['term'],ta98Kind='exact',
                 indices=append(f[selected].astype('<u4').tobytes()),indexCount=int(3*len(selected)),
                 bounds=[world[np.unique(f[selected])].min(0).astype(float).tolist(),world[np.unique(f[selected])].max(0).astype(float).tolist()],
                 nativeVenousRegionRole='forearm-source-region')
        newparts.append(q)
        whole=concepts[pid];whole.update(elements=[pid,rid],componentElements=[pid,rid],ta98Term=region['source_exact_term'],ta98Kind='exact',
            description='Native outer-surface source representation with a separately selectable forearm course; lumen and named-boundary precision remain provisional.')
        newconcepts.extend([dict(id=rid,name=q['name'],elements=[rid],ta98Term=region['term'],ta98Kind='exact',partOf=pid,
            description='Literal forearm portion of the native outer venous surface. Boundary stations are provisional radial-head/distal-radius landmarks; no new finite wall or blood lumen claimed.'),
            dict(id=remainder_concept,name=f'Other source course ({region["side"]})',elements=[pid],partOf=pid,
            description='Remaining literal source surface outside the separately selectable forearm course; this disconnected source selection is not an additional anatomical vessel.')])
        compact=np.asarray(source['source_vertices']);original_compact=Mesh.load(ROOT/source['mesh'])
        assert np.array_equal(compact[original_compact.faces],f)
        loops=boundary_loops(f[selected]);assert len(loops)==2
        shared=[]
        for interface in region['shared_interfaces']:
            ids=compact[np.asarray(interface['source_vertex_ids'])]
            assert set(ids).issubset(set(np.unique(f[selected]))) and set(ids).issubset(set(np.unique(f[remainder])))
            shared.append(dict(source_ring_station=interface['source_ring_station'],source_vertex_ids=ids.tolist(),position_words_sha256=hashlib.sha256(world[ids].tobytes()).hexdigest(),
                normal_words_sha256=hashlib.sha256(np.frombuffer(normals,dtype='<i2').reshape(-1,3)[ids].tobytes()).hexdigest(),same_bank_attributes=True,added_cap_triangles=0))
        atomic_npz(out/(region['term']+'-'+region['side']+'-ancestry.npz'),source_face_ids=np.arange(len(f)),forearm_source_faces=selected,remainder_source_faces=remainder,
                   source_raw_faces=f,source_raw_vertex_ids=compact)
        groups.append(dict(wholeConceptId=pid,wholeTerm=region['source_exact_term'],remainderPartId=pid,remainderConceptId=remainder_concept,
            regionPartId=rid,regionConceptId=rid,regionTerm=region['term'],side=region['side'],sourceId=source['record']['sourceId'],
            sourceBankSha256=source['bank_sha256'],sourceGeometrySha256=source['geometry_sha256'],wholeName=whole['name'],
            sharedInterfaces=shared,sourceFaceCount=len(f),forearmFaceCount=len(selected),remainderFaceCount=len(remainder)))
        proofs.append(dict(source_part_id=pid,source_face_count=len(f),forearm_faces=len(selected),remainder_faces=len(remainder),
            source_positions_sha256=hashlib.sha256(positions).hexdigest(),source_normals_sha256=hashlib.sha256(normals).hexdigest(),source_uvs_sha256=hashlib.sha256(uvs).hexdigest(),
            source_face_partition_exact=True,source_positions_normal_uv_words_literal=True,shared_boundary_loops=shared,
            no_added_removed_or_moved_source_triangles=True,all_original_contacts_identical=True,lumen_patency_unverified=True))
    bank_path.write_bytes(bank);bank_sha=digest_file(bank_path)
    atlas['chunks'].append(dict(url='/local-models/'+str(bank_path.relative_to(ROOT/'.local-models')),bytes=len(bank),sha256=bank_sha))
    atlas['parts'].extend(newparts);atlas['concepts'].extend(newconcepts);atlas['partsCount']=len(atlas['parts'])
    assert atlas['triangles']==original['triangles']
    atlas['nativeVenousRegionIntegration']=dict(revision='male-native-forearm-source-regions-v1',sourceCatalogueSha256=EXPECTED,bankSha256=bank_sha,groups=groups,
        referenceOnly=True,lumenPatencyUnverified=True,anatomicalAcceptance=False,inventory=str(inventory_file.relative_to(ROOT)))
    atlas.setdefault('ta98TermBindings',[]).extend(dict(conceptId=p['conceptId'],partIds=[p['id']],term=p['ta98Term'],kind='exact') for p in newparts)
    candidate=out/'proposed-atlas.json';atomic_json(candidate,atlas)
    # Independently decode actual emitted records. Replay each original face's
    # full corner position/normal/UV words, not merely topology statistics.
    emitted=json.loads(candidate.read_bytes());emittedparts={p['id']:p for p in emitted['parts']};ebank=bank_path.read_bytes()
    for region,proof,group in zip(regions,proofs,groups):
        pid=region['source_part_id'];source=plan['sources'][pid]['record'];raw=reader.asset(original['chunks'][source['chunk']]['url']).read_bytes()
        f=np.frombuffer(raw,dtype='<u4',count=source['indexCount'],offset=source['indices']).reshape(-1,3)
        selected=np.asarray(region['source_face_ids']);remain=np.setdiff1d(np.arange(len(f)),selected)
        for part,faceids in [(emittedparts[pid],remain),(emittedparts[group['regionPartId']],selected)]:
            ef=np.frombuffer(ebank,dtype='<u4',count=part['indexCount'],offset=part['indices']).reshape(-1,3);assert np.array_equal(ef,f[faceids])
            for attr,width,bytesper in [('positions',3,4),('normals',3,2),('uvs',2,4)]:
                count=source['vertexCount']*width*bytesper
                assert ebank[part[attr]:part[attr]+count]==raw[source[attr]:source[attr]+count]
            m=reader.read_part(emitted,part,frame='model:local-ta98:assembled-metres',catalogue_path=candidate,catalogue_sha=digest_file(candidate));assert np.isfinite(m.vertices).all()
        proof['actual_emitted_bank_replay_passed']=True
    unchanged=[p for p in original['parts'] if p['id'] not in {r['source_part_id'] for r in regions}]
    assert all(p==emittedparts[p['id']] for p in unchanged)
    assert digest_file(SOURCE)==EXPECTED
    receipt=dict(schema='native-forearm-venous-source-partition-proof/v1',source_catalogue_sha256=EXPECTED,candidate=str(candidate.relative_to(ROOT)),candidate_sha256=digest_file(candidate),
        geometry_bank=str(bank_path.relative_to(ROOT)),bank_sha256=bank_sha,bank_bytes=len(bank),proofs=proofs,added_distinct_terms=['A12.3.08.021','A12.3.08.022'],
        unchanged_records=len(unchanged),changed_existing_records=4,added_records=4,physical_triangle_count_unchanged=True,
        all_skeleton_heart_lingual_cystic_other_parts_literal=True,azygos_unchanged_held=True,installed=False,controls_validation_pending=True,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'literal-partition-proof.json',receipt)
    print(json.dumps({k:receipt[k] for k in ['candidate','candidate_sha256','bank_bytes','unchanged_records','added_distinct_terms','elapsed_seconds']}))


if __name__=='__main__':main()
