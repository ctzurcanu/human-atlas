"""Prepare a guarded female install job only from complete actual target checks.

Never run or apply a job here. A failed common fit cannot enter this workflow.
The ordinary geometry job independently rechecks actual transformed anatomy;
its subsequent surface/section review remains mandatory before installation.
"""
from __future__ import annotations
import argparse
import copy
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from tongue_basal_contract import require

spec=importlib.util.spec_from_file_location('female_tongue_producer',ROOT/'server/build-tongue-female-transfer-candidate.py')
producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)


def read(path):return json.loads(Path(path).read_text())
def bound(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=digest_file(p))


def target_reports(paths):
    """Accept actual finite pair evidence, never a visual or summary assertion."""
    require(paths,'Actual target-context audit reports are required')
    result=[]
    for path in paths:
        doc=read(path)
        require(doc.get('schema')=='female-common-tongue-context-audit/v1'
                and doc.get('all_pairs_complete') is True
                and doc.get('skeleton_modified') is False
                and doc.get('source_geometry_modified') is False,
                'Actual target audit is incomplete or changed fixed/source anatomy')
        require(doc.get('proper_crossing_pairs')==0 and doc.get('coplanar_overlap_pairs')==0
                and not doc.get('failed_pairs'),'Failed female target fit cannot be prepared for installation')
        require(doc.get('fit') and digest_file(doc['fit'])==doc.get('fit_sha256'),'Actual target fit binding changed')
        for row in doc.get('pairs',[]):
            s=row.get('summary',{})
            require(s.get('complete_within_pair_budget') is True and s.get('truncated') is False
                    and s.get('crossing_face_pairs')==0 and s.get('coplanar_overlap_face_pairs')==0
                    and s.get('degenerate_pair_count')==0,
                    'A failed/incomplete actual target triangle pair blocks installation')
            require(row.get('first_geometry_sha256') and row.get('second_geometry_sha256')
                    and row.get('constraint',{}).get('kind')=='disjoint'
                    and isinstance(row['constraint'].get('parameters'),dict),
                    'Actual target pair lacks explicit geometry/constraint bindings')
            require(digest_file(row['witnesses'])==row['witnesses_sha256'],'Actual target witness array changed')
        result.append((Path(path).resolve(),doc))
    return result


def prepare(candidate,paths,output):
    began=time.monotonic();reports=target_reports(paths)
    candidate=candidate.resolve();candidate=candidate/'candidate-atlas.json' if candidate.is_dir() else candidate
    require(candidate.is_relative_to(producer.WORK) and candidate.is_file(),'Use an explicit prepared female candidate')
    require(not output.exists() and output.resolve().is_relative_to(producer.WORK),'Use a new local job preparation stage')
    evidence_path=candidate.with_name('transfer-evidence.json');evidence=read(evidence_path)
    require(evidence.get('schema')=='tongue-female-transfer-geometry/v1' and evidence['candidate']==bound(candidate)
            and evidence.get('fixed_skeleton') is True and evidence.get('shared_charts_and_raw_bank_normals_literal') is True,
            'Female candidate lacks actual literal-bank/world-coordinate evidence')
    require(evidence['source']==bound(producer.MALE) and evidence['target']==bound(producer.FEMALE)
            and digest_file(producer.MALE)==producer.MALE_SHA and digest_file(producer.FEMALE)==producer.FEMALE_SHA,
            'An installed donor or original female source changed')
    before,after=read(producer.FEMALE),read(candidate);bp={p['id']:p for p in before['parts']};ap={p['id']:p for p in after['parts']}
    marker=after['finiteTongueIntegration'];active=set(marker['partIds']);mapping=evidence['part_mapping']
    require(marker.get('sex')=='female' and marker.get('stage')=='partial-local-review'
            and marker.get('anatomicalAcceptance') is False and marker.get('tissueOwnershipClaim') is False
            and len(active)==len(marker['partIds'])==12 and len(marker['termGroups'])==17,
            'Keep the exact partial 17-control/12-carrier scope')
    require(set(ap)-set(bp)==active and not(set(bp)-set(ap)),'Unexpected new or removed female anatomy')
    for pid,part in bp.items():
        require(ap[pid]==(dict(part,suppressed=True) if pid in producer.SUPPRESS else part),'Unrelated female source changed: '+pid)
    require(all(after['materials'][k]==v for k,v in before['materials'].items()),'Unrelated female material changed')
    require(after['chunks'][:len(before['chunks'])]==before['chunks'],'Native female bank changed')
    reader=AtlasReader();meshes={pid:reader.read_part(after,ap[pid],frame=producer.FRAME) for pid in active|set(producer.REQUIRED_FIXED)}
    for row in evidence['actual_world_meshes']:
        require(digest_file(row['mesh']['path'])==row['mesh']['sha256']
                and meshes[row['part_id']].geometry_sha==row['geometry_sha256'],'Emitted world geometry changed')
    expected={(pid,fixed) for pid in active for fixed in producer.REQUIRED_FIXED};seen=set()
    bindings={str(candidate):digest_file(candidate),str(evidence_path):digest_file(evidence_path),
              str(producer.MALE):producer.MALE_SHA,str(producer.FEMALE):producer.FEMALE_SHA}
    for path,doc in reports:
        bindings[str(path)]=digest_file(path);bindings[doc['fit']]=doc['fit_sha256']
        for row in doc['pairs']:
            pid=mapping.get(row['partId'],row['partId']);target=row['targetPartId'];pair=(pid,target)
            require(pair in expected and pair not in seen,'Duplicate or unexpected actual target pair')
            require(row['first_geometry_sha256']==meshes[pid].geometry_sha
                    and row['second_geometry_sha256']==meshes[target].geometry_sha,
                    'Audit must measure the actual emitted geometryInstance world coordinates, not a stale fitted NPZ')
            seen.add(pair);bindings[row['witnesses']]=row['witnesses_sha256']
    require(seen==expected,'All twelve carriers and eleven fixed target contexts require actual checks')
    constraints=[dict(kind='fixed',first=pid,parameters=dict(geometry_sha256=meshes[pid].geometry_sha)) for pid in producer.REQUIRED_FIXED]
    topo={pid:topology(meshes[pid]) for pid in meshes}
    for first,second in sorted(expected):
        params=dict(occupancy='solid') if topo[first]['closed_oriented_manifold'] and topo[second]['closed_oriented_manifold'] else {}
        constraints.append(dict(kind='disjoint',first=first,second=second,parameters=params))
    intrinsics=[pid for pid in active if ':intrinsic:' in pid]
    constraints.extend(dict(kind='disjoint',first=a,second=b,parameters=dict(occupancy='solid')) for a,b in combinations(sorted(intrinsics),2))
    for pid in sorted(active):
        for group in ap[pid].get('tissueDrawGroups',[]):
            other=group.get('neighbor')
            if not other or pid>other:continue
            require(other in active,'A literal reciprocal chart has no transferred peer')
            def face_ids(a,b):return [f for g in ap[a].get('tissueDrawGroups',[]) if g.get('neighbor')==b for f in range(g['start']//3,(g['start']+g['count'])//3)]
            constraints.append(dict(kind='shared-interface',first=pid,second=other,
                parameters=dict(first_faces=face_ids(pid,other),second_faces=face_ids(other,pid))))
    inventory_path=ROOT/'reports/ta98-tongue-inventory.json';inventory=read(inventory_path)
    require(inventory['term_count']==len(inventory['terms'])==36,'Complete classified authority inventory required')
    authority=inventory['authority'];groups={g['term']:g for g in marker['termGroups']}
    inventory.update(schema='ta98-geometry-inventory/v1',model='local-female-ta98',model_sex='female',
        root_terms=['A05.1.04.001','A04.1.05.001'],source_catalogue=str(producer.FEMALE),source_catalogue_sha256=producer.FEMALE_SHA,
        authority_path=authority['path'],authority_sha256=authority['sha256'],
        authority_corrections_path=authority['corrections_path'],authority_corrections_sha256=authority['corrections_sha256'])
    for row in inventory['terms']:
        group=groups.get(row['authority']['id']);row.update(source_meshes=row['source_meshes']['female'],
            applicable_to_model=row.get('required_in_typical_adult',True),candidate_control=group,
            candidate_control_absent=group is None,candidate_part_ids=group['partIds'] if group else [],anatomical_acceptance=False)
    for path in [inventory_path,Path(__file__),ROOT/'server/build-tongue-female-transfer-candidate.py',ROOT/'ATTRIBUTION.md',Path(authority['path']),Path(authority['corrections_path'])]:bindings[str(path.resolve())]=digest_file(path)
    bindings.update(evidence['source_bank_bindings'])
    require(all(digest_file(path)==sha for path,sha in bindings.items()),'An explicit transfer/audit input changed')
    output.mkdir(parents=True);atomic_json(output/'inventory.json',inventory)
    spec=dict(model='local-female-ta98',proposed_catalogue=str(candidate),expected_catalogue_sha256=producer.FEMALE_SHA,
        changed_parts=sorted(active|set(producer.SUPPRESS)),closed_parts=sorted(active),require_closed=False,
        max_pairs=250000,max_tests=2000000,inventory=str(output/'inventory.json'),interface_checks=constraints,
        evidence_files=sorted(bindings),explicit_input_bindings=bindings,stage='partial-local-review',external_deployment=False,
        assembly=dict(whole_organ=marker['wholeConceptId'],part_ids=marker['partIds'],
            fixed=['Entire skeleton','Every unrelated female source','Entire installed male model'],
            adjoining_interfaces=[dict(kind='source-bound-proper-similarity',scope='One explicit placement of twelve literal-bank carriers; no target cortical seam is inferred'),
                dict(kind='actual-female-context',scope='All132 explicit fixed-target pairs are independently rechecked by this geometry job')],
            remaining_limits=[producer.LIMIT,'Actual emitted-bank surface/section review is required before apply'],
            representation_mode=marker['representationMode'],tissue_ownership_claim=False,anatomical_acceptance=False))
    atomic_json(output/'install-spec.json',spec);atomic_json(output/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec))
    print(json.dumps(dict(output=str(output.resolve()),interface_checks=len(constraints),target_pairs=len(seen),
                         controls=17,physical_meshes=12,installed=False,job_started=False,elapsed_seconds=time.monotonic()-began)))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--candidate',type=Path,required=True)
    ap.add_argument('--context-audit',type=Path,action='append',required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();prepare(args.candidate,args.context_audit,args.output)
