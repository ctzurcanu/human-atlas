"""Bind exact post-forearm catalogue after literal native-family comparison."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader


def main():
    began=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=base/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'source-binding-v19';assert not O.exists();O.mkdir()
    old=base/'male-cystic-artery-correction-proposal-v43/proposed-atlas.json';current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';prior='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69';expected='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
    assert digest_file(old)==prior and digest_file(current)==expected;a=json.loads(old.read_bytes());raw=current.read_bytes();b=json.loads(raw);ap={p['id']:p for p in a['parts']};bp={p['id']:p for p in b['parts']};assert b['chunks'][:len(a['chunks'])]==a['chunks']
    plan=json.loads((base/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());inputs=json.loads((ROOT/plan['source_inputs']).read_bytes());e=json.loads((D/'feasibility-v2/evidence.json').read_bytes());wanted=set(r['part_id']for r in inputs['fixed_sources']if not r.get('missing_native_source'));wanted.update(r['part_id']for r in e['context']);assert all(ap[p]==bp[p]for p in wanted)
    snapshot=O/'atlas-snapshot.json';snapshot.write_bytes(raw);assert digest_file(snapshot)==expected;reader=AtlasReader();geometries={r['part_id']:r['target_geometry_sha256']for r in e['context']};measured=[]
    for pid,sha in geometries.items():
        m=reader.read_part(b,bp[pid],frame=plan['frame'],catalogue_path=snapshot,catalogue_sha=expected);assert m.geometry_sha==sha;measured.append(dict(part_id=pid,geometry_sha256=sha,source_record_literal=True))
    assert digest_file(current)==expected;result=dict(schema='literal-facial-native-source-rebase/v19',previous_catalogue=str(old.relative_to(ROOT)),previous_catalogue_sha256=prior,current_catalogue=str(current.relative_to(ROOT)),current_catalogue_sha256=expected,snapshot=str(snapshot.relative_to(ROOT)),snapshot_sha256=expected,literal_native_family_and_context_records=sorted(wanted),literal_record_count=len(wanted),literal_context_geometries=measured,unchanged_chunk_prefix=True,candidate_geometry_unchanged=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);print(json.dumps({k:v for k,v in result.items()if k not in{'literal_native_family_and_context_records','literal_context_geometries'}}))


if __name__=='__main__':main()
