"""One bounded metadata-only registry batch for the installed female partial.

Run under geometry_tools.runtime --seconds25. This command never writes an
anatomy catalogue, source bank, frontend file or hierarchy. It indexes at most
1000 already-existing source records in the dedicated local atlas_geom tables.
"""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
OUT=WORK/'female-registry-reimport-v1'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
MALE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
BODY=WORK/'female-core-only-staged-validation-v1/final-classified/served-body.json'
EXPECTED={SOURCE:'8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e',
 MALE:'1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39',
 BODY:'55aabb689f1b66e3f22435fd0fd7d81ba50e54e66cd36de18a55ec87620b5ae1'}


def snapshot():
    result={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in EXPECTED}
    assert all(result[str(p)]==sha for p,sha in EXPECTED.items()),'Source, served body or male catalogue changed'
    for p,count in ((SOURCE,9255),(BODY,3279)):
        a=json.loads(p.read_bytes());m=a['finiteTongueIntegration']
        assert len(a['parts'])==count and m['partialMode']=='female-core-only'
        assert len(m['termGroups'])==14 and len(m['partIds'])==6
    return result


def main(batch):
    assert 6<=batch<=10,'Use an explicit remaining bounded stage'
    began=time.monotonic();report=OUT/f'batch-{batch:02d}.json';proof=OUT/f'batch-{batch:02d}-source-proof.json'
    assert not report.exists()and not proof.exists(),'Never overwrite a completed stage'
    before=snapshot()
    limit=255 if batch==10 else 1000
    result=subprocess.run([sys.executable,'-B','-m','geometry_tools','--wall-seconds','25','--report',str(report),
        'registry-import',str(SOURCE),'--model','local-female-ta98','--max-parts',str(limit)],cwd=ROOT)
    if result.returncode:raise SystemExit(result.returncode)
    after=snapshot();assert before==after
    saved=json.loads(report.read_bytes());assert not saved.get('errors')
    assert saved['parts']==min(batch*1000,9255) and saved['parts_imported_this_run']<=limit
    assert saved['json_sha256']==EXPECTED[SOURCE]
    assert saved['complete']==(batch==10) and saved['current_catalogue_updated']==(batch==10)
    record=dict(schema='female-registry-source-ui-unchanged/v1',batch=batch,source_and_saved_ui_hashes_before=before,
        source_and_saved_ui_hashes_after=after,controls=14,tongue_physical_meshes=6,existing_source_records=9255,saved_served_records=3279,
        registry_only_metadata_index=True,catalogue_frontend_hierarchy_writes=False,male_unchanged=True,batch_record_limit=limit,
        revision_sha256=saved['revision_sha256'],records_indexed=saved['parts'],complete=saved['complete'],
        elapsed_seconds=time.monotonic()-began)
    with proof.open('x')as f:json.dump(record,f,indent=2);f.write('\n')
    print(json.dumps(record))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--batch',type=int,required=True)
    main(p.parse_args().batch)
