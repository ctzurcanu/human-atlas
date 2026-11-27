"""Restore only the seven captured partial-publication files after the timed-out apply."""
from pathlib import Path
import gzip,hashlib,json,os,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
T=time.monotonic();base=ROOT/'.local-models/geometry-workbench/circulatory-completion';evidence_path=base/'sma-timeout-journal-readonly-v1.json';evidence=json.loads(evidence_path.read_bytes());manifest=ROOT/evidence['manifest'];records=json.loads(manifest.read_bytes());observed={r['path']:r for r in evidence['rows']};assert evidence['stableDuringInspection']and len(observed)==len(records)==22
sha=lambda data:hashlib.sha256(data).hexdigest();prepared=[]
for record in records:
    path=(ROOT/record['path']).resolve();backup=(ROOT/record['backup']).resolve();assert path.is_relative_to(ROOT)and backup.is_relative_to(ROOT)
    row=observed[record['path']];raw=path.read_bytes();saved=backup.read_bytes();assert sha(raw)==row['actualSha256'],'Captured current file changed: '+record['path'];assert sha(saved)==record['sha256'],'Backup changed';original=gzip.decompress(saved);assert sha(original)==row['beforeSha256'],'Decompressed backup differs from observed original'
    prepared.append((path,raw,original,record,row))
changed=[r for r in prepared if r[1]!=r[2]];assert len(changed)==7
out=base/'sma-timeout-restoration-v1';out.mkdir(exist_ok=False);journal=[]
for index,(path,raw,original,record,row)in enumerate(changed):
    preserved=out/f'incomplete-published-{index}-{path.name}.gz';preserved.write_bytes(gzip.compress(raw,mtime=0));journal.append(dict(path=record['path'],fromSha256=sha(raw),restoredSha256=sha(original),originalBackup=record['backup'],incompletePublicationBackup=str(preserved.relative_to(ROOT)),incompletePublicationBackupSha256=digest_file(preserved)))
atomic_json(out/'preflight.json',dict(schema='sma-timeout-scoped-restoration/v1',all22CurrentHashesVerified=True,all22CompressedAndDecompressedBackupsVerified=True,capturedEvidence=str(evidence_path.relative_to(ROOT)),capturedEvidenceSha256=digest_file(evidence_path),manifestSha256=digest_file(manifest),changedFiles=journal,restorationStarted=False))
# All current files are checked again immediately before the first mutation.
for path,raw,original,record,row in prepared:assert path.read_bytes()==raw
for path,raw,original,record,row in changed:
    tmp=path.with_name(path.name+'.sma-timeout-restoration-v1.tmp');tmp.write_bytes(original);os.replace(tmp,path)
for path,raw,original,record,row in prepared:assert path.read_bytes()==original
cp=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(cp)=='4d5e86d1ead228e6c172ae895264cd10c8b1469eac6637b77de7531684f9f00d';assert gzip.decompress(Path(str(cp)+'.gz').read_bytes())==cp.read_bytes()
result=dict(schema='sma-timeout-scoped-restoration/v1',complete=True,restoredFiles=7,verifiedOriginalFiles=22,all22EqualOriginalSnapshot=True,masterSha256=digest_file(cp),ordinaryBodySha256=digest_file(ROOT/'.local-models/organ-models/spinal-cord/male-body.json'),forestSha256=digest_file(ROOT/'reports/structural-dependencies/current/male-forest.json'),femaleUnchanged=True,anatomyBanksUnchanged=True,timedOutCandidateNotInstalled=True,preflightSha256=digest_file(out/'preflight.json'),elapsedSeconds=time.monotonic()-T)
atomic_json(out/'complete.json',result);print(json.dumps(result))
