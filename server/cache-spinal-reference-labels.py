"""Cache only selected official PAM50 data, with commit/hash and byte limits.

No third-party code is run and no atlas geometry is installed. The local data
remain reference material until registration, source rights and anatomy review.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,os,time,urllib.request
from geometry_tools.common import ROOT,atomic_json

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--reference-manifest');ap.add_argument('--white-matter',action='store_true');args=ap.parse_args();start=time.monotonic()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
headers={'User-Agent':'Human-Atlas-local-anatomy-reference-review','Accept':'application/vnd.github+json'}
if args.reference_manifest:commit=json.loads(Path(args.reference_manifest).read_text())['commit']
else:
    req=urllib.request.Request('https://api.github.com/repos/spinalcordtoolbox/PAM50/commits/master',headers=headers)
    with urllib.request.urlopen(req,timeout=10)as response:commit=json.load(response)['sha']
assert len(commit)==40 and all(c in '0123456789abcdef'for c in commit)
files=['README.md','CHANGES.md','atlas/info_label.txt','template/info_label.txt',
    'template/PAM50_spinal_levels.nii.gz','template/PAM50_spinal_midpoint.nii.gz',
    'template/PAM50_cord.nii.gz','template/PAM50_rootlets.nii.gz',
    'atlas/PAM50_atlas_30.nii.gz','atlas/PAM50_atlas_32.nii.gz','atlas/PAM50_atlas_34.nii.gz']
if args.white_matter:files=['README.md','atlas/info_label.txt',*('atlas/PAM50_atlas_'+str(i).zfill(2)+'.nii.gz'for i in range(30))]
def fetch(relative):
    url='https://raw.githubusercontent.com/spinalcordtoolbox/PAM50/'+commit+'/'+relative
    request=urllib.request.Request(url,headers={'User-Agent':headers['User-Agent']})
    with urllib.request.urlopen(request,timeout=10)as response:data=response.read(32*1024*1024+1)
    assert len(data)<=32*1024*1024,'Reference file exceeds the declared byte budget'
    if relative.endswith('.nii.gz'):assert data[:2]==b'\x1f\x8b','Expected actual compressed NIfTI data, not a pointer or HTML'
    p=out/relative;p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.part');tmp.write_bytes(data);os.replace(tmp,p)
    return dict(path=relative,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
with ThreadPoolExecutor(max_workers=3)as pool:records=list(pool.map(fetch,files))
assert sum(r['bytes']for r in records)<=64*1024*1024
atomic_json(out/'manifest.json',dict(repository='https://github.com/spinalcordtoolbox/PAM50',commit=commit,files=records,
    raw_data_only=True,installed=False,local_reference_only=True,source_rights_for_redistribution_verified=False,
    elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(commit=commit,files=len(records),bytes=sum(r['bytes']for r in records),elapsed_seconds=time.monotonic()-start,installed=False)),flush=True)
