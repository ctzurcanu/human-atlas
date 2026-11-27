"""Write an isolated position/index storage preview using actual atlas origins.

python -m geometry_tools.commands.preview_tube_storage FIT_DIRECTORY --output NEW_DIRECTORY
Normals/materials and the full catalogue are deliberately absent: never install
this geometry-only review preview as an atlas. Source profiles are replayed.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from geometry_tools.course_io import load_arrays,load_course_spec,sha,write_json
from geometry_tools.model import AtlasReader


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fit_directory');parser.add_argument('--output',required=True)
    args=parser.parse_args();started=time.monotonic()
    folder=Path(args.fit_directory).resolve();out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new immutable preview directory')
    fit=json.loads((folder/'run.json').read_text())
    spec,g,source,*_=load_course_spec(fit['specPath'])
    if fit['bindings']['specSha256']!=sha(fit['specPath']) or fit['displacementsSha256']!=sha(folder/fit['displacements']):
        raise ValueError('Stale fit specification or field')
    d=load_arrays(folder/fit['displacements'])['delta']
    if d.shape!=g['centres'].shape or not np.isfinite(d).all():raise ValueError('Finite complete field required')
    atlas=json.loads(Path(spec['catalogue']).read_text());records={r['id']:r for r in atlas['parts']}
    proposed=g['vertices']+d[g['cells']];actual=proposed.copy();rows=[];out.mkdir();bank=bytearray();preview=[]
    for entry in source['records']:
        pid=entry['partId']
        if pid not in spec['moving_parts']:continue
        record=records[pid]
        if record.get('geometryInstance') or record.get('geometryScale',1.)!=1. or np.any(record.get('sourceOffset',[0,0,0])):
            raise ValueError('An explicit alternate encoding pose is required for this source')
        off,n=entry['offset'],entry['vertices']
        q=proposed[off:off+n];origin=np.asarray(record.get('geometryOrigin',[0.,0.,0.]))
        local=np.asarray(q-origin,dtype='<f4')
        f=g['faces'][entry['faceOffset']:entry['faceOffset']+entry['faces']]-off
        pos=len(bank);bank.extend(local.tobytes());indices=len(bank);bank.extend(np.asarray(f,dtype='<u4').tobytes())
        preview.append(dict(id=pid,chunk=0,positions=pos,indices=indices,vertexCount=n,indexCount=int(f.size),geometryOrigin=origin.tolist()))
    binary=out/'positions-indices.bin';binary.write_bytes(bank)
    reader=AtlasReader();url='/local-models/'+str(binary.relative_to(reader.root/'.local-models'))
    catalogue=dict(chunks=[dict(url=url,sha256=sha(binary))],parts=preview,geometryOnlyReview=True,notInstallable=True)
    write_json(out/'preview.json',catalogue)
    byid={r['partId']:r for r in source['records']}
    max_error=0.
    for record in preview:
        entry=byid[record['id']];off,n=entry['offset'],entry['vertices']
        m=reader.read_part(catalogue,record,frame='atlas-world-metres',catalogue_path=out/'preview.json')
        # AtlasReader compacts unused controls; this command requires all raw
        # source controls to be indexed, so correspondence cannot drift silently.
        if len(m.vertices)!=n or m.source['source_vertex_ids']!=list(range(n)):
            raise ValueError('Unused source controls need an explicit raw-storage map')
        actual[off:off+n]=m.vertices
        error=float(np.linalg.norm(m.vertices-proposed[off:off+n],axis=1).max(initial=0));max_error=max(max_error,error)
        path=out/f'part-{len(rows):03}.npz'
        np.savez_compressed(path,vertices=m.vertices,faces=m.faces,source_vertices=g['vertices'][off:off+n])
        rows.append(dict(partId=record['id'],geometryFile=path.name,sha256=sha(path),maximumPositionRoundoff_m=error))
        if time.monotonic()-started>15:raise TimeoutError('15s native storage preview deadline')
    failures=[];max_pair=0.;max_area=0.
    for p in source['profiles']:
        raw=np.asarray(p['rawIds'])+p['offset'];a=g['vertices'][raw];q=actual[raw]
        pair=float(abs(np.linalg.norm(a[:,None]-a[None,:],axis=2)-np.linalg.norm(q[:,None]-q[None,:],axis=2)).max())
        def area(v):
            v=v-v.mean(0);return np.linalg.norm(np.cross(v,np.roll(v,-1,axis=0)).sum(0))*.5
        old=area(a);rel=float(abs(area(q)/old-1))if old else None
        max_pair=max(max_pair,pair)
        if rel is not None:max_area=max(max_area,rel)
        if rel is None or pair>5e-8 or rel>.005:failures.append(dict(partId=p['partId'],cell=p['cell'],pairError_m=pair,relativeAreaError=rel))
    report=dict(schema='tube-storage-preview/v1',specPath=fit['specPath'],fitReportSha256=sha(folder/'run.json'),
        rows=rows,storageCatalogueSha256=sha(out/'preview.json'),storageBankSha256=sha(binary),
        actualAtlasReaderRoundTrip=True,originalSourceProfilesChecked=len(source['profiles']),
        sourceGaugeFailures=failures,maximumSourcePairDistanceError_m=max_pair,
        maximumSourceAreaRelativeError=max_area,maximumPositionRoundoff_m=max_error,
        normalsAndMaterialsAbsent=True,geometryOnlyNotInstallable=True,
        allSurfaceAndBodyChecksPending=True,installed=False,anatomicalAcceptance=False,
        sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items()if k not in ['rows','sourceGaugeFailures']}))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
