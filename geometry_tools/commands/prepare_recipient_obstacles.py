"""Retain full fixed recipient surfaces implicated by an actual contact audit.

This preserves original source identities and topology classifications. Open
surfaces are not silently converted into solids or anatomical interfaces.
"""
import argparse,json,time
from pathlib import Path
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
from geometry_tools.course_io import sha,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['catalogue','contacts','output']:p.add_argument('--'+k,required=True)
    p.add_argument('--part',action='append',required=True)
    args=p.parse_args();start=time.monotonic();cp,ap=Path(args.catalogue).resolve(),Path(args.contacts).resolve()
    atlas,contacts=json.loads(cp.read_text()),json.loads(ap.read_text());out=Path(args.output).resolve()
    if contacts['catalogueSha256']!=sha(cp)or not contacts['complete']:raise ValueError('Recipient audit is stale or incomplete')
    if out.exists():raise ValueError('Use a new obstacle extraction stage')
    by={r['id']:r for r in atlas['parts']};ids=set()
    for r in contacts['rows']:
        if r['partId']in args.part:ids.update(r['regressions'])
    if not ids:raise ValueError('No measured recipient regression for this source scope')
    out.mkdir();reader=AtlasReader();rows=[]
    for pid in sorted(ids):
        row=by[pid]
        if row.get('suppressed')or pid in atlas.get('defaultHidden',[]):raise ValueError('A retired surface was included in the recipient audit')
        if row['system']in ['arterial','venous','cardiac','skin']:continue
        mesh=reader.read_part(atlas,row,frame='model:local-ta98:assembled-metres',catalogue_path=cp)
        file=out/f'context-{len(rows):02}.npz';mesh.save(file);t=topology(mesh)
        rows.append(dict(partId=pid,name=row['name'],system=row['system'],path=str(file),sha256=sha(file),
            geometrySha256=mesh.geometry_sha,sourceRecord=row,closedOrientedManifold=t['closed_oriented_manifold'],
            components=len(t['components']),openEdges=t['open_edges'],
            role='Fixed existing recipient surface; tissue/cavity ownership and intended vascular entry require classification'))
        if time.monotonic()-start>18:raise TimeoutError('18-second recipient extraction stage')
    write_json(out/'report.json',dict(rows=rows,sourceCatalogue=str(cp),sourceCatalogueSha256=sha(cp),
        contactReport=str(ap),contactReportSha256=sha(ap),sourceParts=args.part,
        anatomicalInterfacesInferred=False,sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(fixedSurfaces=len(rows),closedSurfaces=sum(r['closedOrientedManifold']for r in rows),
                         elapsedSeconds=time.monotonic()-start)))


if __name__=='__main__':main()
