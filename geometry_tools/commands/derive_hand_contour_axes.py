"""Recover curved axes of recorded hand source components without changing meshes."""
import argparse,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,primary_axis
from geometry_tools.algorithms.geodesic_contour_axis import geodesic_contour_axis


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source','profiles','binding','output']:p.add_argument('--'+k,required=True)
    p.add_argument('--parts',required=True,help='Comma-separated exact source part IDs')
    args=p.parse_args();start=time.monotonic();paths={k:Path(getattr(args,k)).resolve()for k in ['source','profiles','binding']}
    data={k:json.loads(v.read_text())for k,v in paths.items()};out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new source-axis inventory stage')
    selected=set(args.parts.split(','));rows=[r for r in data['source']['rows']if r['partId']in selected]
    if {r['partId']for r in rows}!=selected:raise ValueError('Unknown source part')
    out.mkdir();records=[]
    for row in rows:
        if sha(row['path'])!=row['sha256']:raise ValueError('Original source arrays changed')
        m=Mesh.load(row['path']);bank=next(r for r in data['profiles']['rows']if r['partId']==row['partId'])
        if sha(bank['path'])!=bank['sha256']:raise ValueError('Original section pairing changed')
        with np.load(bank['path'],allow_pickle=False)as z:pairing=z['retained_source_triangle_pairs'].copy()
        graph=source_profile_graph(m.vertices,m.faces,triangle_pairing=pairing,extend_annuli=True)
        for owner in data['binding']['rows']:
            if owner['partId']!=row['partId']:continue
            if not owner['anatomicallyReviewed']:raise ValueError('Explicit source ownership required')
            fids=np.array(owner['sourceFaceIds']);raw=np.unique(m.faces[fids]);remap=np.full(len(m.vertices),-1,int);remap[raw]=np.arange(len(raw))
            old,prior=primary_axis(graph,m.faces[fids]);mouth=graph['boundary_raw'][graph['boundary_nodes'].index(prior['proximalNode'])]
            try:
                axis,proof=geodesic_contour_axis(m.vertices[raw],remap[m.faces[fids]],remap[mouth],old[-1],spacing_m=.00075,wall_seconds=3.)
            except ValueError as exc:
                records.append(dict(partId=row['partId'],component=owner['physicalComponent'],complete=False,
                    status='SourceAxisUnresolved',error=str(exc),sourceMeshSha256=row['sha256'],geometryModified=False))
                if time.monotonic()-start>20:raise TimeoutError('20-second native source-axis inventory')
                continue
            dest=out/('component-'+str(len(records)).zfill(2));dest.mkdir();file=dest/'axis.npz'
            np.savez(file,axis=axis,prior_axis=old,source_raw_ids=raw,source_face_ids=fids)
            record=dict(complete=True,partId=row['partId'],component=owner['physicalComponent'],actualMouthRaw=mouth.tolist(),
                sourceMesh=row['path'],sourceMeshSha256=row['sha256'],axisPath=str(file),axisSha256=sha(file),
                oldAxisStations=len(old),newAxisStations=len(axis),originalPrimaryAxisProof=prior,
                proof=proof,geometryModified=False,installed=False)
            write_json(dest/'report.json',record);records.append(dict(complete=True,partId=row['partId'],component=owner['physicalComponent'],
                path=str(dest/'report.json'),sha256=sha(dest/'report.json'),oldAxisStations=len(old),newAxisStations=len(axis)))
            if time.monotonic()-start>20:raise TimeoutError('20-second native source-axis inventory')
    report=dict(complete=all(r['complete']for r in records),rows=records,allSourceMeshesLiteral=True,installed=False,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json',report);print(json.dumps(report))


if __name__=='__main__':main()
