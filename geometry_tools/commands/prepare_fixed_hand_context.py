"""Export every eligible active fixed hand neighbour, including noncontacting ones.

Bounds only shortlist source records. Full original triangles remain the query
surface; no neighbour is omitted merely because an earlier candidate missed it.
Skin and canonical bones are separate, and insertion-site overlays are recorded
as reference annotations rather than additional physical blocking tissues.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.course_io import sha,write_json
from geometry_tools.open3d_backend import library


def closed_hint(mesh):
    v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces]
    directed=f[:,[0,1,1,2,2,0]].reshape(-1,2)
    edges,owner,counts=np.unique(np.sort(directed,axis=1),axis=0,return_inverse=True,return_counts=True)
    winding=np.bincount(owner,weights=np.where(directed[:,0]<directed[:,1],1.,-1.))
    closed=bool(len(counts)and np.all(counts==2)and np.all(winding==0))
    if closed:
        o3d=library();native=o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v),o3d.utility.Vector3iVector(f))
        closed=bool(native.is_vertex_manifold())
    return closed,dict(boundaryEdges=int(np.sum(counts==1)),nonmanifoldEdges=int(np.sum(counts>2)),
                       inconsistentWindingEdges=int(np.sum((counts==2)&(winding!=0))),sourceSelfIntersectionsChecked=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['catalogue','moving-source','bones','skin-reference','output']:p.add_argument('--'+k,required=True)
    p.add_argument('--side',choices=['left','right'],required=True)
    args=p.parse_args();start=time.monotonic();paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in
        ['catalogue','moving-source','bones','skin-reference']};data={k:json.loads(v.read_text())for k,v in paths.items()}
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new fixed-context inventory stage')
    out.mkdir();atlas=data['catalogue'];hidden=set(atlas['defaultHidden']);moving={r['partId']for r in data['moving-source']['rows']}
    skin=set(data['skin-reference']['sourcePartIds']);bone=set()
    for r in data['bones']['rows']:
        if r['partId'].endswith('.l'if args.side=='left'else'.r'):bone.update(r.get('activePartIds',r['sourcePartIds']))
    low=np.array([.21,.73,-.02]);high=np.array([.39,.95,.15])
    if args.side=='right':low[0],high[0]=-high[0],-low[0]
    reader=AtlasReader();rows=[];exclusions=[]
    for part in atlas['parts']:
        if part['id']in hidden or part.get('suppressed')or part.get('defaultHidden')or not part.get('indexCount'):continue
        box=np.asarray(part.get('bounds',[]))
        if box.shape!=(2,3)or np.any(box[1]<low)or np.any(box[0]>high):continue
        reason=None
        if part['id']in moving:reason='Explicitly moving arterial assembly; independent whole-network replay required'
        elif part['id']in skin:reason='Actual active body skin is queried separately at its independent clearance'
        elif part['id']in bone:reason='Actual canonical whole-bone envelope is queried separately'
        elif part.get('system')=='attachments'and '2: Muscular insertions'in part.get('groups',[]):reason='Source muscular insertion/origin reference overlay; not an extra physical tissue'
        if reason:
            exclusions.append(dict(partId=part['id'],reason=reason));continue
        m=reader.read_part(atlas,part,frame='model:local-ta98:assembled-metres',catalogue_path=paths['catalogue'])
        closed,proof=closed_hint(m);file=out/('fixed-'+str(len(rows)).zfill(3)+'.npz');m.save(file)
        rows.append(dict(partId=part['id'],name=part['name'],system=part['system'],path=str(file),sha256=sha(file),
            geometrySha256=m.geometry_sha,sourceRecord=part,closedOrientedManifold=closed,topologyHint=proof,
            clearanceSidePolicy='initial-control-side'if closed else'unsigned-sheet',
            role='Fixed active recipient tissue; specific vascular entry not inferred from proximity'))
        if time.monotonic()-start>20:raise TimeoutError('20-second native full hand context export')
    report=dict(complete=True,side=args.side,rows=rows,sourceCatalogue=str(paths['catalogue']),
        sourceCatalogueSha256=sha(paths['catalogue']),regionBounds_m=[low.tolist(),high.tolist()],
        shortlistingOnly=True,fullSourceTrianglesRetained=True,exclusions=exclusions,
        earlierContactRequired=False,geometryChanged=False,installed=False,anatomicalAcceptance=False,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json',report)
    print(json.dumps(dict(rows=len(rows),closedHints=sum(r['closedOrientedManifold']for r in rows),
                         exclusions=len(exclusions),elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
