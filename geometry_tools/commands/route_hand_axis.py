"""Search a fixed-context corridor for a hash-bound original hand artery axis.

Only a guide is produced. Free endpoints impose a parent recoupling obligation;
they never authorize detached vessel installation or removal of source detail.
"""
import argparse,json,time
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.commands.check_hand_network import combine
from geometry_tools.open3d_backend import SurfaceIndex,library
from geometry_tools.algorithms.feature_normals import indexed_feature_normals,inward_sheet_clearance
from geometry_tools.algorithms.layered_tube_corridor import layered_tube_corridor
from geometry_tools.algorithms.recipient_clearance_field import RecipientClearanceField
from geometry_tools.algorithms.hybrid_surface_index import HybridSurfaceIndex
from geometry_tools.algorithms.tube_sections import arclength


@lru_cache(maxsize=1)
def fixed_context_indices(bones,skin_record,recipients,margin):
    # One cached context per bounded worker. It is shared by independent
    # finger routes, and cannot accumulate a catalogue history in memory.
    bone=SurfaceIndex(combine([Mesh.load(path)for path,digest in bones]))
    skin=SurfaceIndex(Mesh.load(skin_record[0]))
    context=RecipientClearanceField([(Mesh.load(path),closed,label)for path,digest,closed,label in recipients],margin_m=margin)
    scene=HybridSurfaceIndex(combine([bone.mesh,skin.mesh,context.surface.mesh]))
    return bone,skin,context,scene


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['source-axis-report','course-match-report','bones','skin-reference','recipient-context','output']:
        p.add_argument('--'+k,required=True)
    p.add_argument('--free-endpoints',action='store_true')
    p.add_argument('--grid-mm',type=float,default=.5)
    p.add_argument('--half-width-mm',type=float,default=4.)
    p.add_argument('--recipient-clearance-mm',type=float,default=.02)
    p.add_argument('--transverse-neighbour-cells',type=int,default=1)
    p.add_argument('--native-route-seconds',type=float,default=16.)
    p.add_argument('--allow-lateral-detours',action='store_true')
    args=p.parse_args();start=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in
        ['source-axis-report','course-match-report','bones','skin-reference','recipient-context']}
    data={k:json.loads(v.read_text())for k,v in paths.items()}
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new saved corridor stage')
    ref=data['source-axis-report']
    if sha(ref['axisPath'])!=ref['axisSha256'] or sha(ref['sourceMesh'])!=ref['sourceMeshSha256']:
        raise ValueError('Original source or contour axis changed')
    source_axis=np.load(ref['axisPath'])['axis']
    rec=next(r for r in data['course-match-report']['courseMatches']if
             r['sourcePartId']==ref['partId']and r['sourceComponent']==ref['component'])
    target=np.array(rec['targetAxis']);side='left'if '-left_'in ref['partId']else'right'
    control=target if args.free_endpoints else np.r_[source_axis[:3],target]
    s=arclength(control);t=np.linspace(0,s[-1],int(np.ceil(s[-1]/.001))+1)
    guide=PchipInterpolator(s,control)(t)
    source=Mesh.load(ref['sourceMesh']);mouth=source.vertices[ref['actualMouthRaw']]
    mouth_radius=np.linalg.norm(mouth-mouth.mean(0),axis=1).max()
    # The final cap-centre contour tends to zero; keep its actual source radius
    # from the last recovered level, rather than inventing a thinner endpoint.
    radii=np.array([mouth_radius,*[r['maximumContourRadius_m']for r in ref['proof']['contourWitnesses']]])
    radii=np.r_[radii,radii[-1]]
    radius=np.interp(np.linspace(0,1,len(guide)),np.linspace(0,1,len(radii)),radii)
    bones=[r for r in data['bones']['rows']if r['partId'].endswith('.l'if side=='left'else'.r')]
    fixed=[*bones,data['skin-reference'],*data['recipient-context']['rows']]
    for r in fixed:
        if sha(r['path'])!=r['sha256']:raise ValueError('Fixed recipient arrays changed')
    if not 0<args.recipient_clearance_mm<=1:raise ValueError('Explicit positive recipient margin within 1 mm required')
    bone,skin,context,scene=fixed_context_indices(tuple((r['path'],r['sha256'])for r in bones),
        (data['skin-reference']['path'],data['skin-reference']['sha256']),
        tuple((r['path'],r['sha256'],r['closedOrientedManifold'],r['partId'])for r in data['recipient-context']['rows']),
        args.recipient_clearance_mm*.001)
    queries=dict(points=0,skinRimWitnesses=0)
    def occupied(index,q):
        tensor=library().core.Tensor(((q-index.origin)/index.scale).astype(np.float32))
        return index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
    def gap(points):
        queries['points']+=len(points);result=[]
        for first in range(0,len(points),24000):
            q=points[first:first+24000];hit,_=bone.closest_points(q)
            value=hit['distance_m']*np.where(occupied(bone,q),-1.,1.)-.0002
            hit,_=skin.closest_points(q);normal,rim=indexed_feature_normals(skin,hit)
            _,g,_=inward_sheet_clearance(q,hit,normal,rim)
            queries['skinRimWitnesses']+=int(rim.sum());value=np.minimum(value,g-.0002)
            value=np.minimum(value,context.clearance(q))
            result.append(value)
        return np.concatenate(result)
    def segments(a,b):
        return scene.segment_clear(a,b)
    endpoint_gap=gap(guide[[0,-1]])
    route,proof=layered_tube_corridor(guide,radius,gap,segment_clear=segments,
        grid_m=args.grid_mm*.001,half_width_m=args.half_width_mm*.001,
        fixed_prefix=0 if args.free_endpoints else 2,fixed_distal=not args.free_endpoints,wall_seconds=args.native_route_seconds,
        transverse_neighbour_cells=args.transverse_neighbour_cells,allow_lateral_detours=args.allow_lateral_detours)
    if time.monotonic()-start>22:raise TimeoutError('22-second native source-axis routing stage')
    out.mkdir();file=out/'route.npz'
    np.savez(file,guide=guide,radii=radius,route=np.empty((0,3))if route is None else route,source_axis=source_axis)
    report=dict(proof=proof,partId=ref['partId'],component=ref['component'],
        routePath=str(file),routeSha256=sha(file),endpointGap_m=endpoint_gap.tolist(),
        endpointRadii_m=radius[[0,-1]].tolist(),queryCounts=queries,
        parentsRequireRecoupling=bool(args.free_endpoints),skeletonChanged=False,
        sourceGeometryChanged=False,installed=False,installationAllowed=False,
        openRecipientVolumesNotEstablished=True,continuousContainmentNotEstablished=True,
        bodyClearance_m=.0002,recipientClearance_m=args.recipient_clearance_mm*.001,
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
        sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-start)
    if 'partialRoute_m'in proof:
        last=np.array(proof['partialRoute_m'][-1]);near=np.r_[last[None],np.array(proof['nextFreeStates_m'])]
        report['failureWitnesses']=context.witnesses(near)
    report['queryAcceleration']=dict(nativeFaces=len(context.surface.native_ids),float64FallbackFaces=len(context.surface.fallback_ids),allSourceFacesRetained=True)
    write_json(out/'report.json',report)
    print(json.dumps(dict(component=ref['component'],proof={k:v for k,v in proof.items()if k not in ['partialRoute_m','nextFreeStates_m']},
        endpointGap_m=report['endpointGap_m'],elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
