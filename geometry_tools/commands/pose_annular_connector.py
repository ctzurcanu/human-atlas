"""Transport an existing complete annular conduit along a declared guide prefix.

Both actual source mouths remain literal. Original faces, section shapes and
attributes are retained; this produces a proposal requiring triangle review.
"""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.algorithms.source_profile_axis import source_profile_graph,primary_axis
from geometry_tools.algorithms.tube_sections import arclength,unit,rotation_minimizing_frames
from geometry_tools.commands.couple_hand_arteries import rotation


def guided_course(guide,start,end,start_tangent,end_tangent,stations,*,retain_guide_corners=False):
    s=arclength(guide);s/=s[-1]
    if retain_guide_corners:
        corners=[guide[0]]
        for a,b,c in zip(guide[:-2],guide[1:-1],guide[2:]):
            reach=min(.002,np.linalg.norm(b-a)*.25,np.linalg.norm(c-b)*.25)
            corners.extend([b-unit(b-a)*reach,b,b+unit(c-b)*reach])
        corners.append(guide[-1]);points=np.array(corners)
        knots=arclength(points);knots/=knots[-1]
    else:
        knots=np.linspace(0.,1.,5)
        points=np.column_stack([np.interp(knots,s,guide[:,k])for k in range(3)])
    points[0],points[-1]=start,end
    dense_t=unit(np.gradient(guide,s,axis=0))
    tangent=unit(np.column_stack([np.interp(knots,s,dense_t[:,k])for k in range(3)]))
    if retain_guide_corners:tangent=unit(np.gradient(points,knots,axis=0))
    tangent[0],tangent[-1]=start_tangent,end_tangent
    length=arclength(guide)[-1]+np.linalg.norm(start-guide[0])+np.linalg.norm(end-guide[-1])
    curve=CubicHermiteSpline(knots,points,tangent*length)
    t=np.linspace(0.,1.,513);sample=curve(t);arc=arclength(sample);arc/=arc[-1]
    at=np.interp(stations,arc,t)
    path,derivative=curve(at),unit(curve.derivative()(at))
    path[0],path[-1]=start,end
    derivative[0],derivative[-1]=start_tangent,end_tangent
    return path,derivative


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','candidate','part','guide-report','guide-part','output']:
        p.add_argument('--'+key,required=True)
    p.add_argument('--guide-chain',type=int,required=True)
    p.add_argument('--guide-end-m',type=float,required=True)
    p.add_argument('--free-corridor-bones')
    p.add_argument('--free-corridor-skin')
    p.add_argument('--avoid-unjoined-arteries',action='store_true')
    p.add_argument('--recipient-context',help='Hash-bound full fixed recipient surfaces from the actual contact audit')
    p.add_argument('--fixed-prefix-sections',type=int,default=1,
                   help='Retain this many original proximal sections as the existing source junction collar')
    args=p.parse_args();started=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in ['source','candidate','guide-report']}
    data={k:json.loads(v.read_text())for k,v in paths.items()}
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new annular connector stage')
    src=next(r for r in data['source']['rows']if r['partId']==args.part)
    row=next(r for r in data['candidate']['rows']if r['partId']==args.part)
    guide=next(r for r in data['guide-report']['rows']if r['partId']==args.guide_part)
    for r in [src,row,guide]:
        if sha(r['path'])!=r['sha256']:raise ValueError('Source or declared guide changed')
    original,current=Mesh.load(src['path']),Mesh.load(row['path'])
    if not np.array_equal(original.faces,current.faces):raise ValueError('Source face ancestry changed')
    g=source_profile_graph(original.vertices,original.faces,
          triangle_pairing=np.arange(len(original.faces)//2*2).reshape(-1,2),extend_annuli=True)
    old,axis_proof=primary_axis(g,original.faces);nodes=axis_proof['sourceNodes']
    if set(nodes)!=set(range(len(g['profiles'])))or len(g['boundary_nodes'])!=2:
        raise ValueError('This adapter requires one complete unbranched annular source conduit')
    rings=[g['profiles'][i]for i in nodes]
    _,rep,phys=np.unique(original.vertices,axis=0,return_index=True,return_inverse=True)
    if sorted(np.concatenate([phys[s]for s in rings]).tolist())!=list(range(len(rep))):
        raise ValueError('Every original source control must belong to exactly one real section')
    normals=[]
    old_t=unit(np.gradient(old,axis=0))
    for section,t in zip(rings,old_t):
        v=original.vertices[section];v=v-v.mean(0)
        n=unit(np.cross(v,np.roll(v,-1,axis=0)).sum(0))
        normals.append(n if n@t>=0 else -n)
    normals=np.array(normals)
    before=rotation_minimizing_frames(old,normals,
            initial_reference=original.vertices[rings[0][0]]-old[0])
    ends=[rotation(original.vertices[s],current.vertices[s])for s in [rings[0],rings[-1]]]
    target_ends=[current.vertices[s].mean(0)for s in [rings[0],rings[-1]]]
    if not 1<=args.fixed_prefix_sections<len(rings)-1:
        raise ValueError('Keep a finite original collar and at least two movable sections')
    first=args.fixed_prefix_sections-1
    if first:
        if not np.array_equal(current.vertices[rings[0]],original.vertices[rings[0]]):
            raise ValueError('The declared fixed original source collar has moved')
        target_ends[0]=old[first];ends[0]=np.eye(3)
    with np.load(guide['path'],allow_pickle=False)as z:
        route=z['centres'][guide['chains'][args.guide_chain]].copy()
    if route[0,1]<route[-1,1]:route=route[::-1]
    arc=arclength(route)
    if not 0<args.guide_end_m<arc[-1]:raise ValueError('The declared guide prefix must be finite and internal')
    end=np.array([np.interp(args.guide_end_m,arc,route[:,k])for k in range(3)])
    route=np.r_[route[arc<args.guide_end_m],end[None]]
    corridor_proof=None
    if args.free_corridor_bones or args.free_corridor_skin:
        if not args.free_corridor_bones or not args.free_corridor_skin:
            raise ValueError('The free corridor needs both fixed bone and actual skin references')
        from geometry_tools.algorithms.free_space_route import free_space_route
        from geometry_tools.open3d_backend import SurfaceIndex
        from geometry_tools.commands.check_hand_network import combine
        bp,sp=Path(args.free_corridor_bones).resolve(),Path(args.free_corridor_skin).resolve()
        br,sr=json.loads(bp.read_text()),json.loads(sp.read_text())
        bone_rows=[v for v in br['rows']if v['partId'].endswith('.l'if row['side']=='left'else'.r')]
        for v in [*bone_rows,sr]:
            if sha(v['path'])!=v['sha256']:raise ValueError('Fixed corridor source changed')
        if not sr['defaultHiddenExcluded']:raise ValueError('Actual visible skin required')
        radius=max(float(np.linalg.norm(original.vertices[s]-old[i],axis=1).max())for i,s in enumerate(rings))
        obstacles=[];obstacle_ids=[]
        if args.avoid_unjoined_arteries:
            joined={peer for pair in data['source']['pairs']if args.part in [pair['firstPart'],pair['secondPart']]
                         for peer in [pair['firstPart'],pair['secondPart']]}
            others=[v for v in data['candidate']['rows']if v['side']==row['side']and v['partId']not in joined]
            for v in others:
                if sha(v['path'])!=v['sha256']:raise ValueError('Fixed unjoined artery changed')
            if others:obstacles=[SurfaceIndex(combine([Mesh.load(v['path'])for v in others]))]
            obstacle_ids=[v['partId']for v in others]
        if args.recipient_context:
            rp=Path(args.recipient_context).resolve();context=json.loads(rp.read_text())
            if sha(context['sourceCatalogue'])!=context['sourceCatalogueSha256']:
                raise ValueError('Recipient catalogue changed after context extraction')
            for v in context['rows']:
                if sha(v['path'])!=v['sha256']:raise ValueError('Recipient obstacle changed')
            obstacles.append(SurfaceIndex(combine([Mesh.load(v['path'])for v in context['rows']])))
            paths['recipientContext']=rp
        route,corridor_proof=free_space_route(*target_ends,
            SurfaceIndex(combine([Mesh.load(v['path'])for v in bone_rows])),
            SurfaceIndex(Mesh.load(sr['path'])),route,radius_m=radius,surface_obstacles=obstacles)
        corridor_proof['unjoinedArteryIds']=obstacle_ids
        if args.recipient_context:corridor_proof['recipientPartIds']=[v['partId']for v in context['rows']]
        paths.update(corridorBones=bp,corridorSkin=sp)
    station=arclength(old[first:]);station/=station[-1]
    path,tangents=guided_course(route,*target_ends,ends[0]@normals[first],ends[1]@normals[-1],station,
                               retain_guide_corners=corridor_proof is not None)
    after=rotation_minimizing_frames(path,tangents,
          initial_reference=ends[0]@before[first,:,0],final_reference=ends[1]@before[-1,:,0])
    q=current.vertices.copy()
    for i,section in enumerate(rings):
        if i<first:
            for raw in section:q[phys==phys[raw]]=original.vertices[raw]
            continue
        R=after[i-first]@before[i].T
        for raw in section:
            same=phys==phys[raw]
            q[same]=(original.vertices[same]-old[i])@R.T+path[i-first]
    # Literal endpoint replay takes precedence over floating rotation arithmetic.
    for s in [rings[0],rings[-1]]:
        for raw in s:q[phys==phys[raw]]=current.vertices[raw]
    error=max(float(np.max(np.abs(np.linalg.norm(q[s,None]-q[s][None],axis=2)-
                    np.linalg.norm(original.vertices[s,None]-original.vertices[s][None],axis=2))))for s in rings)
    if error>1e-12:raise ValueError('An original source section changed calibre')
    n=np.zeros((len(rep),3));t=q[current.faces]
    np.add.at(n,phys[current.faces.ravel()],np.repeat(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),3,axis=0))
    n=unit(n)[phys]
    attrs=dict(current.attributes,normals=np.rint(n*32767).astype('<i2'),normalSigns=np.ones(len(q),'i1'))
    result=Mesh(q,current.faces,current.frame,attrs,current.source)
    out.mkdir();file=out/'connector.npz';result.save(file)
    report=copy.deepcopy(data['candidate']);target=next(r for r in report['rows']if r['partId']==args.part)
    target.update(path=str(file),sha256=sha(file),geometrySha256=result.geometry_sha)
    report.update(complete=False,installed=False,installationAllowed=False,skeletonChanged=False,
          annularGuideTransport=dict(partId=args.part,sourceSections=len(rings),sourceFacesLiteral=True,
              sourceSectionGaugeError_m=error,mouthCoordinatesLiteral=True,
              guidePart=args.guide_part,guideChain=args.guide_chain,guideEnd_m=args.guide_end_m,
              freeSpaceCorridor=corridor_proof,routeWaypoints=route.tolist(),
              fixedOriginalPrefixSections=args.fixed_prefix_sections,
              oldAxis=old.tolist(),targetAxis=np.r_[old[:first],path].tolist(),
              targetTangents=np.r_[normals[:first],tangents].tolist()),
          inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
          sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    for pair in report['pairs']:
        if args.part in [pair['firstPart'],pair['secondPart']]:
            pair.update(coordinatesLiteral=True,requiresRecoupling=False,
                        literalCoordinatesAndEffectiveNormals=False,effectiveNormalWordsRequireSync=True)
    if time.monotonic()-started>18:raise TimeoutError('18-second annular guide transport')
    write_json(out/'report.json',report)
    print(json.dumps(report['annularGuideTransport']|dict(oldAxis='recorded',targetAxis='recorded',targetTangents='recorded')))


if __name__=='__main__':main()
