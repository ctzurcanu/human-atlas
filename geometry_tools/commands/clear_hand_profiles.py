"""Local section-preserving clearance proposals on reviewed hand source meshes.

Actual source section controls and strict annular face links define the unknowns.
No proximity assignment across digits, radial shrink, new faces or installation.
Open skin normals guide proposals only; actual triangles require separate checks.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.open3d_backend import SurfaceIndex,library
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores,profile_translation_weights
from geometry_tools.algorithms.profile_graph_transport import finite_diagnostic,transport_linked_profiles
from geometry_tools.algorithms.section_course_socp import solve_weighted_section_socp


def section_proposal(vertices,faces,profiles,links,bone,skin,*,trust_m=.002,physical_ids=None,whole_source=False):
    """Translate whole sections; account for section vertices and strip centroids."""
    count=len(profiles);centres=np.array([vertices[p].mean(0)for p in profiles])
    owner=np.full(len(vertices),-1,int)
    for i,p in enumerate(profiles):owner[p]=i
    if whole_source:
        influence=profile_translation_weights(vertices,faces,profiles,physical_ids=physical_ids)
        points=[vertices,vertices[faces].mean(1)];weights=[influence,influence[faces].mean(1)]
        for a,b in [(0,1),(1,2),(2,0)]:
            points.append((vertices[faces[:,a]]+vertices[faces[:,b]])/2)
            weights.append((influence[faces[:,a]]+influence[faces[:,b]])/2)
        strip=faces
    else:
        sample_ids=np.concatenate(profiles);points=[vertices[sample_ids]]
        weights=[np.eye(count)[owner[sample_ids]]]
        eligible=np.all(owner[faces]>=0,axis=1)
        strip=faces[eligible];points.append(vertices[strip].mean(1))
        weights.append(np.eye(count)[owner[strip]].mean(1))
    points=np.concatenate(points);W=np.concatenate(weights)
    L=np.zeros((count,count))
    edge=np.zeros((len(links),count))
    for k,(a,b)in enumerate(links):
        L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
        edge[k,a]=-1;edge[k,b]=1
    K=np.linalg.inv(np.eye(count)+2*L.T@L)
    K=(K+K.T)*.5
    direction=centres[links[:,1]]-centres[links[:,0]] if len(links)else np.empty((0,3))
    length=np.linalg.norm(direction,axis=1)
    if np.any(length<=1e-10):raise ValueError('Coincident linked section centres')
    normal=direction/length[:,None];delta=np.zeros((count,3));runs=[]
    prior_w=[];prior_n=[];prior_b=[]
    for iteration in range(2):
        p=points+W@delta;near,_=bone.closest_points(p)
        inside=bone.scene.compute_occupancy(library().core.Tensor(((p-bone.origin)/bone.scale).astype(np.float32)),nthreads=2,nsamples=3).numpy()>.5
        d=p-near['points'];d[inside]*=-1;dn=np.linalg.norm(d,axis=1)
        if np.any(dn<1e-12):raise ValueError('Ambiguous bone witness')
        bn=d/dn[:,None];signed=near['distance_m']*np.where(inside,-1.,1.)
        sn,_=skin.closest_points(p);t=skin.mesh.triangles[sn['face_ids']]
        outward=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);outward/=np.linalg.norm(outward,axis=1)[:,None]
        inward=-outward;inside_skin=np.einsum('ij,ij->i',p-sn['points'],inward)
        prior_w.extend([W,W]);prior_n.extend([bn,inward])
        prior_b.extend([.0002-signed+np.einsum('ij,ij->i',bn,W@delta),
                        .0002-inside_skin+np.einsum('ij,ij->i',inward,W@delta)])
        w=np.concatenate([edge,*prior_w]);n=np.concatenate([normal,*prior_n]);rhs=np.concatenate([-.4*length,*prior_b])
        proposed,evidence=solve_weighted_section_socp(K,w,n,rhs,maximum_motion_m=trust_m,
            residual_m=1e-8,time_limit_seconds=.4,maximum_iterations=60)
        runs.append(finite_diagnostic(evidence))
        if proposed is None:
            return centres,dict(complete=False,runs=runs,actualSourceStripCentroidSamples=len(strip),
                                status='ConstraintProposalFailed',sourceRetainedLiteral=True)
        delta=proposed
    return centres+delta,dict(complete=True,runs=runs,status='Proposed',
        actualSourceStripCentroidSamples=len(strip),maximumSectionTranslation_m=float(np.linalg.norm(delta,axis=1).max()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['candidate','profiles','binding','bones','skin-reference','output']:parser.add_argument('--'+name,required=True)
    parser.add_argument('--whole-source-constraints',action='store_true')
    parser.add_argument('--split-components',action='store_true')
    parser.add_argument('--align-section-frames',action='store_true')
    parser.add_argument('--core-method',choices=['harmonic','arap'],default='harmonic')
    parser.add_argument('--maximum-section-motion-mm',type=float,default=2.)
    args=parser.parse_args();start=time.monotonic()
    if not 0<args.maximum_section_motion_mm<=5:parser.error('Section motion must be within (0,5] mm')
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in ['candidate','profiles','binding','bones','skin-reference']}
    data,bank,binding,bones,sr=[json.loads(p.read_text())for p in paths.values()]
    if binding['sourceBankSha256']!=sha(paths['profiles']):raise ValueError('Profile ownership declaration mismatch')
    if not sr['defaultHiddenExcluded']or sha(sr['path'])!=sr['sha256']:raise ValueError('Visible skin reference changed')
    side=data['side'];suffix='.l'if side=='left'else'.r';fixed=[]
    for r in bones['rows']:
        if r['sourceObject'].endswith(suffix):
            if sha(r['path'])!=r['sha256']:raise ValueError('Fixed bone changed')
            fixed.append(Mesh.load(r['path']))
    t=np.concatenate([m.triangles for m in fixed]);bone=SurfaceIndex(Mesh(t.reshape(-1,3),np.arange(t.size//3).reshape(-1,3),fixed[0].frame))
    skin=SurfaceIndex(Mesh.load(sr['path']));out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new candidate stage')
    out.mkdir();rows=[]
    for row in data['rows']:
        if sha(row['path'])!=row['sha256']:raise ValueError('Input pose changed')
        mesh=Mesh.load(row['path']);source=next(r for r in bank['rows']if r['partId']==row['partId'])
        if sha(source['path'])!=source['sha256']or source['sourceSha256']!=row['sourceSha256']:raise ValueError('Source profile bank changed')
        owners=[r for r in binding['rows']if r['side']==side and r['partId']==row['partId']]
        if not owners or not all(r['anatomicallyReviewed']for r in owners):raise ValueError('Unreviewed anatomical component')
        with np.load(source['path'],allow_pickle=False)as z:
            raw=z['profiles'].copy();physical=z['physical_ids'].copy();links=z['profile_links'].copy()
        kept=[];mapping={};owned=set()
        for i,p in enumerate(raw):
            identities=set(map(int,physical[p]))
            if len(identities)!=len(p)or owned.intersection(identities):continue
            if sum(set(map(int,p))<=set(r['rawVertexIds'])for r in owners)!=1:raise ValueError('Profile crosses declared anatomical components')
            mapping[i]=len(kept);kept.append(p);owned.update(identities)
        links=np.array([[mapping[a],mapping[b]]for a,b in links if a in mapping and b in mapping],int).reshape(-1,2)
        v=mesh.vertices.copy();gauge=0.;proofs=[]
        groups=owners if args.split_components else [dict(rawVertexIds=list(range(len(v))),sourceFaceIds=list(range(len(mesh.faces))),physicalComponent='whole')]
        for group in groups:
            ids=np.array(group['rawVertexIds'],int);fids=np.array(group['sourceFaceIds'],int)
            remap=np.full(len(v),-1,int);remap[ids]=np.arange(len(ids));local_faces=remap[mesh.faces[fids]]
            if np.any(local_faces<0):raise ValueError('Component does not own its faces')
            selected=[k for k,p in enumerate(kept)if np.all(remap[p]>=0)];positions={k:i for i,k in enumerate(selected)}
            local_profiles=[remap[kept[k]]for k in selected]
            local_links=np.array([[positions[a],positions[b]]for a,b in links if a in positions and b in positions],int).reshape(-1,2)
            if not selected:raise ValueError('No source profiles in reviewed component')
            target,part_proof=section_proposal(mesh.vertices[ids],local_faces,local_profiles,local_links,bone,skin,
                physical_ids=physical[ids],whole_source=args.whole_source_constraints,trust_m=args.maximum_section_motion_mm*.001)
            part_proof['sourceComponent']=group['physicalComponent'];proofs.append(part_proof)
            if part_proof['complete']:
                rotations=np.repeat(np.eye(3)[None],len(local_profiles),axis=0)
                if args.align_section_frames:
                    centres=np.array([mesh.vertices[ids][p].mean(0)for p in local_profiles])
                    rotations,frames=transport_linked_profiles(centres,target,rotations,local_links)
                    part_proof['sectionFrameTransport']=frames
                moved=transport_profile_cores(mesh.vertices[ids],local_faces,local_profiles,target,
                    rotations,physical_ids=physical[ids],core_method=args.core_method,
                    core_maximum_steps=32,core_wall_seconds=.3)
                part_proof['coreTransport']=dict(method=args.core_method,arap=moved['arap_evidence'])
                v[ids]=moved['vertices']
            if time.monotonic()-start>22:raise TimeoutError('22-second native section clearance stage')
        proof=dict(complete=all(p['complete']for p in proofs),components=proofs,
            wholeSourceVertexEdgeAndFaceSamples=args.whole_source_constraints,maximumSectionMotionBound_m=args.maximum_section_motion_mm*.001)
        for p in kept:
            gauge=max(gauge,float(np.max(abs(np.linalg.norm(v[p][:,None]-v[p][None,:],axis=2)-np.linalg.norm(mesh.vertices[p][:,None]-mesh.vertices[p][None,:],axis=2)))))
        if gauge>1e-12:raise ValueError('Section gauge changed')
        face_n=np.cross(v[mesh.faces[:,1]]-v[mesh.faces[:,0]],v[mesh.faces[:,2]]-v[mesh.faces[:,0]])
        _,inv=np.unique(physical,return_inverse=True);normals=np.zeros((inv.max()+1,3));np.add.at(normals,inv[mesh.faces.ravel()],np.repeat(face_n,3,axis=0));normals=normals[inv]
        length=np.linalg.norm(normals,axis=1);old=mesh.attributes.get('normals',np.zeros_like(normals)).astype(float);normals=np.divide(normals,length[:,None],out=old.copy(),where=length[:,None]>0)
        attrs=dict(mesh.attributes,normals=normals,normalSigns=np.ones(len(v)))
        candidate=Mesh(v,mesh.faces,mesh.frame,attrs,dict(sourcePartId=row['partId']));path=out/f'part-{len(rows):02}.npz';candidate.save(path)
        rows.append(dict(row,path=str(path),sha256=sha(path),geometrySha256=candidate.geometry_sha,
            parentPosePath=row['path'],parentPoseSha256=row['sha256'],sourceComponentShapesPreserved=False,
            rigidProfiles=len(kept),gaugeMaximumError_m=gauge,clearanceProposal=proof))
        if time.monotonic()-start>22:raise TimeoutError('22-second native section clearance stage')
    report=dict(schema='local-hand-section-clearance/v1',side=side,rows=rows,
        inputs={k:dict(path=str(p),sha256=sha(p))for k,p in paths.items()},sourceCodeSha256=sha(__file__),
        skeletonChanged=False,sourceFacesRetained=True,sectionGaugeRetained=True,skinContainmentEstablished=False,
        installed=False,installationAllowed=False,elapsedSeconds=time.monotonic()-start,
        complete=all(r['clearanceProposal']['complete']for r in rows))
    write_json(out/'report.json',report);print(json.dumps(dict(parts=len(rows),proposed=sum(r['clearanceProposal']['complete']for r in rows),elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
