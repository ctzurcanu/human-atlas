"""Couple a reviewed parent/digit pose through the complete source wrist graph.

All source section closures and real shared mouth controls move together.
Fixed forearm controls and the explicit parent/digit targets remain literal.
No anatomical acceptance or installation is performed by this constructor.
"""
import argparse
import copy
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
from geometry_tools.commands.couple_hand_arteries import rotation


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','candidate','profiles','fixed-parts','output']:
        p.add_argument('--'+key,required=True)
    p.add_argument('--side',choices=['left','right'],required=True)
    p.add_argument('--fixed-proximal-y',type=float,required=True)
    p.add_argument('--fixed-component-binding',help='Explicit reviewed source component inventory')
    p.add_argument('--fixed-component-part')
    p.add_argument('--fixed-component',type=int)
    p.add_argument('--arap-initialization',action='store_true',
                   help='Preserve branch core shape during the common wrist initialization')
    p.add_argument('--common-spatial-flow',action='store_true',
                   help='Use one bounded smooth spatial flow through the explicit passed source targets')
    args=p.parse_args(); started=time.monotonic()
    if args.common_spatial_flow and args.arap_initialization:raise ValueError('Choose one explicit initialization')
    paths={k:Path(getattr(args,k)).resolve()for k in ['source','candidate','profiles']}
    data={k:json.loads(v.read_text())for k,v in paths.items()}
    component_owner=None
    if args.fixed_component_binding:
        paths['component-binding']=Path(args.fixed_component_binding).resolve()
        owners=json.loads(paths['component-binding'].read_text())['rows']
        owners=[r for r in owners if r['partId']==args.fixed_component_part and r['physicalComponent']==args.fixed_component]
        if len(owners)!=1 or not owners[0]['anatomicallyReviewed']:
            raise ValueError('Explicit reviewed component ownership required')
        component_owner=owners[0]
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new cohort stage')
    banks={r['partId']:r for r in data['profiles']['rows']}
    target={r['partId']:r for r in data['candidate']['rows']}
    rows=[r for r in data['source']['rows']if r['side']==args.side]
    fixed=set(args.fixed_parts.split(','))
    if not fixed.issubset({r['partId']for r in rows}):raise ValueError('Fixed part absent')
    meshes={};offsets={};offset=0;faces=[];raw=[];groups=[]
    for row in rows:
        if sha(row['path'])!=row['sha256']:raise ValueError('Source changed')
        m=Mesh.load(row['path']);meshes[row['partId']]=m;offsets[row['partId']]=offset
        graph=source_profile_graph(m.vertices,m.faces,
            triangle_pairing=load_source_pairing(banks.get(row['partId']),m.faces),extend_annuli=True)
        groups.extend(np.array(s)+offset for s in graph['profiles']+graph['boundary_raw'])
        faces.append(m.faces+offset);raw.append(m.vertices);offset+=len(m.vertices)
    raw=np.concatenate(raw);faces=np.concatenate(faces)
    points,representative,physical=np.unique(raw,axis=0,return_index=True,return_inverse=True)
    merged=[]
    for ids in groups:
        members=set(map(int,physical[ids]));touch=[i for i,g in enumerate(merged)if g.intersection(members)]
        for i in reversed(touch):members.update(merged.pop(i))
        merged.append(members)
    sections=[representative[np.array(sorted(g))]for g in merged]
    prescribed={int(k):points[k].copy()for k in np.flatnonzero(points[:,1]>=args.fixed_proximal_y)}
    fixed_requests=[(pid,None)for pid in fixed]
    if component_owner:
        fixed_requests.append((args.fixed_component_part,np.array(component_owner['sourceFaceIds'],int)))
    for pid,face_ids in fixed_requests:
        row=target[pid]
        if sha(row['path'])!=row['sha256']:raise ValueError('Target changed')
        m=Mesh.load(row['path']);old=meshes[pid]
        if not np.array_equal(m.faces,old.faces):raise ValueError('Target source faces changed')
        selected_raw=np.arange(len(m.vertices))if face_ids is None else np.unique(old.faces[face_ids])
        ids=physical[offsets[pid]+selected_raw]
        for k,q in zip(ids,m.vertices[selected_raw]):
            if int(k)in prescribed and np.linalg.norm(prescribed[int(k)]-q)>1e-12:
                raise ValueError('Fixed targets disagree at an actual shared source control')
            prescribed[int(k)]=q.copy()
    keys=np.array(sorted(prescribed));controls=representative[keys];values=np.array([prescribed[int(k)]for k in keys])
    flow_proof=None
    if args.common_spatial_flow:
        from geometry_tools.algorithms.landmark_flow import landmark_flow
        moving=[];targets=[];stationary=[]
        for s in sections:
            if all(int(k)in prescribed for k in physical[s]):
                old=raw[s].mean(0);new=np.array([prescribed[int(k)]for k in physical[s]]).mean(0)
                if np.linalg.norm(new-old)>1e-7:moving.append(old);targets.append(new)
                elif old[1]<args.fixed_proximal_y:stationary.append(old)
        if not moving:raise ValueError('No explicit moved source section targets')
        anchors=np.array([[x,y,z]for x in [raw[:,0].min()-.02,raw[:,0].max()+.02]
                          for y in [args.fixed_proximal_y+.01,args.fixed_proximal_y+.05]
                          for z in [raw[:,2].min()-.02,raw[:,2].max()+.02]])
        if stationary:
            # Bounded farthest-point coverage of explicitly fixed source
            # sections prevents a local wrist pose from dragging other digits.
            stationary=np.array(stationary);keep=[int(np.argmin(stationary[:,1]))]
            distance=np.linalg.norm(stationary-stationary[keep[0]],axis=1)
            while len(keep)<min(32,len(stationary)):
                k=int(np.argmax(distance))
                if distance[k]<1e-10:break
                keep.append(k);distance=np.minimum(distance,np.linalg.norm(stationary-stationary[k],axis=1))
            anchors=np.r_[anchors,stationary[keep]]
        old=np.r_[np.array(moving),anchors];new=np.r_[np.array(targets),anchors]
        desired,_,flow_proof=landmark_flow(raw,old,new,wall_seconds=7.,
            fixed_above=args.fixed_proximal_y,full_below=args.fixed_proximal_y-.04)
        first=dict(vertices=desired)
    else:
        first=transport_profile_cores(raw,faces,[],np.empty((0,3)),np.empty((0,3,3)),
            physical_ids=physical,prescribed_vertices=controls,prescribed_positions=values,
            core_method='arap'if args.arap_initialization else'harmonic',core_wall_seconds=3.)
    desired=first['vertices']; centres=[];rotations=[]
    for s in sections:
        local=np.array([i for i,k in enumerate(physical[s])if int(k)in prescribed])
        if len(local)>=3:
            source=raw[s[local]];posed=np.array([prescribed[int(k)]for k in physical[s[local]]])
            R=rotation(source,posed);centre=posed.mean(0)+(raw[s].mean(0)-source.mean(0))@R.T
        else:
            R=rotation(raw[s],desired[s]);centre=desired[s].mean(0)
            if len(local):
                anchors=np.array([prescribed[int(k)]for k in physical[s[local]]])
                centre+=(anchors-((raw[s[local]]-raw[s].mean(0))@R.T+centre)).mean(0)
        q=(raw[s]-raw[s].mean(0))@R.T+centre
        if len(local)and np.linalg.norm(q[local]-np.array([prescribed[int(k)]for k in physical[s[local]]]),axis=1).max()>1e-10:
            raise ValueError('Fixed target conflicts with a complete rigid shared source section')
        centres.append(centre);rotations.append(R)
    if args.common_spatial_flow:
        # Preserve the common spatially posed branch cores, correcting only
        # residual section restoration, rather than re-solving the complete
        # displacement harmonically on a different material graph.
        boundary={int(k):v.copy()for k,v in zip(keys,values)}
        for s,c,R in zip(sections,centres,rotations):
            pose=(raw[s]-raw[s].mean(0))@R.T+c
            for k,v in zip(physical[s],pose):
                if int(k)in boundary and np.linalg.norm(boundary[int(k)]-v)>1e-10:
                    raise ValueError('Shared spatial section conflicts with fixed target')
                boundary[int(k)]=v
        ids=np.array(sorted(boundary));posed=np.array([boundary[int(k)]for k in ids])
        second=transport_profile_cores(desired,faces,[],np.empty((0,3)),np.empty((0,3,3)),
            physical_ids=physical,prescribed_vertices=representative[ids],prescribed_positions=posed)
    else:
        second=transport_profile_cores(raw,faces,sections,np.array(centres),np.array(rotations),
            physical_ids=physical,prescribed_vertices=controls,prescribed_positions=values)
    q=second['vertices']
    for k,value in zip(keys,values):q[physical==k]=value
    gauge=max(float(np.abs(np.linalg.norm(q[s][:,None]-q[s][None,:],axis=2)-
                            np.linalg.norm(raw[s][:,None]-raw[s][None,:],axis=2)).max())for s in sections)
    if gauge>1e-10:raise ValueError('Complete source section gauge changed')
    if np.max(np.linalg.norm(q[controls]-values,axis=1))>1e-12:raise ValueError('Fixed controls moved')
    # One physical normal field across actual joined source mouths.
    n=np.zeros_like(points);t=q[faces];fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    np.add.at(n,physical[faces.ravel()],np.repeat(fn,3,axis=0));n=n[physical]
    n/=np.linalg.norm(n,axis=1)[:,None]
    out.mkdir();result=copy.deepcopy(data['candidate']);emitted=[]
    for row in result['rows']:
        pid=row['partId']
        if pid not in meshes:continue
        old=meshes[pid];a=offsets[pid];sl=slice(a,a+len(old.vertices))
        attrs=dict(old.attributes,normals=np.rint(n[sl]*32767).astype('<i2'),normalSigns=np.ones(len(old.vertices),'i1'))
        m=Mesh(q[sl],old.faces,old.frame,attrs,old.source);file=out/f'part-{len(emitted):02}.npz';m.save(file)
        row.update(path=str(file),sha256=sha(file),geometrySha256=m.geometry_sha,
                   sourcePath=next(r['path']for r in rows if r['partId']==pid),
                   sourceSha256=next(r['sha256']for r in rows if r['partId']==pid),
                   sourceTrianglePairing=dict(path=banks[pid]['path'],sha256=banks[pid]['sha256']))
        emitted.append(dict(partId=pid,maximumMotion_m=float(np.linalg.norm(q[sl]-raw[sl],axis=1).max())))
    if time.monotonic()-started>22:raise TimeoutError('22-second complete cohort stage')
    construction_complete=not args.arap_initialization or bool(first['arap_evidence']['converged'])
    result.update(complete=False,constructionComplete=construction_complete,installed=False,installationAllowed=False,
        skeletonChanged=False,portsRequireRecoupling=False,
        wristCohort=dict(side=args.side,parts=emitted,fixedParts=sorted(fixed),rigidControlClosures=len(sections),
                        fixedComponent=None if component_owner is None else dict(partId=args.fixed_component_part,
                            component=args.fixed_component,sourceFaceIds=component_owner['sourceFaceIds']),
                        maximumGaugeError_m=gauge,actualSharedSourceControls=True,
                        initialization=flow_proof if args.common_spatial_flow else first.get('arap_evidence'),
                        commonSpatialFlow=args.common_spatial_flow,
                        fixedProximalY_m=args.fixed_proximal_y),
        pending=['Complete source-network crossings, bone/skin/tissue checks, packed seam and Blender review'],
        inputs={k:dict(path=str(v),sha256=sha(v))for k,v in paths.items()},
        sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json',result)
    print(json.dumps(dict(parts=len(emitted),rigidControlClosures=len(sections),gaugeError_m=gauge,
                          elapsedSeconds=result['elapsedSeconds'])))


if __name__=='__main__':main()
