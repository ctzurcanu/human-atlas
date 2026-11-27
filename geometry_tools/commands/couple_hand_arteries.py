"""Explicit source-port coupling for the bilateral hand arterial assemblies.

Only the hash-bound, uniquely paired, oppositely wound source mouths may share
controls. A source-only seam pass can anchor already defective source parts;
the fitting mode uses distal hand targets and a common palm/wrist carrier.
Both modes are diagnostic until independent geometry/storage/render checks.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
from geometry_tools.algorithms.source_tube_topology import recover_regular_source_annuli
from geometry_tools.algorithms.tube_continuity import boundary_loops


def rotation(a,b):
    u,_,vt=np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)));R=vt.T@u.T
    if np.linalg.det(R)<0:vt[-1]*=-1;R=vt.T@u.T
    return R


def rank(pid):
    if 'proper_'in pid or 'dorsal_digital' in pid:return 5
    if any(s in pid for s in ['common_palmar','princeps','radialis']):return 4
    if 'metacarpal'in pid:return 3
    if 'arch'in pid:return 2
    return 1


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','ports','pairs','profiles','left','right','left-pose','right-pose','source-self','output']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--source-only',action='store_true');args=p.parse_args();start=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in ['source','ports','pairs','profiles','left','right','left-pose','right-pose','source-self']}
    data={k:json.loads(q.read_text())for k,q in paths.items()}
    if data['pairs']['sourceInventorySha256']!=sha(paths['ports'])or data['pairs']['ambiguousMouths']:raise ValueError('Port declaration changed or ambiguous')
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new coupled stage')
    out.mkdir();native={};base={};profiles={};physical={};source_rows={r['partId']:r for r in data['source']['rows']}
    defect={r['partId']for r in data['source-self']['rows']if r['pairs']}
    seed={r['partId']:r for side in ['left','right']for r in data[side]['rows']}
    palms={}
    for side in ['left','right']:
        q=next(r for r in data[side+'-pose']['rows']if 'common_palmar_digital_arteries'in r['partId'])
        palms[side]=np.array(q['componentTransforms'][0]['matrix'])
    for pid,row in source_rows.items():
        if sha(row['path'])!=row['sha256']:raise ValueError('Source changed')
        mesh=Mesh.load(row['path']);native[pid]=mesh
        if args.source_only:q=mesh.vertices.copy()
        elif pid in seed:
            r=seed[pid]
            if sha(r['path'])!=r['sha256']:raise ValueError('Distal target changed')
            q=Mesh.load(r['path']).vertices.copy()
        else:
            M=palms[row['side']];target=mesh.vertices@M[:3,:3].T+M[:3,3]
            t=np.clip((1.07-mesh.vertices[:,1])/.13,0,1);fade=t*t*t*(10-15*t+6*t*t)
            q=mesh.vertices+fade[:,None]*(target-mesh.vertices)
        base[pid]=q
        old=next((r for r in data['profiles']['rows']if r['partId']==pid),None)
        if old:
            if sha(old['path'])!=old['sha256']:raise ValueError('Profile bank changed')
            with np.load(old['path'],allow_pickle=False)as z:
                if not np.array_equal(z['source_vertices'],mesh.vertices)or not np.array_equal(z['source_faces'],mesh.faces):raise ValueError('Profile arrays do not replay current whole source')
                profiles[pid]=z['profiles'].copy();physical[pid]=z['physical_ids'].copy()
        else:
            paired=np.arange(len(mesh.faces)//2*2).reshape(-1,2)
            recovered=recover_regular_source_annuli(mesh.vertices,mesh.faces,paired,ring_size=6)
            inv=recovered['raw_to_physical'];rep={int(k):int(np.flatnonzero(inv==k)[0])for k in np.unique(inv)}
            profiles[pid]=[np.array([rep[k]for k in ring],int)for ring in recovered['rings']];physical[pid]=inv
        if time.monotonic()-start>22:raise TimeoutError('22-second coupled hand construction')
    prescribed={pid:{}for pid in native};mouths=data['ports']['mouths'];pairs=[]
    # Every unmatched boundary retains its own carrier position; these are not
    # silently classified as newly constructed physiological openings.
    for pid,m in native.items():
        _,rep,inv=np.unique(physical[pid],return_index=True,return_inverse=True)
        for loop in boundary_loops(inv[m.faces]):
            for raw in rep[loop]:prescribed[pid][int(physical[pid][raw])]=base[pid][raw].copy()
    for link in data['pairs']['rows']:
        a,b=[mouths[link[k]]for k in ['first','second']];pa,pb=a['partId'],b['partId'];ia=np.array(a['rawVertices']);ib=np.array(b['rawVertices'])[link['permutation']]
        if not np.array_equal(native[pa].vertices[ia],a['points'])or not np.array_equal(native[pb].vertices[np.array(b['rawVertices'])],b['points']):raise ValueError('Source mouth indexing changed')
        perm=np.array(link['permutation'])
        if not np.all((np.roll(perm,-1)-perm)%len(perm)==len(perm)-1):raise ValueError('Source mouths are not oppositely wound')
        if args.source_only and pa in defect:target=base[pa][ia].copy();canonical=pa
        elif args.source_only and pb in defect:target=base[pb][ib].copy();canonical=pb
        elif rank(pa)>=rank(pb):target=base[pa][ia].copy();canonical=pa
        else:target=base[pb][ib].copy();canonical=pb
        for pid,ids in [(pa,ia),(pb,ib)]:
            for raw,point in zip(ids,target):prescribed[pid][int(physical[pid][raw])]=point.copy()
        pairs.append(dict(firstPart=pa,secondPart=pb,firstRaw=ia.tolist(),secondRaw=ib.tolist(),canonicalPart=canonical,sourceMaximumGap_m=link['maximumGap_m']))
    result={};rows=[]
    for pid,m in native.items():
        ph=physical[pid];rep={int(k):int(np.flatnonzero(ph==k)[0])for k in np.unique(ph)}
        boundaries=sorted(prescribed[pid]);raw=np.array([rep[k]for k in boundaries],int);targets=np.array([prescribed[pid][k]for k in boundaries])
        if args.source_only and pid in defect:
            q=m.vertices.copy();kept=[];gauge=0.
            if not np.array_equal(q[raw],targets):raise ValueError('Source defect anchor would move')
        else:
            first=transport_profile_cores(base[pid],m.faces,[],np.empty((0,3)),np.empty((0,3,3)),physical_ids=ph,prescribed_vertices=raw,prescribed_positions=targets)
            kept=[];seen=set();fixed=set(boundaries)
            for section in profiles[pid]:
                own=set(map(int,ph[section]))
                if len(own)!=len(section)or own&seen or own&fixed:continue
                seen|=own;kept.append(section)
            centres=np.array([first['vertices'][s].mean(0)for s in kept]).reshape(-1,3)
            rotations=np.array([rotation(base[pid][s],first['vertices'][s])for s in kept]).reshape(-1,3,3)
            moved=transport_profile_cores(base[pid],m.faces,kept,centres,rotations,physical_ids=ph,prescribed_vertices=raw,prescribed_positions=targets)
            q=moved['vertices'];gauge=0.
            for section in kept:
                gauge=max(gauge,float(np.max(abs(np.linalg.norm(q[section][:,None]-q[section][None,:],axis=2)-np.linalg.norm(base[pid][section][:,None]-base[pid][section][None,:],axis=2)))))
            if gauge>1e-12:raise ValueError('Interior section gauge changed')
        attrs={k:v.copy()for k,v in m.attributes.items()}
        if not args.source_only:
            normal=np.zeros((len(rep),3));_,inv=np.unique(ph,return_inverse=True);t=q[m.faces];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);np.add.at(normal,inv[m.faces.ravel()],np.repeat(n,3,axis=0));normal=normal[inv];length=np.linalg.norm(normal,axis=1)
            old=attrs['normals'].astype(float)*attrs.get('normalSigns',np.ones(len(q)))[:,None];old/=np.linalg.norm(old,axis=1)[:,None]
            normal=np.divide(normal,length[:,None],out=old.copy(),where=length[:,None]>0);attrs['normals']=np.rint(normal*32767).astype('<i2');attrs['normalSigns']=np.ones(len(q),'i1')
        result[pid]=Mesh(q,m.faces,m.frame,attrs,dict(sourcePartId=pid))
        rows.append(dict(source_rows[pid],sourcePath=source_rows[pid]['path'],sourceSha256=source_rows[pid]['sha256'],retainedRigidProfiles=len(kept),gaugeError_m=gauge,sourceDefectAnchoredLiteral=args.source_only and pid in defect))
        if time.monotonic()-start>22:raise TimeoutError('22-second coupled hand construction')
    final_attributes={pid:{k:a.copy()for k,a in mesh.attributes.items()}for pid,mesh in result.items()}
    for link in pairs:
        pa,pb=link['firstPart'],link['secondPart'];a=np.array(link['firstRaw']);b=np.array(link['secondRaw']);ma,mb=result[pa],result[pb]
        if not np.array_equal(ma.vertices[a],mb.vertices[b]):raise ValueError('Shared port coordinates failed literal replay')
        canonical=link['canonicalPart'];mc=result[canonical];ids=a if canonical==pa else b;ca=final_attributes[canonical]
        shared=ca['normals'][ids].astype(np.int32)*ca.get('normalSigns',np.ones(len(mc.vertices),int))[ids,None]
        for pid,raw in [(pa,a),(pb,b)]:
            mesh=result[pid]
            if args.source_only and pid in defect:continue
            for k,word in zip(physical[pid][raw],shared):
                mask=physical[pid]==k;final_attributes[pid]['normals'][mask]=word;final_attributes[pid]['normalSigns'][mask]=1
        link['literalCoordinatesAndEffectiveNormals']=True
    for i,row in enumerate(rows):
        old=result[row['partId']];m=Mesh(old.vertices,old.faces,old.frame,final_attributes[row['partId']],old.source);path=out/f'part-{i:02}.npz';m.save(path);row.update(path=str(path),sha256=sha(path),geometrySha256=m.geometry_sha,maximumSourceMotion_m=float(np.linalg.norm(m.vertices-native[row['partId']].vertices,axis=1).max()))
    report=dict(schema='source-bound-coupled-hand-arteries/v1',rows=rows,pairs=pairs,sourceOnly=args.source_only,
        inputs={k:dict(path=str(q),sha256=sha(q))for k,q in paths.items()},sourceCodeSha256=sha(__file__),
        skeletonChanged=False,sourceFacesLiteral=True,newFaces=0,complete=True,installed=False,installationAllowed=False,elapsedSeconds=time.monotonic()-start)
    write_json(out/'report.json',report);print(json.dumps(dict(parts=len(rows),pairedPorts=len(pairs),sourceOnly=args.source_only,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
