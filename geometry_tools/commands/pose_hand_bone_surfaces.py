"""Carry the complete detailed arterial hand by paired bone surface frames.

Source and recipient bones are immutable guides. Measured width/length ratios
place section centres; every actual vessel section keeps its original gauge.
This proposal requires independent full tissue, crossing and render checks.
"""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha,write_json
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.algorithms.paired_bone_transport import surface_frame_carriage
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.commands.check_hand_network import combine


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['source','profiles','bone-reference','carpal-reference','output']:
        parser.add_argument('--'+key,required=True)
    parser.add_argument('--maximum-carrier-scale',type=float,default=2.5)
    args=parser.parse_args();started=time.monotonic()
    paths={k:Path(getattr(args,k.replace('-','_'))).resolve()for k in
           ['source','profiles','bone-reference','carpal-reference']}
    data={k:json.loads(p.read_text())for k,p in paths.items()}
    out=Path(args.output).resolve()
    if out.exists():raise ValueError('Use a new bounded source carriage stage')
    out.mkdir();banks={r['partId']:r for r in data['profiles']['rows']}
    refs=[r for r in data['bone-reference']['rows']if r['side']=='left']
    carpal=data['carpal-reference']['rows'][-1]['namedBoneCorrespondences']
    def points(row,source):
        if source:
            path=row.get('sourcePath',row.get('sourceMesh',{}).get('path'))
            digest=row.get('sha256',row.get('sourceMesh',{}).get('sha256'))
            if sha(path)!=digest:raise ValueError('Donor bone guide changed')
            return Mesh.load(path)
        pieces=row['targetPieces']
        for r in pieces:
            if sha(r['path'])!=r['sha256']:raise ValueError('Fixed recipient bone changed')
        return combine([Mesh.load(r['path'])for r in pieces])
    frames=[]
    for source in [True,False]:
        centres=[];directions=[]
        for r in refs:
            if r['segment']!='metacarpal':continue
            v=np.unique(points(r,source).vertices,axis=0);centre=v.mean(0)
            _,_,vt=np.linalg.svd(v-centre,full_matrices=False);axis=vt[0]
            if axis[1]>0:axis=-axis
            centres.append((r['digit'],centre));directions.append(axis)
        z=np.mean(directions,axis=0);z/=np.linalg.norm(z)
        centres=dict(centres);x=centres[1]-centres[5];x-=z*(x@z);x/=np.linalg.norm(x)
        frames.append(np.column_stack([x,np.cross(z,x),z]))
    carriers=[];evidence=[]
    for r in refs+carpal:
        meshes=[points(r,s)for s in [True,False]];local=[]
        for m,palm in zip(meshes,frames):
            v=np.unique(m.vertices,axis=0);centre=v.mean(0);frame=palm.copy()
            if r in refs:
                _,_,vt=np.linalg.svd(v-centre,full_matrices=False);z=vt[0]
                if z@palm[:,2]<0:z=-z
                y=palm[:,1]-z*(z@palm[:,1]);y/=np.linalg.norm(y)
                frame=np.column_stack([np.cross(y,z),y,z])
            width=np.diff(np.quantile((v-centre)@frame,[.05,.95],axis=0),axis=0)[0]
            local.append((centre,frame,width))
        a,z=local;scale=z[2]/a[2]
        if np.any(scale<.5)or np.any(scale>args.maximum_carrier_scale):
            write_json(out/'failure.json',dict(complete=False,installed=False,geometryEmitted=False,
                bone=r.get('targetOwner',r['segment']),measuredScales=scale.tolist(),
                sourceWidths_m=a[2].tolist(),targetWidths_m=z[2].tolist(),
                maximumCarrierScale=args.maximum_carrier_scale,
                reason='Measured width ratio outside the explicit guide bound'))
            raise ValueError('Measured bone width ratio exceeds the explicit carrier bound: '+str(r.get('targetOwner',r['segment']))+' '+str(scale.tolist()))
        carriers.append(dict(sourceFrame=a[1],targetFrame=z[1],sourceCentre=a[0],targetCentre=z[0],
                             scales=scale,index=SurfaceIndex(meshes[0])))
        evidence.append(dict(name=r.get('targetOwner',r['segment']),sourceBoneSha256=meshes[0].geometry_sha,
            targetBoneSha256=meshes[1].geometry_sha,scales=scale.tolist(),sourceFrame=a[1].tolist(),targetFrame=z[1].tolist()))
    selected=[r for r in data['source']['rows']if r['side']=='left']
    meshes=[];raw=[];sections=[];offset=0
    for r in selected:
        if sha(r['path'])!=r['sha256']:raise ValueError('Source vessel changed')
        m=Mesh.load(r['path']);meshes.append(m);raw.append(m.vertices)
        g=source_profile_graph(m.vertices,m.faces,triangle_pairing=load_source_pairing(banks[r['partId']],m.faces),extend_annuli=True)
        sections.extend(np.array(s)+offset for s in g['profiles']+g['boundary_raw']);offset+=len(m.vertices)
    raw=np.concatenate(raw);_,rep,physical=np.unique(raw,axis=0,return_index=True,return_inverse=True)
    closures=[]
    for section in sections:
        group=set(map(int,physical[section]));touch=[i for i,g in enumerate(closures)if g.intersection(group)]
        for i in reversed(touch):group.update(closures.pop(i))
        closures.append(group)
    sections=[rep[np.array(sorted(g))]for g in closures]
    q,R,proof=surface_frame_carriage(raw,sections,carriers,physical,
                                   maximum_carrier_scale=args.maximum_carrier_scale)
    fade=np.clip((1.08-raw[:,1])/.12,0,1);fade=fade**3*(10-15*fade+6*fade**2)
    q=raw+fade[:,None]*(q-raw)
    for section in sections:
        centre=raw[section].mean(0);w=np.clip((1.08-centre[1])/.12,0,1);w=w**3*(10-15*w+6*w*w)
        rotation=Rotation.from_rotvec(Rotation.from_matrix(R[section[0]]).as_rotvec()*w).as_matrix()
        placed=centre+(q[section].mean(0)-centre)
        for i in section:
            aliases=physical==physical[i];q[aliases]=(raw[aliases]-centre)@rotation.T+placed
    # Alias and gauge replay are independent of the shape-placement proof.
    if not np.array_equal(q,q[rep][physical]):raise ValueError('Shared source mouth controls separated')
    gauge=max(float(np.abs(np.linalg.norm(q[s][:,None]-q[s][None,:],axis=2)-np.linalg.norm(raw[s][:,None]-raw[s][None,:],axis=2)).max())for s in sections)
    if gauge>1e-10:raise ValueError('Source section gauge changed')
    faces=np.concatenate([m.faces+sum(len(x.vertices)for x in meshes[:i])for i,m in enumerate(meshes)])
    n=np.zeros((len(rep),3));t=q[faces];fn=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);np.add.at(n,physical[faces.ravel()],np.repeat(fn,3,axis=0));n=n[physical];norm=np.linalg.norm(n,axis=1)
    if np.any(norm<1e-20):raise ValueError('Undefined emitted normal')
    n=np.rint(n/norm[:,None]*32767).astype('<i2');emitted={};offset=0
    for row,m in zip(selected,meshes):
        end=offset+len(m.vertices);attrs=dict(m.attributes,normals=n[offset:end],normalSigns=np.ones(len(m.vertices),'i1'))
        final=Mesh(q[offset:end],m.faces,m.frame,attrs,m.source);file=out/f'part-{len(emitted):02}.npz';final.save(file)
        emitted[row['partId']]=dict(row,path=str(file),sha256=sha(file),geometrySha256=final.geometry_sha,
                                   sourcePath=row['path'],sourceSha256=row['sha256']);offset=end
    if time.monotonic()-started>18:raise TimeoutError('18-second named-surface carriage stage')
    report=copy.deepcopy(data['source']);report.update(rows=[emitted.get(r['partId'],r)for r in report['rows']],
        complete=False,constructionComplete=True,installationAllowed=False,installed=False,skeletonChanged=False,
        surfaceCarriage=proof,carriers=evidence,rigidSourceClosures=len(sections),sectionGaugeError_m=gauge,
        pending=['Whole vascular crossings and actual bone/skin/tissue replay','Packed reciprocal mouths and native Blender sections'],
        inputs={k:dict(path=str(p),sha256=sha(p))for k,p in paths.items()},sourceCodeSha256=sha(__file__),elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json',report);print(json.dumps(dict(carriers=len(carriers),rigidSourceClosures=len(sections),gaugeError_m=gauge,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
