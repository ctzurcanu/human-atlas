"""Transport explicitly selected shared atlas banks without duplicating tissue.

All coordinates and covector normals in a shared position/normal block receive
one mapping. Unselected parts keep their original banks. Topology, ownership,
corner UVs, textures and normal signs are preserved.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from .model import AtlasReader,Mesh
from .common import atomic_bytes,digest_file

def transport(atlas,part_ids,field,output,*,source_scale=1.,source_translation=(0.,0.,0.)):
    reader=AtlasReader();output=Path(output);output.mkdir(parents=True,exist_ok=False)
    selected=[p for p in atlas['parts']if p['id']in set(part_ids)]
    bychunk={};records=[];rewrites={};T=np.array(source_translation,float);S=float(source_scale)
    for p in selected:bychunk.setdefault(p['chunk'],[]).append(p)
    for ci,rows in sorted(bychunk.items()):
        oldurl=atlas['chunks'][ci]['url'];old=reader.asset(oldurl).read_bytes();data=bytearray(old);seen=set();normal_seen=set();local_records=[]
        for p in rows:
            n=p['vertexCount'];o=np.array(p.get('geometryOrigin',[0.,0.,0.]))+np.array(p.get('sourceOffset',[0.,0.,0.]));scale=float(p.get('geometryScale',1.));key=(p['positions'],n,scale,tuple(o))
            original=reader.view(old,p['positions'],n,'<f4',3)
            source=(original.astype(float)*scale+o-T)/S
            mapped,J=field(source);det=np.linalg.det(J);changed=np.any(mapped!=source,axis=1)
            if det.min()<=.05:raise ValueError('Tissue transport folds or nearly collapses')
            if key not in seen:
                local=original.copy();local[changed]=((mapped[changed]*S+T-o)/scale).astype('<f4')
                data[p['positions']:p['positions']+local.nbytes]=local.tobytes();seen.add(key)
                local_records.append(dict(positions=p['positions'],vertices=n,minimum_jacobian=float(det.min()),maximum_displacement_m=float(np.linalg.norm(mapped-source,axis=1).max()*S),fixed_controls=int((~changed).sum())))
            if 'normals'in p:
                nk=(p['normals'],n,key)
                if nk not in normal_seen:
                    words=reader.view(old,p['normals'],n,'<i2',3).copy();normal=np.linalg.solve(J[changed].transpose(0,2,1),words[changed].astype(float)[:,:,None])[:,:,0]
                    normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-30);words[changed]=np.rint(np.clip(normal,-1,1)*32767).astype('<i2')
                    data[p['normals']:p['normals']+words.nbytes]=words.tobytes();normal_seen.add(nk)
        path=output/('bank-'+str(ci)+'.bin');atomic_bytes(path,bytes(data));url='/local-models/'+str(path.relative_to(reader.root/'.local-models'));newci=len(atlas['chunks']);atlas['chunks'].append(dict(url=url,bytes=len(data)))
        for p in rows:
            p['chunk']=newci; m=reader.read_part(atlas,p,frame='transport:assembled-metres');p['bounds']=m.bounds.tolist()
            if 'hepaticAttributeSignature'in p:
                payload=b''.join(bytes(data[p[k]:p[k]+p['vertexCount']*width])for k,width in [('positions',12),('normals',6),('normalSigns',1)])
                p['hepaticAttributeSignature']=hashlib.sha256(payload+json.dumps([p.get('geometryOrigin',[0.,0.,0.]),p.get('geometryScale',1.)]).encode()).hexdigest()
        rewrites[oldurl]=dict(url=url,sha256=digest_file(path),old_chunk=ci,new_chunk=newci,path=str(path));records.append(dict(source=oldurl,source_sha256=hashlib.sha256(old).hexdigest(),destination=url,blocks=local_records,parts=[p['id']for p in rows]))
    return rewrites,records
