"""Select between explicit overlapping outer-envelope disks using source fields.

This does not identify tissue or infer anatomy. Both caller-declared alternatives
must independently pass the complete redundant-disk proof and share one rim.
Exactly one retained disk must agree with every original outward corner field;
an ambiguous field is rejected instead of selecting a topology-only winner.
"""
import numpy as np
from .redundant_surface_disk import remove_redundant_surface_disk

def source_disk_field(mesh,faces):
    if not {'normals','normalSigns'}<=mesh.attributes.keys():raise ValueError('Original signed normal words required')
    f=mesh.faces[np.asarray(faces,int)];p=mesh.vertices[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);length=np.linalg.norm(n,axis=1)
    w=mesh.attributes['normals'][f].astype(float)*mesh.attributes['normalSigns'][f,None];wl=np.linalg.norm(w,axis=2)
    if np.any(length==0)or np.any(wl==0):raise ValueError('Nondegenerate source faces and fields required')
    cosine=np.einsum('tki,ti->tk',w/wl[:,:,None],n/length[:,None])
    return dict(corners=int(cosine.size),minimumOutwardCosine=float(cosine.min()),nonoutwardCorners=int((cosine<=0).sum()))

def select_outward_source_disk(mesh,choices,*,representation,minimum_outward_cosine):
    if representation!='outer-source-envelope':raise ValueError('Explicit outer source-envelope classification required')
    if len(choices)!=2 or not 0<float(minimum_outward_cosine)<1:raise ValueError('Two explicit alternatives and a positive outward margin required')
    choices=[np.asarray(c,int) for c in choices]
    if set(choices[0])&set(choices[1]):raise ValueError('Alternative source disks must not share faces')
    results=[remove_redundant_surface_disk(mesh,c) for c in choices]
    rings=[{tuple(p) for p in proof['attachmentRingPhysicalPoints_m']} for _,proof in results]
    if rings[0]!=rings[1]:raise ValueError('Alternatives must span the same actual source rim')
    fields=[source_disk_field(mesh,c) for c in choices]
    valid=[i for i in range(2) if fields[1-i]['minimumOutwardCosine']>=minimum_outward_cosine and fields[i]['nonoutwardCorners']>0]
    if len(valid)!=1:raise ValueError('Source outward fields do not establish one unique retained disk')
    i=valid[0];keep,index=results[i]
    return keep,dict(complete=True,removedAlternative=i,retainedAlternative=1-i,sourceAlternatives=[c.tolist() for c in choices],sourceFields=fields,minimumRequiredOutwardCosine=float(minimum_outward_cosine),sourceFieldsUnmodified=True,anatomicalAcceptance=False),index
