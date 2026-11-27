"""Smooth an assembled exterior without blending it into dissected interfaces.

Coordinates and triangles are unchanged. Every reciprocal interface uses
one normal chart with opposite signs. A sharp/re-entrant corner that cannot
use the averaged outward hemisphere retains its actual facet normal.
"""
from collections import defaultdict
import numpy as np


def shared_tissue_corner_normals(positions,cells,interface_channels=None):
    p=np.asarray(positions,float);cells={name:np.asarray(f,np.int64)for name,f in cells.items()}
    channels={}if interface_channels is None else {tuple(sorted(k)):v for k,v in interface_channels.items()}
    owners=defaultdict(list)
    for name,faces in cells.items():
        for i,t in enumerate(faces):owners[tuple(sorted(map(int,t)))].append((name,i,t))
    sums={};assignments=[];reciprocal_faces=0
    for values in owners.values():
        if len(values)==1:channel='outer';canonical=values[0];signs=[1]
        elif len(values)==2:
            values=sorted(values,key=lambda r:r[0]);a,b=values
            if a[0]==b[0]or not all((int(y),int(x))in list(zip(b[2],np.roll(b[2],-1)))for x,y in zip(a[2],np.roll(a[2],-1))):
                raise ValueError('A shared triangle must have exactly two opposed tissue owners')
            channel=channels.get((a[0],b[0]),'interface:'+a[0]+':'+b[0]);canonical=a;signs=[1,-1];reciprocal_faces+=1
        else:raise ValueError('Three or more tissues own one triangle')
        tri=p[canonical[2]];cross=np.cross(tri[1]-tri[0],tri[2]-tri[0])
        if not np.linalg.norm(cross)>0:raise ValueError('Zero-area source triangle')
        if channel not in sums:sums[channel]=np.zeros_like(p)
        for vertex in canonical[2]:sums[channel][vertex]+=cross
        assignments.append((channel,canonical,values,signs,cross))
    for channel,n in sums.items():
        length=np.linalg.norm(n,axis=1);sums[channel]=np.divide(n,length[:,None],out=np.zeros_like(n),where=length[:,None]>0)
    result={name:np.empty((len(f),3,3))for name,f in cells.items()};fallbacks=0
    for channel,canonical,values,signs,cross in assignments:
        ids=canonical[2];n=sums[channel][ids].copy();bad=n@cross<=0
        if bad.any():n[bad]=cross/np.linalg.norm(cross);fallbacks+=int(bad.sum())
        by_vertex={int(i):normal for i,normal in zip(ids,n)}
        for (name,fi,tri),sign in zip(values,signs):result[name][fi]=np.array([by_vertex[int(i)]*sign for i in tri])
    for name,f in cells.items():
        tri=p[f];cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        if np.any(np.einsum('fki,fi->fk',result[name],cross)<=0):raise ValueError('A generated corner normal is not outward')
    return result,dict(shared_interface_triangles=reciprocal_faces,normal_charts=len(sums),sharp_corner_facet_fallbacks=fallbacks,
        all_emitted_corner_normals_outward=True,reciprocal_normals_literal_sign_inversion=True,geometry_changed=False,anatomical_acceptance=False)
