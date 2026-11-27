"""Source-face side portal and a rounded, section-preserving arterial daughter.

The parent keeps every original control and all faces except one explicitly
replaced portal face. No disk spans the daughter mouth. This is an outer
vascular envelope, not a claim of finite histological vessel walls.
"""
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from shapely import Polygon,constrained_delaunay_triangles
from geometry_tools.model import Mesh
from geometry_tools.commands.filter_hand_components import encoded_attribute
from .tube_sections import rotation_minimizing_frames,unit

def make_side_branch(parent,face_id,end_point,end_tangent,*,radius_m=.00015,spacing_m=.0007):
    if not 0<radius_m<=.0005 or not 0<spacing_m<=.002:raise ValueError('Explicit bounded branch dimensions required')
    old_v=parent.vertices;old_f=parent.faces;tri=old_v[old_f[face_id]];winding_n=unit(np.cross(tri[1]-tri[0],tri[2]-tri[0]))
    lengths=np.linalg.norm(tri[[1,2,0]]-tri[[2,0,1]],axis=1);centre=np.einsum('i,ij->j',lengths,tri)/lengths.sum();inradius=np.linalg.norm(np.cross(tri[1]-tri[0],tri[2]-tri[0]))/lengths.sum()
    if radius_m>inradius*.8:raise ValueError('Requested branch calibre does not fit this source portal face')
    words=encoded_attribute('normals',parent.attributes['normals'],np.dtype('<i2'));signed=words.astype(float)*parent.attributes['normalSigns'][:,None];outward=unit(signed[old_f[face_id]].mean(0));orientation=1 if outward@winding_n>=0 else -1;n=winding_n*orientation
    end=np.asarray(end_point,float);terminal=unit(end_tangent);distance=np.linalg.norm(end-centre)
    if distance<radius_m*5 or distance>.08 or np.dot(end-centre,n)<.25*distance:raise ValueError('A finite outward daughter course is required')
    u=unit(tri[1]-tri[0]);v=np.cross(n,u);basis=np.column_stack((u,v));theta=np.arange(6)*np.pi/3
    profile=radius_m*np.c_[np.cos(theta),np.sin(theta)];inner=centre+profile@basis.T;outer=(tri-centre)@basis
    region=Polygon(outer,[profile]);assert region.is_valid
    coords=np.r_[outer,profile];global_ids=np.r_[old_f[face_id],np.arange(len(old_v),len(old_v)+6)];annulus=[]
    for piece in constrained_delaunay_triangles(region).geoms:
        xy=np.array(piece.exterior.coords)[:3];idx=[]
        for p in xy:
            d=np.linalg.norm(coords-p,axis=1);j=int(d.argmin());assert d[j]<1e-14;idx.append(int(global_ids[j]))
        points=np.r_[old_v,inner][idx]
        if np.cross(points[1]-points[0],points[2]-points[0])@winding_n<0:idx[1],idx[2]=idx[2],idx[1]
        annulus.append(idx)
    annulus=np.array(annulus,int)
    assert abs(sum(Polygon(np.array([coords[np.flatnonzero(global_ids==j)[0]]for j in f])).area for f in annulus)-region.area)<1e-15
    spline=CubicHermiteSpline([0.,1.],np.array([centre,end]),np.array([n,terminal])*distance*.65)
    count=min(160,max(12,int(np.ceil(distance/spacing_m))+1));ts=np.linspace(0,1,count);path=spline(ts);tangents=unit(spline(ts,1));frames=rotation_minimizing_frames(path,tangents,initial_reference=u)
    turns=np.arccos(np.clip(np.einsum('ij,ij->i',tangents[:-1],tangents[1:]),-1,1))
    if turns.max()>np.deg2rad(15):raise ValueError('Branch station turn exceeds 15 degrees')
    rings=path[:,None,:]+np.einsum('nij,kj->nki',frames[:,:,:2],profile);rings[0]=inner
    all_v=np.r_[old_v,rings.reshape(-1,3)];branch_faces=[];offset=len(old_v)
    for k in range(count-1):
        for j in range(6):
            a=offset+k*6+j;b=offset+k*6+(j+1)%6;c=b+6;d=a+6;branch_faces.extend([[a,b,c],[a,c,d]])
    branch_faces=np.array(branch_faces,int)
    if orientation<0:branch_faces=branch_faces[:,[0,2,1]]
    parent_faces=np.r_[old_f[np.arange(len(old_f))!=face_id],annulus];all_f=np.r_[parent_faces,branch_faces]
    # Angle weighted physical normals on the actual joined envelope. Original
    # source words stay literal; new portal words are shared by both pieces.
    physical,inv=np.unique(all_v,axis=0,return_inverse=True);f=inv[all_f];t=physical[f];fn=unit(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]))*orientation;normal=np.zeros_like(physical)
    for k in range(3):
        a=t[:,(k+1)%3]-t[:,k];b=t[:,(k+2)%3]-t[:,k];angle=np.arctan2(np.linalg.norm(np.cross(a,b),axis=1),np.einsum('ij,ij->i',a,b));np.add.at(normal,f[:,k],fn*angle[:,None])
    new_normals=unit(normal[inv[len(old_v):]]);attrs={}
    for key,value in parent.attributes.items():
        if key=='normals':attrs[key]=np.r_[words,np.rint(new_normals*32767).astype('<i2')]
        elif key=='normalSigns':attrs[key]=np.r_[value,np.ones(len(all_v)-len(old_v),np.int8)]
        else:attrs[key]=np.concatenate([value,np.repeat(value[old_f[face_id][:1]],len(all_v)-len(old_v),axis=0)])
    joined=Mesh(all_v,all_f,parent.frame,attrs,parent.source);patched=Mesh(all_v[:len(old_v)+6],parent_faces,parent.frame,{k:a[:len(old_v)+6]for k,a in attrs.items()},parent.source)
    used,local=np.unique(branch_faces,return_inverse=True);daughter=Mesh(all_v[used],local.reshape(-1,3),parent.frame,{k:a[used]for k,a in attrs.items()},parent.source)
    areas=np.linalg.norm(np.cross(rings-rings.mean(1)[:,None],np.roll(rings-rings.mean(1)[:,None],-1,axis=1)).sum(1),axis=1)/2
    assert np.max(abs(areas/areas[0]-1))<1e-10
    proof=dict(sourceFaceId=int(face_id),sourceFacesRetained=len(old_f)-1,sourceControlsLiteral=True,portalCentre_m=centre.tolist(),portalRadius_m=radius_m,portalInradius_m=float(inradius),parentPortalRawIds=list(range(len(old_v),len(old_v)+6)),daughterPortalRawIds=list(range(6)),stations=count,maximumStationTurn_deg=float(np.rad2deg(turns.max())),maximumRelativeSectionAreaError=float(np.max(abs(areas/areas[0]-1))),endPoint_m=end.tolist(),endTangent=terminal.tolist(),lumenDisksCreated=0,finiteHistologicalWallEstablished=False)
    return patched,daughter,joined,proof
