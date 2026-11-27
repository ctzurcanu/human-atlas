"""Unfold an exactly planar internal interface after coordinate storage.

Only explicitly mutable interfaces are eligible. Their complete directed
boundary, holes and coordinates stay literal; no anatomical surface is filled
or projected. A connected embedded planar patch has one normal orientation.
Mixed orientations indicate a folded triangulation of the same boundary.
"""
from collections import Counter
import numpy as np
import manifold3d as md
from .annular_pipe_envelopes import boundary_loops,face_components


def edge_chain(faces):
    count=Counter((int(a),int(b))for t in faces for a,b in zip(t,np.roll(t,-1)))
    return {e:n-count[e[::-1]]for e,n in count.items()if n!=count[e[::-1]]}


def retriangulate_axis_interfaces(positions,surfaces,mutable_surfaces):
    from geometry_tools.model import Mesh
    from geometry_tools.contact import pair_contacts
    v=np.asarray(positions,float);out={sid:np.array(f,copy=True)for sid,f in surfaces.items()};changes=[]
    for sid in sorted(set(mutable_surfaces)):
        if sid not in out:raise ValueError('Unknown explicitly mutable interface: '+sid)
        faces=out[sid];tri=v[faces];remove=set();replacement=[]
        for axis in range(3):
            planar=np.all(tri[:,:,axis]==tri[:,:1,axis],axis=1)
            for level in np.unique(tri[planar,0,axis]):
                ids=np.flatnonzero(planar&np.all(tri[:,:,axis]==level,axis=1));part=faces[ids]
                count,labels=face_components(part)
                for component in range(count):
                    selected=ids[labels==component];old=faces[selected];p=v[old]
                    normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])[:,axis]
                    if not(np.any(normal>0)and np.any(normal<0)):continue
                    loops=boundary_loops(old)
                    if not loops:raise ValueError('A folded plane has no declared boundary')
                    axes=[(axis+1)%3,(axis+2)%3];origin=v[np.unique(old)].mean(0)
                    polygons=[(v[loop]-origin)[:,axes]for loop in loops]
                    signed=sum(float((q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)).sum()/2)for q in polygons)
                    if signed==0:raise ValueError('A whole planar interface collapsed to zero signed area')
                    sign=1 if signed>0 else -1
                    ordered=loops if sign==1 else[loop[::-1]for loop in loops]
                    polygons=[(v[loop]-origin)[:,axes]for loop in ordered]
                    indices=np.concatenate(ordered);new=indices[np.asarray(md.triangulate(polygons,epsilon=0.,allow_convex=False),int)]
                    if sign<0:new=new[:,[0,2,1]]
                    if edge_chain(new)!=edge_chain(old):raise ValueError('Planar retriangulation changed a literal directed boundary edge')
                    q=v[new];area=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])[:,axis]
                    if np.any(sign*area<=0):raise ValueError('The declared boundary has no consistently oriented planar triangulation')
                    _,check=pair_contacts(Mesh(v,new,'stored:metres'),Mesh(v,new,'stored:metres'),0.,selfcheck=True,max_pairs=1000000)
                    if check['truncated']or check['crossing_face_pairs']or check['coplanar_overlap_face_pairs']:
                        raise ValueError('The replacement planar interface is not embedded')
                    remove.update(map(int,selected));replacement.append(new)
                    changes.append(dict(surface=sid,axis=axis,coordinate_m=float(level),old_faces=len(old),new_faces=len(new),
                        boundary_vertices=sum(map(len,loops)),boundary_loops=len(loops),signed_area_m2=signed,
                        original_coordinates_literal=True,directed_boundary_literal=True,self_intersections=check))
        if remove:out[sid]=np.concatenate([faces[[i for i in range(len(faces))if i not in remove]],*replacement])
    return out,dict(changes=changes,coordinates_changed=False,only_declared_exact_axis_interfaces=True)
