"""Literal half-graphs for controlled apposition of a tube's own lining.

The proposed displacement is a shear within each occupied half-space. Its
continuous Jacobian determinant is one; emitted triangles still need their own
inversion, crossing, normal and occupied-volume checks. Full closure additionally
requires separate apposing tissue cells and a state-specific fluid-boundary proof.
"""
import numpy as np
from geometry_tools.model import Mesh


def clip_indexed(vertices, faces, axis, limit, keep_above, parents=None):
    """Clip triangles with shared edge controls; introduce no lumen cap."""
    vertices=np.asarray(vertices,float);faces=np.asarray(faces,int)
    if axis not in (0,1,2)or not np.isfinite(limit):raise ValueError('Explicit Cartesian clipping plane required')
    points=list(vertices);distance=vertices[:,axis]-limit
    if not keep_above:distance=-distance
    cache={};out=[];source=[]
    parents=np.arange(len(faces))if parents is None else np.asarray(parents,int)
    def crossing(a,b):
        if distance[a]==0:return int(a)
        if distance[b]==0:return int(b)
        key=tuple(sorted((int(a),int(b))))
        if key not in cache:
            lo,hi=key;fraction=distance[lo]/(distance[lo]-distance[hi]);q=vertices[lo]+fraction*(vertices[hi]-vertices[lo]);q[axis]=limit
            cache[key]=len(points);points.append(q)
        return cache[key]
    for face,parent in zip(faces,parents):
        polygon=[]
        for a,b in zip(face,np.roll(face,-1)):
            inside_a=distance[a]>=0;inside_b=distance[b]>=0
            if inside_a:polygon.append(int(a))
            if inside_a!=inside_b:polygon.append(crossing(a,b))
        clean=[]
        for q in polygon:
            if not clean or clean[-1]!=q:clean.append(q)
        if len(clean)>1 and clean[0]==clean[-1]:clean.pop()
        for k in range(1,len(clean)-1):
            tri=[clean[0],clean[k],clean[k+1]];v=np.asarray([points[i]for i in tri])
            if np.linalg.norm(np.cross(v[1]-v[0],v[2]-v[0]))>0:out.append(tri);source.append(int(parent))
    if not out:return np.empty((0,3)),np.empty((0,3),int),np.empty(0,int)
    out=np.asarray(out,int);used,inverse=np.unique(out,return_inverse=True)
    return np.asarray(points)[used],inverse.reshape(-1,3),np.asarray(source,int)


def lining_half_graphs(vertices,faces,basis,origin,axial_limits,frame,*,partition="median"):
    """Retain and certify each physical half-lining as a projected graph."""
    from geometry_tools.contact import pair_contacts
    from geometry_tools.audit import topology
    vertices=np.asarray(vertices,float);faces=np.asarray(faces,int);basis=np.asarray(basis,float);origin=np.asarray(origin,float)
    if basis.shape!=(3,3)or origin.shape!=(3,)or not np.allclose(basis@basis.T,np.eye(3),atol=1e-14,rtol=0)or np.linalg.det(basis)<=0:
        raise ValueError('A proper explicit orthonormal chart is required')
    low,high=map(float,axial_limits)
    if not np.isfinite([low,high]).all()or low>=high:raise ValueError('Finite ordered axial bounds required')
    local=(vertices-origin)@basis.T;tri=local[faces];keep=(tri[:,:,1].max(1)>low)&(tri[:,:,1].min(1)<high)
    used,index=np.unique(faces[keep],return_inverse=True);v=local[used];f=index.reshape(-1,3);source_indices=np.flatnonzero(keep)
    if partition not in ('median','normal-facing'):raise ValueError('Unknown lining partition')
    orientation=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])[:,0]
    if partition=='normal-facing' and np.any(orientation==0):
        raise ValueError('A source facet is vertical in the proposed graph chart')
    result={}
    for name,positive in [('positive',True),('negative',False)]:
        if partition=='median':q,g,parent=clip_indexed(v,f,0,0.,positive,source_indices)
        else:
            select=orientation<0 if positive else orientation>0
            q=v;g=f[select];parent=source_indices[select]
        q,g,parent=clip_indexed(q,g,1,low,True,parent)
        q,g,parent=clip_indexed(q,g,1,high,False,parent)
        if not len(g):raise ValueError('A half-lining is empty')
        projected=q.copy();projected[:,0]=0.;cross=np.cross(projected[g[:,1]]-projected[g[:,0]],projected[g[:,2]]-projected[g[:,0]])[:,0]
        expected=-1 if positive else 1
        if np.any(expected*cross<=0):raise ValueError(name+' lining has a vertical or folded projection')
        mesh=Mesh(projected,g,frame);top=topology(mesh)
        contacts,summary=pair_contacts(mesh,mesh,0.,selfcheck=True,max_pairs=20000,max_tests=1000000)
        forbidden=np.isin(contacts['relation'],['crossing','coplanar-overlap'])
        if summary['truncated']or forbidden.any():
            error=ValueError(name+' projected lining overlaps or its check is incomplete')
            error.witness=dict(side=name,summary=summary,positions=q,faces=g,source_faces=parent,
                               contacts={k:np.asarray(v)[forbidden]for k,v in contacts.items()})
            raise error
        result[name]=dict(mesh=mesh,positions=q,faces=g,source_faces=parent,
                          report=dict(vertices=len(q),triangles=len(g),topology=top,projected_intersections=summary,
                                      projected_crossings=int(forbidden.sum()),passed=True))
    return result


def quintic_window(values,center,half_contact,half_support):
    values=np.asarray(values,float)
    if not 0<=half_contact<half_support:raise ValueError('A contact core inside a wider support is required')
    t=np.clip((abs(values-center)-half_contact)/(half_support-half_contact),0,1)
    return (1-t)**3*(1+3*t+6*t*t)


def quintic_window_derivative(values,center,half_contact,half_support):
    values=np.asarray(values,float)
    if not 0<=half_contact<half_support:raise ValueError('A contact core inside a wider support is required')
    t=np.clip((abs(values-center)-half_contact)/(half_support-half_contact),0,1)
    return -30*t*t*(1-t)**2*np.sign(values-center)/(half_support-half_contact)
