"""Exact directional support planes for literal rigid source-profile groups.

When every group incident to two triangles satisfies its paired support
inequality, the full triangles are separated along that direction. Route
choice and complete anatomical packing still require independent review.
"""
import numpy as np


def minimum_tangent_support_axis(first, second, normal, preferred, *, angular_samples=360):
    """Choose the least separating motion for complete cells in a tangent plane.

    This is a finite, reproducible geometric proposal. It does not infer the
    anatomical side of either vessel. All points of both literal cells enter
    the interval calculation; no face, branch or profile is discarded.
    """
    a,b,n,u=map(lambda x:np.asarray(x,float),(first,second,normal,preferred))
    if a.ndim!=2 or b.ndim!=2 or a.shape[1]!=3 or b.shape[1]!=3 or not len(a) or not len(b):
        raise ValueError('Two nonempty 3D source cells required')
    if not all(np.isfinite(x).all()for x in [a,b,n,u])or angular_samples<4:
        raise ValueError('Finite cells and a bounded angular inventory required')
    n=n/np.linalg.norm(n);u=u-n*(u@n)
    if np.linalg.norm(u)<1e-10:raise ValueError('Preferred axis cannot be normal to the tangent plane')
    u/=np.linalg.norm(u);v=np.cross(n,u)
    angles=np.arange(angular_samples)*2*np.pi/angular_samples
    axes=np.cos(angles)[:,None]*u+np.sin(angles)[:,None]*v
    required=(b@axes.T).max(0)-(a@axes.T).min(0)
    chosen=int(np.argmin(required))
    return axes[chosen],dict(angularSamples=angular_samples,requiredTranslation_m=float(required[chosen]),
                             chosenAngleDegrees=float(angles[chosen]*180/np.pi),anatomicalAcceptance=False)


def split_contact_roles(parts,contacts):
    preserve=[];separate=[]
    for contact in contacts:
        a,b=[parts[contact[k]]['system']for k in ['first','second']]
        if {a,b}=={'arterial','venous'}:separate.append(contact)
        elif a==b and a in ['arterial','venous']:preserve.append(contact)
        else:raise ValueError('Blood contact roles require explicit arterial/venous source systems')
    return preserve,separate


def support_rows(meshes,groups,mappings,routes,*,clearance_m,columns):
    weights=[];normals=[];offsets=[];records=[]
    for route in routes:
        a,b=route['artery'],route['vein'];u=np.array(route['direction']);assert abs(np.linalg.norm(u)-1)<1e-12
        for i,j in route['groupPairs']:
            first=meshes[a].vertices[groups[a][i]]@u;second=meshes[b].vertices[groups[b][j]]@u;w=np.zeros(columns)
            if mappings[a][i]>=0:w[mappings[a][i]]+=1
            if mappings[b][j]>=0:w[mappings[b][j]]-=1
            bound=float(second.max()-first.min()+clearance_m);weights.append(w);normals.append(u);offsets.append(bound);records.append(dict(artery=a,vein=b,arterialGroup=i,venousGroup=j))
    return np.asarray(weights).reshape(-1,columns),np.asarray(normals).reshape(-1,3),np.array(offsets),records


def minimum_support_gap(meshes,groups,routes):
    minimum=np.inf
    for route in routes:
        a,b=route['artery'],route['vein'];u=np.array(route['direction'])
        for i,j in route['groupPairs']:
            minimum=min(minimum,float((meshes[a].vertices[groups[a][i]]@u).min()-(meshes[b].vertices[groups[b][j]]@u).max()))
    return float(minimum)
