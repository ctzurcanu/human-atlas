"""Occupied union for a declared source terminal-to-side interface.

Exactly one prealigned internal terminal ring is explicit; other ports remain literal.
No interface is inferred from proximity. Temporary closures cannot survive.
"""
import numpy as np
import manifold3d as mf
from scipy.spatial import cKDTree
from geometry_tools.model import Mesh
from .tube_continuity import boundary_loops,continuity_audit

def union_terminal_side(mesh,cell_face_ids,root_ids):
    assert len(cell_face_ids)==2 and len(root_ids)==1
    _,rep,physical=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True)
    roots=[np.asarray(r,int)for r in root_ids];rootsets=[{tuple(v)for v in mesh.vertices[r]}for r in roots]
    centers=[mesh.vertices[r].mean(0)for r in roots];origin=np.mean(centers,axis=0);scale=.001
    values=[mesh.vertices];channels={};at=3
    if 'normals'in mesh.attributes:
        n=mesh.attributes['normals'].astype(float)*mesh.attributes['normalSigns'][:,None];n/=np.linalg.norm(n,axis=1)[:,None];values.append(n);channels['normals']=(at,at+3);at+=3
    for key in ['uvs','colors']:
        if key in mesh.attributes:
            v=mesh.attributes[key].astype(float);values.append(v);channels[key]=(at,at+v.shape[1]);at+=v.shape[1]
    props0=np.column_stack(values);cells=[];rootcaps=set();portcaps=set();ports=[];removed=[];internal=0
    for ci,ids in enumerate(cell_face_ids):
        ids=np.asarray(ids,int);removed.extend(ids.tolist());faces=mesh.faces[ids];loops=[rep[loop]for loop in boundary_loops(physical[faces])]
        if len(loops)!=2:raise ValueError('Each declared collar cell must have two source boundaries')
        props=props0.copy();ff=faces.tolist();lineage=ids.tolist();rootends=0
        for li,loop in enumerate(loops):
            coords={tuple(v)for v in mesh.vertices[loop]};match=[i for i,root in enumerate(rootsets)if root==coords];assert len(match)<=1
            capid=1000000+ci*2+li
            if match:center=centers[match[0]];rootcaps.add(capid);rootends+=1
            else:center=mesh.vertices[loop].mean(0);portcaps.add(capid);ports.append(loop)
            index=len(props);props=np.r_[props,np.r_[center,props0[loop,3:].mean(0)][None]]
            for a,b in zip(loop,np.roll(loop,-1)):ff.append([index,int(b),int(a)]);lineage.append(capid)
        assert rootends<=1
        internal+=rootends==2
        used,inv=np.unique(np.array(ff,int),return_inverse=True);p=props[used].copy();p[:,:3]=(p[:,:3]-origin)/scale;f=inv.reshape(-1,3)
        tri=p[f,:3];local=tri-tri.mean((0,1));volume=np.einsum('ij,ij->i',local[:,0],np.cross(local[:,1],local[:,2])).sum()/6
        if volume<0:f=f[:,[0,2,1]]
        _,reps,inverse=np.unique(p[:,:3],axis=0,return_index=True,return_inverse=True);canonical=reps[inverse];seams=np.flatnonzero(canonical!=np.arange(len(p)))
        raw=mf.Mesh64(np.ascontiguousarray(p),np.ascontiguousarray(f,dtype=np.uint64),merge_from_vert=seams.astype(np.uint64),merge_to_vert=canonical[seams].astype(np.uint64),face_id=np.asarray(lineage,dtype=np.uint64),tolerance=0)
        solid=mf.Manifold(raw)
        if solid.status()!=mf.Error.NoError or solid.is_empty():raise ValueError('Source annulus computation solid failed: '+str(solid.status()))
        cells.append(solid)
    union=mf.Manifold.batch_boolean(cells,mf.OpType.Add)
    if union.status()!=mf.Error.NoError:raise ValueError('Coupled source-cell union failed')
    raw=union.to_mesh64();props=np.asarray(raw.vert_properties).copy();props[:,:3]=props[:,:3]*scale+origin;faces=np.asarray(raw.tri_verts,int);lineage=np.asarray(raw.face_id,int)
    parents=np.arange(len(props))
    def find(i):
        while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
        return i
    for a,b in zip(raw.merge_from_vert,raw.merge_to_vert):parents[find(int(a))]=find(int(b))
    canonical=np.array([find(i)for i in range(len(props))]);error=float(np.linalg.norm(props[:,:3]-props[canonical,:3],axis=1).max(initial=0));assert error<1e-11
    props[:,:3]=props[canonical,:3]
    if set(lineage)&rootcaps:raise ValueError('A shared computational root closure survived')
    take=~np.isin(lineage,list(portcaps));faces=faces[take];lineage=lineage[take];assert len(faces)and np.all(lineage<1000000)
    removed=np.unique(removed);controls=np.unique(mesh.faces[removed]);tree=cKDTree(mesh.vertices[controls]);distance,near=tree.query(props[:,:3]);mapping=np.full(len(props),-1,int)
    for i in np.flatnonzero(distance<=1e-11):
        choices=controls[tree.query_ball_point(props[i,:3],1e-11)];compatible=[j for j in choices if np.allclose(props0[j,3:],props[i,3:],atol=1e-7,rtol=1e-7)]
        if compatible:mapping[i]=compatible[0]
    vertices=mesh.vertices.tolist();attributes={k:v.tolist()for k,v in mesh.attributes.items()}
    for i in np.unique(faces):
        if mapping[i]>=0:continue
        mapping[i]=len(vertices);vertices.append(props[i,:3].tolist())
        for key in attributes:
            if key=='normalSigns':value=1
            elif key=='normals':
                lo,hi=channels[key];n=props[i,lo:hi];value=np.rint(n/np.linalg.norm(n)*32767).astype('<i2').tolist()
            elif key in channels:
                lo,hi=channels[key];value=props[i,lo:hi].astype(mesh.attributes[key].dtype).tolist()
            else:raise ValueError('Unsupported source attribute: '+key)
            attributes[key].append(value)
    faces=mapping[faces];_,physical=np.unique(np.asarray(vertices),axis=0,return_inverse=True)
    expected={(int(physical[a]),int(physical[b]))for loop in ports for a,b in zip(loop,np.roll(loop,-1))}
    for sign in [1,-1]:
        f=faces if sign==1 else faces[:,[0,2,1]]
        try:loops=boundary_loops(physical[f])
        except ValueError:continue
        actual={(int(a),int(b))for loop in loops for a,b in zip(loop,np.roll(loop,-1))}
        if actual==expected:faces=f;break
    else:raise ValueError('Coupled carrier lost or altered a literal oriented external port')
    cycles=0;assert internal==0
    euler=2-2*cycles-len(ports);audit=continuity_audit(np.asarray(vertices),faces,expected_mouths=len(ports),expected_components=1,expected_euler=euler)
    if not audit['topological_continuity_pass']:raise ValueError('Coupled carrier changed the declared branch topology')
    candidate=Mesh(np.asarray(vertices),faces,mesh.frame,{k:np.asarray(v,dtype=mesh.attributes[k].dtype)for k,v in attributes.items()},mesh.source)
    return candidate,dict(complete=True,sourceRoots=[r.tolist()for r in roots],sourcePortRawIds=[r.tolist()for r in ports],removedSourceFaceIds=removed.tolist(),newFaceSourceLineage=lineage.tolist(),rootClosuresRetained=0,emittedLumenDisks=0,declaredGraphCycleRank=cycles,carrierAudit=audit,conditioningOrigin_m=origin.tolist(),conditioningScale_m=scale,propertySeamError_m=error,finiteVascularWallEstablished=False,anatomicalAcceptance=False)
