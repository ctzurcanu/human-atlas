"""Open an explicitly matched branch root, reconciling one rounded control.

Mother and branch sheets must be complete source charts. Only the declared
six-point mother disk is removed; no branch reconstruction or cap is emitted.
"""
import numpy as np
from geometry_tools.model import Mesh
from .shared_boundary_arc import topology
from .tube_continuity import boundary_loops

def reconcile_source_root_hole(mesh,mother_faces,branch_faces,disk_faces,*,maximum_motion_m=1e-6):
    if not 0<maximum_motion_m<=1e-6:raise ValueError('Source-rounding budget must be at most one micrometre')
    mother,branch,disk=[np.asarray(x,int) for x in [mother_faces,branch_faces,disk_faces]]
    if not np.array_equal(np.sort(np.r_[mother,branch]),np.arange(len(mesh.faces))):raise ValueError('Complete disjoint mother/branch charts required')
    if not set(disk)<=set(mother) or len(np.unique(disk))!=4:raise ValueError('Explicit four-face mother disk required')
    v,alias=np.unique(mesh.vertices,axis=0,return_inverse=True);f=alias[mesh.faces]
    mt,bt,dt=[topology(mesh.vertices,mesh.faces[x]) for x in [mother,branch,disk]]
    if any(t['nonmanifoldEdges'] or t['edgeConnectedComponents']!=1 for t in [mt,bt,dt]):raise ValueError('Input charts must be regular connected surfaces')
    if (mt['eulerCharacteristic'],mt['boundaryLoops'],bt['eulerCharacteristic'],bt['boundaryLoops'])!=(-1,3,0,2):raise ValueError('Mother must have three ends; branch must be one annulus')
    if dt['eulerCharacteristic']!=1 or dt['boundaryLoopLengths']!=[6]:raise ValueError('Mother patch must be one six-point simple disk')
    root=boundary_loops(f[disk])[0]
    roots=[ring for ring in boundary_loops(f[branch]) if len(ring)==6 and len(set(root)&set(ring))==5]
    if len(roots)!=1:raise ValueError('Exactly one branch root with five literal shared controls required')
    br=roots[0];missing=list(set(br)-set(root));target=list(set(root)-set(br))
    if len(missing)!=1 or len(target)!=1:raise ValueError('Exactly one corresponding rounded control required')
    missing,target=missing[0],target[0];distance=float(np.linalg.norm(v[missing]-v[target]))
    if not 0<distance<=maximum_motion_m:raise ValueError('Source-root discrepancy exceeds the fixed rounding budget')
    mapped=np.where(br==missing,target,br)
    if not any(np.array_equal(np.roll(mapped,i),root) for i in range(6)):raise ValueError('Branch and disk contours do not have compatible winding')
    moving=np.flatnonzero(alias==missing)
    if np.isin(mesh.faces[mother],moving).any():raise ValueError('Moving source control belongs to fixed mother')
    # Real mother mouths and the remote branch mouth must remain literal.
    external=boundary_loops(f[mother])+[x for x in boundary_loops(f[branch]) if not np.array_equal(x,br)]
    external_edges={tuple(sorted(e)) for ring in external for e in zip(ring,np.roll(ring,-1))}
    if any(missing in ring for ring in external):raise ValueError('A true external mouth would move')
    q=mesh.vertices.copy();q[moving]=v[target];keep=np.setdiff1d(np.arange(len(f)),disk)
    candidate=Mesh(q,mesh.faces[keep],mesh.frame,mesh.attributes,mesh.source)
    t=topology(q,candidate.faces)
    if t['nonmanifoldEdges'] or t['edgeConnectedComponents']!=1 or t['eulerCharacteristic']!=-2 or t['boundaryLoops']!=4:raise ValueError('Joined source must be one genus-zero four-ended envelope')
    _,after_alias=np.unique(q,axis=0,return_inverse=True)
    # Compare boundary coordinates rather than re-numbered physical vertices.
    edge_points=lambda vv,edges:{tuple(sorted((tuple(vv[a]),tuple(vv[b])))) for a,b in edges}
    after_loops=boundary_loops(after_alias[candidate.faces]);after_v=np.unique(q,axis=0)
    after_edges=[e for ring in after_loops for e in zip(ring,np.roll(ring,-1))]
    if edge_points(v,external_edges)!=edge_points(after_v,after_edges):raise ValueError('True source mouth coordinates changed')
    return candidate,dict(complete=True,retainedSourceFaceIds=keep.tolist(),removedSourceFaceIds=disk.tolist(),
        motherFaces=mother.tolist(),branchFaces=branch.tolist(),motherDiskTopology=dt,
        movingRawIds=moving.tolist(),targetPoint_m=v[target].tolist(),maximumMotion_m=distance,
        sourceBranchContourPreservedWithinRounding=True,originalExternalMouthsLiteral=True,
        newFaces=0,newCaps=0,topology=t,anatomicalAcceptance=False)
