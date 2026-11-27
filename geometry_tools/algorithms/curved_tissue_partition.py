"""Conforming cell partitions about an explicitly supplied smooth graph.

Anatomical ownership and the graph come from the caller. This module supplies
neither fissure locations nor tissue thicknesses. Original carrier triangles
are partitioned literally; new graph controls refine only the shared interior
interface. Its two tissue owners receive the exact same faces in reverse.
"""
from dataclasses import dataclass
from collections import Counter
import numpy as np
from .shared_closed_cell_clip import SharedCutBank, split_closed_cell, _volume


@dataclass(frozen=True)
class GraphPartition:
    positions: np.ndarray
    cells: dict
    interfaces: dict
    evidence: dict


def _coordinates(origin,basis):
    origin=np.asarray(origin,float);basis=np.asarray(basis,float)
    if origin.shape!=(3,) or basis.shape!=(3,3) or not np.isfinite(np.r_[origin,basis.ravel()]).all():
        raise ValueError('Finite origin and a proper orthonormal graph basis required')
    if not np.allclose(basis@basis.T,np.eye(3),atol=1e-12,rtol=0) or np.linalg.det(basis)<.999999999999:
        raise ValueError('Graph coordinates require a proper rotation, without scale or reflection')
    return origin,basis


class _GraphCutBank(SharedCutBank):
    def __init__(self,positions,origin,basis,height):
        world=np.asarray(positions,float)
        local=(world-origin)@basis.T
        h=np.asarray(height(local[:,:2]),float)
        if h.shape!=(len(world),) or not np.isfinite(h).all():
            raise ValueError('Graph must return one finite height per control')
        local[:,2]-=h
        super().__init__(local)
        self.world=world.tolist()

    def intersection(self,a,b,normal,offset,plane_id):
        before=len(self.positions)
        index=super().intersection(a,b,normal,offset,plane_id)
        if len(self.positions)>before:
            row=self.intersections[-1];x,y=row['edge'];t=row['fraction']
            p=np.asarray(self.world[x])*(1-t)+np.asarray(self.world[y])*t
            assert index==len(self.world)
            self.world.append(p.tolist())
        return index


def _curved_cap(positions,faces,origin,basis,height,maximum_edge_m,max_faces):
    """One bounded CDT of explicit interior samples; no boundary Steiner points.

    Resolution controls interior sampling, not the unaltered carrier rim. Long
    carrier edges stay long and their graph approximation remains in evidence.
    """
    import meshpy.triangle as triangle
    f=np.asarray(faces,int).copy();initial=Counter(tuple(sorted((int(a),int(b))))for t in f for a,b in zip(t,np.roll(t,-1)))
    if any(n>2 for n in initial.values()):raise ValueError('Nonmanifold graph interface')
    boundary={e for e,n in initial.items() if n==1};rim=np.array(sorted({i for edge in boundary for i in edge}),int)
    v=np.asarray(positions,float);xy=((v-origin)@basis.T)[:,:2];patch=xy[f]
    spacing=maximum_edge_m/np.sqrt(2.)
    if spacing<=0:raise ValueError('Curved-interface sampling step underflows')
    lo=xy[rim].min(0);hi=xy[rim].max(0)
    requested=(hi-lo)/spacing
    if not np.isfinite(requested).all()or np.any(requested>max_faces*2):
        raise ValueError('Curved-interface sample budget exceeded before allocation')
    counts=np.maximum(1,np.ceil(requested).astype(int))
    if int(np.prod(counts))>max_faces*2:raise ValueError('Curved-interface sample budget exceeded')
    grid=np.array(np.meshgrid(*[np.linspace(lo[i],hi[i],counts[i]+1)[1:-1]for i in [0,1]],indexing='ij')).reshape(2,-1).T
    def inside(points):
        keep=np.zeros(len(points),bool)
        for a,b,c in patch:
            ab,ac=b-a,c-a;den=ab[0]*ac[1]-ab[1]*ac[0]
            if den<=0:raise ValueError('Graph projection reverses an interface triangle')
            q=points-a;beta=(q[:,0]*(c-a)[1]-q[:,1]*(c-a)[0])/den
            gamma=((b-a)[0]*q[:,1]-(b-a)[1]*q[:,0])/den
            keep|=(beta>=-1e-12)&(gamma>=-1e-12)&(beta+gamma<=1+1e-12)
        return keep
    grid=grid[inside(grid)]
    # Keep samples away from unrefined rim segments and coincident controls.
    keep=np.ones(len(grid),bool)
    for a,b in boundary:
        edge=xy[b]-xy[a];t=np.clip((grid-xy[a])@edge/(edge@edge),0,1)
        keep&=np.linalg.norm(grid-xy[a]-t[:,None]*edge,axis=1)>spacing*.2
    grid=grid[keep];points=np.concatenate([xy[rim],grid])
    if len(points)*2>max_faces+len(rim)+2:raise ValueError('Curved-interface face budget exceeded before triangulation')
    lookup={int(i):j for j,i in enumerate(rim)}
    info=triangle.MeshInfo();info.set_points(points.tolist())
    info.set_facets([[lookup[a],lookup[b]]for a,b in sorted(boundary)])
    result=triangle.build(info,quality_meshing=False,allow_boundary_steiner=False,allow_volume_steiner=False)
    out_points=np.asarray(result.points,float);out_faces=np.asarray(result.elements,int)
    if len(out_faces)>max_faces:raise ValueError('Curved-interface face budget exceeded')
    original_points={tuple(q):j for j,q in enumerate(points)}
    if any(tuple(q)not in original_points for q in out_points):raise ValueError('Unexpected Steiner control')
    if not inside(out_points[out_faces].mean(1)).all():raise ValueError('Triangulation leaves the original interface domain')
    inserted=[];h=np.asarray(height(grid),float)
    if h.shape!=(len(grid),)or not np.isfinite(h).all():raise ValueError('Nonfinite graph refinement')
    for q,w in zip(grid,h):
        point=np.r_[q,w]@basis+origin;inserted.append(len(positions));positions.append(point.tolist())
    world_ids=np.r_[rim,np.array(inserted,int)]
    remap=np.array([world_ids[original_points[tuple(q)]]for q in out_points],int)
    f=remap[out_faces]
    final=Counter(tuple(sorted((int(a),int(b))))for t in f for a,b in zip(t,np.roll(t,-1)))
    if {e for e,n in final.items() if n==1}!=boundary:raise ValueError('Graph refinement changed a carrier rim')
    return f,dict(inserted_controls=inserted,triangles=len(f),original_rim_edges_preserved=True,
        method='bounded interior grid and constrained Delaunay; zero boundary/automatic Steiner points',
        interior_sample_spacing_m=spacing,full_surface_spacing_not_guaranteed_at_fixed_rim=True)


def split_graph_cells(positions,cells,*,origin,basis,height,plane_id,
                      maximum_interface_edge_m=.008,maximum_interface_faces=8000):
    """Split *all* supplied cells through one shared position/edge bank.

    The returned keys are ``(original_name, -1 or +1)``. Nonintersected cells
    retain their original faces. Passing all cells sharing a cut edge is
    mandatory; proximity does not infer that dependency. A graph is sampled
    linearly on each existing carrier edge, so its rim approximation error is
    reported separately. Full triangle, containment and anatomical review are
    required before using the candidate as tissue.
    """
    origin,basis=_coordinates(origin,basis)
    p=np.asarray(positions,float)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite shared world positions required')
    if not np.isfinite(maximum_interface_edge_m)or maximum_interface_edge_m<=0 or maximum_interface_faces<4:
        raise ValueError('Finite positive interface resolution and face budget required')
    bank=_GraphCutBank(p,origin,basis,height);normal=np.array([0.,0.,1.]);outputs={};splits={}
    for name,faces in cells.items():
        faces=np.asarray(faces,int);used=np.unique(faces);d=np.asarray(bank.positions)[used,2]
        if d.min()>=-1e-12:outputs[(name,1)]=faces.copy();continue
        if d.max()<=1e-12:outputs[(name,-1)]=faces.copy();continue
        split=split_closed_cell(bank,faces,normal,0.,plane_id)
        splits[name]=split
        outputs[(name,-1)]=split.negative_faces
        outputs[(name,1)]=split.positive_faces
    world=bank.world;interfaces={};proof={}
    for name,split in splits.items():
        cap,refined=_curved_cap(world,split.interface_faces,origin,basis,height,
            maximum_interface_edge_m,maximum_interface_faces)
        count=len(split.interface_faces)
        outputs[(name,-1)]=np.concatenate([split.negative_faces[:-count],cap])
        outputs[(name,1)]=np.concatenate([split.positive_faces[:-count],cap[:,[0,2,1]]])
        interfaces[name]=cap
        v=np.asarray(world,float);original=np.asarray(cells[name],int)
        original_volume=_volume(p,original)
        child_volumes=[_volume(v,outputs[(name,s)]) for s in [-1,1]]
        if min(child_volumes)<=0:raise ValueError('Curved partition has a nonpositive cell volume')
        if abs(sum(child_volumes)-original_volume)>max(1e-16,abs(original_volume)*1e-9):
            raise ValueError('Curved partition lost whole-cell volume')
        fragments=np.concatenate([split.negative_faces[:-count],split.positive_faces[:-count]])
        owners=np.concatenate([split.negative_ancestry[:-count],split.positive_ancestry[:-count]])
        t=v[fragments];new_area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
        t=p[original];old_area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
        error=np.bincount(owners,weights=new_area,minlength=len(original))-old_area
        if max(abs(error))>1e-14:raise ValueError('Curved split changed the literal outer carrier')
        rim=np.unique(np.array(split.boundary_loops,dtype=object).tolist()[0]) if len(split.boundary_loops)==1 else np.unique(np.concatenate(split.boundary_loops))
        local=(v[rim]-origin)@basis.T
        rim_error=abs(local[:,2]-height(local[:,:2]))
        proof[name]=dict(original_controls_exact=True,whole_volume_m3=original_volume,
            child_volumes_m3=child_volumes,volume_difference_m3=sum(child_volumes)-original_volume,
            maximum_source_face_area_error_m2=float(max(abs(error))),
            maximum_sampled_graph_rim_error_m=float(max(rim_error)),
            reciprocal_interface_literal=True,refinement=refined,
            graph_sampled_linearly_on_carrier_edges=True,anatomical_acceptance=False)
    v=np.asarray(world,float)
    assert np.array_equal(v[:len(p)],p)
    return GraphPartition(v,outputs,interfaces,dict(cells=proof,source_controls=len(p),
        intersection_controls=bank.intersections,all_shared_input_cells_required=True,
        requires_independent_mesh_and_anatomical_review=True))
