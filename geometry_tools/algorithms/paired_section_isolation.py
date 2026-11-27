"""Isolate a returning/branched finite wall at its actual paired section.

The remainder of the organ may cross the same infinite plane elsewhere.
Consequently a half-space component is not the requested anatomical branch.
Only the two section cycles enclosing the explicit fluid query are barriers;
all other plane cuts are reconnected through their original shared edges.
This supplies geometry, never an anatomical branch identity.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .scalar_surface_clip import clip_scalar_surface
from .annular_pipe_envelopes import boundary_loops
from .simple_compartment_section import polygon_contains


def isolate_paired_section(positions,faces,point,normal,*,restore_unrelated=False):
    """Return a source-exact subdivision and two face components.

    `point` must be explicitly inside the desired fluid section, between two
    paired wall contours. No nearest component, cap, glue, or tolerance-based
    vertex welding is introduced. The original position bank remains literal.
    With `restore_unrelated`, retain each original input face whose children
    all belong to the same component. Only the selected section remains
    subdivided; unrelated native interfaces retain their literal triangles.
    The caller owns source-mouth membership, closed-volume and anatomy checks.
    """
    v=np.asarray(positions,float);f=np.asarray(faces)
    point=np.asarray(point,float);normal=np.asarray(normal,float)
    if point.shape!=(3,)or normal.shape!=(3,)or not np.isfinite(np.r_[point,normal]).all()or np.linalg.norm(normal)==0:
        raise ValueError('A finite fluid query and nonzero plane normal are required')
    normal=normal/np.linalg.norm(normal)
    axis=np.eye(3)[np.argmin(abs(normal))];u=np.cross(normal,axis);u/=np.linalg.norm(u)
    basis=np.c_[u,np.cross(normal,u)]
    cut=clip_scalar_surface(v,f,(v-point)@normal,0.)
    contours=boundary_loops(cut['positive_faces'])
    owning=[r for r in contours if polygon_contains((cut['positions'][r]-point)@basis,np.zeros(2))]
    if len(owning)!=2:
        raise ValueError(f'Fluid query must own exactly two wall contours; found {len(owning)}')
    def area(r):
        xy=(cut['positions'][r]-point)@basis;next_xy=np.roll(xy,-1,axis=0)
        return abs((xy[:,0]*next_xy[:,1]-xy[:,1]*next_xy[:,0]).sum())
    owning.sort(key=area)
    barriers={tuple(sorted((int(a),int(b))))for r in owning for a,b in zip(r,np.roll(r,-1))}
    all_faces=np.r_[cut['positive_faces'],cut['negative_faces']]
    parents=np.r_[cut['positive_parents'],cut['negative_parents']]
    edges=np.vstack([all_faces[:,[0,1]],all_faces[:,[1,2]],all_faces[:,[2,0]]])
    owners=np.tile(np.arange(len(all_faces)),3)
    unique,inverse,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    if np.any(count!=2):raise ValueError('The supplied finite wall must remain closed after subdivision')
    keep=np.array([tuple(edge)not in barriers for edge in unique])
    if np.count_nonzero(~keep)!=len(barriers):raise ValueError('Paired section barrier edges are incomplete')
    order=np.argsort(inverse,kind='stable');starts=np.r_[0,np.cumsum(count)[:-1]][keep]
    a=owners[order[starts]];b=owners[order[starts+1]]
    graph=coo_matrix((np.ones(2*len(a)),(np.r_[a,b],np.r_[b,a])),shape=(len(all_faces),len(all_faces))).tocsr()
    components,labels=connected_components(graph,directed=False)
    if components!=2:raise ValueError(f'The paired section must split one wall into two components; found {components}')
    restored=0
    if restore_unrelated:
        lower=np.full(len(f),2,int);upper=np.full(len(f),-1,int)
        np.minimum.at(lower,parents,labels);np.maximum.at(upper,parents,labels)
        if np.any(upper<0):raise ValueError('Subdivision lost an input source face')
        homogeneous=lower==upper;original=np.flatnonzero(homogeneous)
        restored=int(np.count_nonzero(homogeneous&(np.bincount(parents,minlength=len(f))>1)))
        retain=~homogeneous[parents]
        all_faces=np.r_[f[original],all_faces[retain]]
        parents=np.r_[original,parents[retain]];labels=np.r_[lower[original],labels[retain]]
        # Every newly cut edge belongs to both incident faces at the selected
        # section. Reject any unintended T-junction after restoring the rest.
        edge=np.vstack([all_faces[:,[0,1]],all_faces[:,[1,2]],all_faces[:,[2,0]]])
        _,count=np.unique(np.sort(edge,axis=1),axis=0,return_counts=True)
        if np.any(count!=2):raise ValueError('Restoring unrelated faces introduced a boundary or nonmanifold edge')
    return dict(positions=cut['positions'],faces=all_faces,source_face_ids=parents,
        component_labels=labels,root_inner=owning[0],root_outer=owning[1],
        root_point=point,root_normal=normal,new_edge_parents=cut['new_edge_parents'],
        new_edge_fractions=cut['fractions'],other_plane_contours_reconnected=len(contours)-2,
        unrelated_subdivided_parents_restored=restored,
        anatomical_caps_added=0,source_controls_literal=True)
