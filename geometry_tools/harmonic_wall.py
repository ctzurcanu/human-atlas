"""Peelable coats inside a finite wall, preserving its actual outer/inner skin.

Constrained tetrahedra preserve every source boundary facet. Positive graph
harmonic coordinates interpolate between explicitly named outer/inner carrier
faces; end annuli and tissue attachment rims are free, not false lumen caps.
The tetrahedra are construction-only, never loaded by the atlas renderer.
"""
import itertools
from collections import Counter
import numpy as np
import meshpy.tet as tet
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from .implicit_cells import Arrangement


def mesher_coordinates(source):
    """Use a lossless binary scale and a reversible coarse frame shift.

    Decimal scaling followed by restoring literal boundary controls can turn
    a nearly coplanar construction tetrahedron inside out. Zero origins on
    axes that cannot round-trip avoid sacrificing source coordinates merely
    to centre the numerical input.
    """
    source=np.asarray(source,float)
    origin=np.rint((source.min(0)+source.max(0))/2)
    for axis in range(3):
        if not np.array_equal((source[:,axis]-origin[axis])+origin[axis],source[:,axis]):
            origin[axis]=0.
    local=(source-origin)*1024.
    if not np.array_equal(local/1024.+origin,source):
        raise ValueError('The mesher coordinate transform changed a literal source control')
    return local,origin,1024.


def tetrahedralize(vertices,faces,internal_facets=None):
    faces=np.asarray(faces,int);extra=np.asarray(internal_facets,int).reshape(-1,3)if internal_facets is not None else np.empty((0,3),int)
    all_faces=np.r_[faces,extra];used,inv=np.unique(all_faces,return_inverse=True);source=np.asarray(vertices,float)[used];af=inv.reshape(-1,3);f=af[:len(faces)]
    local,origin,scale=mesher_coordinates(source);info=tet.MeshInfo();info.set_points(local);info.set_facets(af)
    options=tet.Options('pYQ');options.epsilon=1e-15;m=tet.build(info,options=options)
    p=np.asarray(m.points)/scale+origin;t=np.asarray(m.elements,int)
    # TetGen may reorder points. Map literal supplied controls back to their
    # exact original values before checking boundary equality.
    from scipy.spatial import cKDTree
    d,index=cKDTree(source).query(p)
    matched=d<5e-14;p[matched]=source[index[matched]]
    from .algorithms.tetrahedral_orientation import orient_complex
    lookup={tuple(q):i for i,q in enumerate(p)}
    literal_boundary=np.array([[lookup[tuple(source[i])]for i in tri]for tri in f],int)
    t,orientation=orient_complex(p,t,literal_boundary)
    if orientation['flat_interface_glue_tetrahedra']:
        raise ValueError('The mesher produced flat interface transitions; repair its boundary triangulation before constructing tissue')
    q=p[t]-origin;det=np.linalg.det(q[:,1:]-q[:,:1])
    facets=np.concatenate([t[:,r]for r in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)]]);key,count=np.unique(np.sort(facets,axis=1),axis=0,return_counts=True);boundary=key[count==1]
    if np.any(count>2):raise ValueError('Nonmanifold tetrahedral adjacency')
    expected={tuple(sorted(tuple(x)for x in source[tri]))for tri in f};actual={tuple(sorted(tuple(x)for x in p[tri]))for tri in boundary}
    if actual!=expected:raise ValueError('Constrained mesher changed a source boundary facet')
    if len(extra):
        all_keys={tuple(sorted(tuple(x)for x in p[tri]))for tri in key[count==2]}
        if not {tuple(sorted(tuple(x)for x in source[tri]))for tri in af[len(faces):]}<=all_keys:raise ValueError('Constrained mesher changed an internal anatomical interface')
    source_t=source[f]-origin;volume=np.einsum('ij,ij->i',source_t[:,0],np.cross(source_t[:,1],source_t[:,2])).sum()/6
    residual=abs(abs(det).sum()/6-volume)
    if residual>max(1e-17,abs(volume)*1e-9):raise ValueError('Tetrahedra do not replay the actual finite wall volume')
    return p,t,used,index,dict(boundary_faces_literal=True,source_volume_m3=float(volume),tet_volume_residual_m3=float(residual),construction_tetrahedra=len(t),
        coordinate_scale=scale,coordinate_origin=origin.tolist(),source_coordinate_roundtrip_literal=True,mesher_epsilon=options.epsilon,
        orientation_verification=orientation)


def graph_values(points,edges,fixed,fixed_values):
    edges=np.unique(np.sort(edges,axis=1),axis=0)
    lengths=np.linalg.norm(points[edges[:,0]]-points[edges[:,1]],axis=1)
    if np.any(lengths==0):raise ValueError('Zero-length wall graph edge')
    w=1/lengths;a,b=edges.T;n=len(points)
    lap=coo_matrix((np.r_[w,w,-w,-w],(np.r_[a,b,a,b],np.r_[a,b,b,a])),shape=(n,n)).tocsr()
    values=np.zeros(n);values[fixed]=fixed_values;free=np.setdiff1d(np.unique(edges),fixed)
    values[free]=spsolve(lap[free][:,free],-lap[free][:,fixed]@values[fixed])
    if not np.isfinite(values).all()or values.min()<-1e-12 or values.max()>1+1e-12:raise ValueError('Harmonic wall coordinate violated the maximum principle')
    residual=lap[free]@values
    return values,dict(free_controls=len(free),maximum_graph_residual=float(max(abs(residual),default=0)),positive_edge_weights=True)


def coordinates(points,tets,outer_vertices,inner_vertices):
    outer=np.unique(outer_vertices);inner=np.unique(inner_vertices)
    if not len(outer)or not len(inner)or np.intersect1d(outer,inner).size:raise ValueError('Distinct nonempty inner and outer carrier controls required')
    edges=np.concatenate([tets[:,[i,j]]for i,j in itertools.combinations(range(4),2)])
    fixed=np.r_[outer,inner];values,evidence=graph_values(points,edges,fixed,np.r_[np.zeros(len(outer)),np.ones(len(inner))])
    return values,dict(**evidence,outer_controls=len(outer),inner_controls=len(inner))


def layers(vertices,faces,outer_faces,inner_faces,names,*,tetrahedra=None,origin=None,free_source_vertices=(),rim_inner_vertices=(),method='distance-ratio',tetra_regions=None,interface_fractions=None):
    v=np.asarray(vertices,float);f=np.asarray(faces,int)
    fractions=np.arange(1,len(names),dtype=float)/len(names)if interface_fractions is None else np.asarray(interface_fractions,float)
    if not names or len(set(names))!=len(names)or fractions.shape!=(len(names)-1,)or not np.isfinite(fractions).all()or np.any(fractions<=0)or np.any(fractions>=1)or np.any(np.diff(fractions)<=0):raise ValueError('Unique coat names and strictly increasing interior wall fractions are required')
    if tetrahedra is None:p,t,used,index,evidence=tetrahedralize(v,f)
    else:
        p,t=tetrahedra;used=np.unique(f);index=np.arange(len(p));evidence={}
    # Boundary equality above establishes exact original control matching.
    lookup={tuple(q):i for i,q in enumerate(p)}
    free=set(map(int,free_source_vertices))
    inner_rim=set(map(int,rim_inner_vertices))
    if free&inner_rim or not inner_rim<=set(map(int,np.unique(inner_faces))):raise ValueError('Inner rim controls must belong to the lining and cannot also be free')
    if inner_rim and method not in ('harmonic','boundary-harmonic'):raise ValueError('Explicit inner-rim boundary roles require a harmonic wall coordinate')
    if method in ('harmonic','boundary-harmonic'):
        outer=[lookup[tuple(v[i])]for i in np.unique(outer_faces)if int(i)not in free|inner_rim];inner=[lookup[tuple(v[i])]for i in np.unique(inner_faces)if int(i)not in free]
        if method=='boundary-harmonic':
            boundary=np.array([[lookup[tuple(v[i])]for i in face]for face in f],int)
            edge=np.concatenate([boundary[:,[i,j]]for i,j in [(0,1),(1,2),(2,0)]])
            fixed=np.r_[outer,inner];surface,boundary_evidence=graph_values(p,edge,fixed,np.r_[np.zeros(len(outer)),np.ones(len(inner))])
            edge=np.concatenate([t[:,[i,j]]for i,j in itertools.combinations(range(4),2)])
            fixed=np.unique(boundary);values,harmonic=graph_values(p,edge,fixed,surface[fixed])
            harmonic.update(method=method,boundary_graph=boundary_evidence,all_boundary_values_fixed=True)
        else:values,harmonic=coordinates(p,t,outer,inner)
        field=np.c_[np.full(len(p),-1.),*[values-fraction for fraction in fractions]]
    elif method in ('distance-ratio','distance-fraction'):
        import point_cloud_utils as pcu
        do=pcu.closest_points_on_mesh(p,v,np.ascontiguousarray(outer_faces,dtype=np.int32))[0]
        di=pcu.closest_points_on_mesh(p,v,np.ascontiguousarray(inner_faces,dtype=np.int32))[0]
        outer=[lookup[tuple(v[i])]for i in np.unique(outer_faces)];inner=[lookup[tuple(v[i])]for i in np.unique(inner_faces)]
        do[outer]=0.;di[inner]=0.
        # Keep the two potentials separate. At a pinched attachment rim both
        # vanish and every coat terminates at the same original curve. A
        # single-valued 0/0 wall coordinate would manufacture a critical point.
        if method=='distance-fraction':
            denominator=do+di
            values=np.divide(do,denominator,out=np.full(len(p),.5),where=denominator>0)
            field=np.c_[np.full(len(p),-1.),*[values-fraction for fraction in fractions]]
        else:
            field=np.c_[np.full(len(p),-1.),*[(1-fraction)*do-fraction*di for fraction in fractions]]
        harmonic=dict(method=method,pinched_rim_controls=len(np.intersect1d(outer,inner)),distance_to_actual_triangles=True)
    else:raise ValueError('Unknown finite-wall coordinate method')
    center=np.asarray(origin if origin is not None else p.mean(0),float)
    def label(x):
        coat=names[int(np.sum(x[1:]>=0))]
        return str(tetra_regions[arrangement.current_tet_index])+':'+coat if tetra_regions is not None else coat
    arrangement=Arrangement(p-center,field,label)
    original_facets=np.array([[lookup[tuple(v[i])]for i in face]for face in f],int)
    arrangement.facet_labels={tuple(sorted(face)):i for i,face in enumerate(original_facets)}
    positions,surfaces,cells=arrangement.build(t,skip_fields=[0]);positions+=center
    return dict(positions=positions,surfaces=surfaces,cells=cells,source_faces=arrangement.surface_source_faces,evidence=dict(**evidence,harmonic=harmonic,interface_fractions=fractions.tolist(),inner_rim_controls_fixed_to_lining=len(inner_rim),source_wall_boundary_preserved=True,renderer_contains_tetrahedra=False))
