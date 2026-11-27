"""Exact-coordinate topology and actual triangle checks; no anatomy inferred."""
from collections import defaultdict,Counter
import numpy as np
from .model import Mesh
from .spatial import BVH
from .contact import pair_contacts,self_contact_summary


def topology(mesh):
    used=np.unique(mesh.faces);v,inv=np.unique(mesh.vertices[used],axis=0,return_inverse=True)
    remap=np.full(len(mesh.vertices),-1,int);remap[used]=inv;f=remap[mesh.faces]
    edges=defaultdict(list);adj=defaultdict(set);links=defaultdict(lambda:defaultdict(set))
    for i,tri in enumerate(f):
        for a,b,c in [tri,np.roll(tri,-1),np.roll(tri,-2)]:
            a,b,c=int(a),int(b),int(c);edges[tuple(sorted((a,b)))].append((i,a,b));adj[a].add(b);adj[b].add(a);links[a][b].add(c);links[a][c].add(b)
    component_of={};groups=[]
    for start in range(len(v)):
        if start in component_of:continue
        stack=[start];ids=[];ci=len(groups)
        while stack:
            x=stack.pop()
            if x in component_of:continue
            component_of[x]=ci;ids.append(x);stack.extend(adj[x])
        groups.append({'vertices':len(ids),'edges':0,'faces':0,'boundary_edges':0})
    for edge,owners in edges.items():
        group=groups[component_of[edge[0]]];group['edges']+=1;group['boundary_edges']+=len(owners)==1
    for tri in f:groups[component_of[int(tri[0])]]['faces']+=1
    badlinks=[]
    for vertex,graph in links.items():
        visited=set();stack=[next(iter(graph))]
        while stack:
            x=stack.pop()
            if x in visited:continue
            visited.add(x);stack.extend(graph[x]-visited)
        degrees=[len(a)for a in graph.values()]
        if len(visited)!=len(graph)or not(all(d==2 for d in degrees)or degrees.count(1)==2 and all(d in [1,2]for d in degrees)):badlinks.append(vertex)
    badedges=sum(len(owners)>2 for owners in edges.values())
    badwinding=sum(len(o)==2 and o[0][1:]==o[1][1:]for o in edges.values())
    for group in groups:
        group['euler']=group['vertices']-group['edges']+group['faces']
        group['closed_orientable_genus']=(2-group['euler'])//2 if not group['boundary_edges']and not badedges and not badwinding and not badlinks and group['euler']<=2 and group['euler']%2==0 else None
    t=mesh.triangles;areas=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    center=t.reshape(-1,3).mean(0)if len(t)else np.zeros(3);q=t-center
    integral=float(np.einsum('ij,ij->i',q[:,0],np.cross(q[:,1],q[:,2])).sum()/6)
    open_edges=sum(len(o)==1 for o in edges.values());duplicates=len(f)-len(set(map(tuple,np.sort(f,axis=1))))
    closed=bool(len(f)and not open_edges and not badedges and not badwinding and not badlinks and not duplicates and np.all(areas>0))
    return dict(triangles=len(f),used_vertices=len(v),components=groups,open_edges=open_edges,nonmanifold_edges=badedges,inconsistent_winding_edges=badwinding,
                invalid_vertex_links=badlinks,duplicate_faces=duplicates,zero_area_faces=int(np.sum(areas==0)),surface_area_m2=float(areas.sum()),
                oriented_triangle_integral_m3=integral,closed_oriented_manifold=closed,exact_coordinate_welding_only=True)


def winding_contains(mesh,points,*,boundary_tolerance_m=0.,distance_backend='bvh'):
    """Full solid-angle sums, valid only for a closed oriented manifold.

    Opposite oriented nested shells describe cavities. Points on/near the
    surface are explicitly ambiguous, never silently classified as inside.
    """
    topo=topology(mesh)
    if not topo['closed_oriented_manifold']:raise ValueError('Containment requires a closed oriented manifold; an open sheet cannot define a compartment')
    points=np.asarray(points,float).reshape(-1,3);bvh=BVH(mesh);inside=[];ambiguous=[];winding=[]
    nearest,distance_evidence=bvh.closest_points(points,backend=distance_backend)
    for p,distance in zip(points,nearest['distance_m']):
        q=mesh.triangles-p
        scale=float(np.max(abs(q)));roundoff=256*np.finfo(float).eps*max(scale,np.finfo(float).tiny)
        if distance<=max(boundary_tolerance_m,roundoff):inside.append(False);ambiguous.append(True);winding.append(None);continue
        u,v,w=q[:,0],q[:,1],q[:,2];a,b,c=np.linalg.norm(u,axis=1),np.linalg.norm(v,axis=1),np.linalg.norm(w,axis=1)
        num=np.einsum('ij,ij->i',u,np.cross(v,w));den=a*b*c+np.einsum('ij,ij->i',u,v)*c+np.einsum('ij,ij->i',v,w)*a+np.einsum('ij,ij->i',w,u)*b
        number=float(np.sum(2*np.arctan2(num,den))/(4*np.pi));winding.append(number);inside.append(abs(number)>.5);ambiguous.append(abs(abs(number)-.5)<1e-8)
    return dict(inside=inside,boundary_ambiguous=ambiguous,winding_number=winding,all_triangles_used=True,distance_query=distance_evidence)


def audit(mesh,*,self_intersections=True,max_pairs=250000,max_tests=2000000,expected_components=None,expected_genus=None,require_closed=False):
    from .geometry_instance import instance_preimage
    instance=instance_preimage(mesh)
    if instance is not None:
        canonical,pose=instance
        result=audit(canonical,self_intersections=self_intersections,max_pairs=max_pairs,max_tests=max_tests,expected_components=expected_components,expected_genus=expected_genus,require_closed=require_closed)
        result.update(geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,frame=mesh.frame,geometry_instance=pose)
        result['topology']['surface_area_m2']*=pose['scale']**2
        result['topology']['oriented_triangle_integral_m3']*=pose['scale']**3
        if self_intersections:
            for key in ['first_candidate_patch_area_upper_bound_m2','second_candidate_patch_area_upper_bound_m2']:
                result['self_intersections'][key]*=pose['scale']**2
        return result
    topo=topology(mesh);out=dict(geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,frame=mesh.frame,topology=topo,anatomical_acceptance=False)
    failures=[]
    for key in ['nonmanifold_edges','inconsistent_winding_edges','invalid_vertex_links','duplicate_faces','zero_area_faces']:
        if topo[key]:failures.append(key)
    if require_closed and not topo['closed_oriented_manifold']:failures.append('not_closed_oriented_manifold')
    if expected_components is not None and len(topo['components'])!=expected_components:failures.append('unexpected_components')
    if expected_genus is not None and [g['closed_orientable_genus']for g in topo['components']]!=expected_genus:failures.append('unexpected_genus')
    if self_intersections:
        summary,_=self_contact_summary(mesh,max_pairs=max_pairs,max_tests=max_tests)
        out['self_intersections']={**summary,'triangle_test_budget':max_tests,'pair_result_budget':max_pairs}
        if not summary['complete_within_pair_budget']or summary['defect_pairs']:failures.append('self_intersection_check')
    else:out['self_intersections']={'checked':False}
    out.update(passed=not failures,failures=failures,scope='Numerical geometry gates only; reference, renders, sections and anatomical review remain separate')
    return out


def add_parser(sub):
    p=sub.add_parser('audit');p.add_argument('mesh');p.add_argument('--model');p.add_argument('--frame');p.add_argument('--skip-self',action='store_true');p.add_argument('--max-pairs',type=int,default=250000);p.add_argument('--max-tests',type=int,default=2000000);p.add_argument('--components',type=int);p.add_argument('--genus',type=int,nargs='+');p.add_argument('--closed',action='store_true');p.add_argument('--open3d',action='store_true',help='Include independent native C++ diagnostics');p.set_defaults(handler=run)


def run(args,db):
    from .cli import load_ref
    row,mesh=load_ref(db,args.mesh,args.model,args.frame)
    out=audit(mesh,self_intersections=not args.skip_self,max_pairs=args.max_pairs,max_tests=getattr(args,'max_tests',2000000),expected_components=args.components,expected_genus=args.genus,require_closed=args.closed)
    if args.open3d:
        from .open3d_backend import diagnostics
        out['open3d_diagnostics']=diagnostics(mesh,self_intersections=not args.skip_self)
    return out
