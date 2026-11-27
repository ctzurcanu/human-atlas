"""One source-witnessed provisional hinge candidate; no anatomical acceptance.

Outer pairwise cusp contacts provide commissure candidates. Each U-shaped
arc follows actual source edges through the cusp's basal extremum. Nothing
is copied from the separate native valve reference or installed as tissue.
"""
from pathlib import Path
import sys,json,time,heapq
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import exact_intersection_segment,barycentric

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-provisional-source-hinges-v1';D.mkdir(exist_ok=False)
source=ROOT/'.local-models/geometry-workbench/jobs/b9c8595514cdc1b8671263ec425dbf7bbad8f059c484638ae4b24211c0b2d9dc/prepared-2/proposal.json';a=json.loads(source.read_text());p={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';ids=['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet'];m=[reader.read_part(a,p[pid],frame=frame)for pid in ids]
centers=np.array([x.vertices.mean(0)for x in m]);center=centers.mean(0);axis=np.cross(centers[1]-centers[0],centers[2]-centers[0]);axis/=np.linalg.norm(axis)
if axis[1]<0:axis=-axis
commissures=[]
for i,j in [(0,1),(1,2),(2,0)]:
    contacts,summary=pair_contacts(m[i],m[j],0.,max_pairs=10000,max_tests=100000)
    assert summary['complete_within_pair_budget'];pairs=contacts['pairs'][contacts['relation']=='crossing'];options=[]
    for fa,fb in pairs:
        hit=exact_intersection_segment(m[i].triangles[fa],m[j].triangles[fb])
        for point in hit['endpoints']:
            displacement=point-center;height=displacement@axis;radius=np.linalg.norm(displacement-height*axis)
            options.append((radius,float(height),point,int(fa),int(fb)))
    # A peripheral, superior source contact is a candidate commissure. The
    # central contacts belong to coaptation and are explicitly excluded.
    cutoff=float(np.median([r[1]for r in options]));eligible=[r for r in options if r[1]>=cutoff];chosen=max(eligible,key=lambda r:r[0])
    radius,height,point,fa,fb=chosen
    first_bary=barycentric(point[None],m[i].triangles[fa:fa+1])[0][0];second_bary=barycentric(point[None],m[j].triangles[fb:fb+1])[0][0]
    commissures.append(dict(cusps=[i,j],point_m=point.tolist(),source_faces=[fa,fb],source_barycentric=[first_bary.tolist(),second_bary.tolist()],
                           radius_mm=radius*1000,height_mm=height*1000,classification='Provisional outer/superior actual pairwise cusp contact'))

def dijkstra(graph,start,end):
    distance={start:0.};previous={};queue=[(0.,start)]
    while queue:
        d,at=heapq.heappop(queue)
        if d!=distance.get(at):continue
        if at==end:
            path=[at]
            while at!=start:at=previous[at];path.append(at)
            return path[::-1]
        for nxt,cost in graph[at].items():
            nd=d+cost
            if nd<distance.get(nxt,np.inf):distance[nxt]=nd;previous[nxt]=at;heapq.heappush(queue,(nd,nxt))
    raise ValueError('No actual source-surface path joins the proposed sites')

def curve_display(points,path):
    vertices=[];faces=[]
    for start,end in zip(points[:-1],points[1:]):
        axis=end-start;length=np.linalg.norm(axis)
        if length<1e-12:continue
        axis/=length;u=np.cross(axis,np.eye(3)[np.argmin(abs(axis))]);u/=np.linalg.norm(u);v=np.cross(axis,u);rings=np.cos(np.arange(8)*np.pi/4)[:,None]*u+np.sin(np.arange(8)*np.pi/4)[:,None]*v;base=len(vertices)
        vertices.extend(start+.00013*rings);vertices.extend(end+.00013*rings)
        for k in range(8):
            j=(k+1)%8;faces.extend([[base+k,base+j,base+8+j],[base+k,base+8+j,base+8+k]])
    Mesh(np.array(vertices),np.array(faces,int).reshape(-1,3),frame).save(path)

arcs=[];entries=[];colors=[[.15,.76,.35,1],[.2,.48,.93,1],[.94,.68,.2,1]]
for i,mesh in enumerate(m):
    related=[r for r in commissures if i in r['cusps']];assert len(related)==2
    radial=centers[i]-center;radial-=axis*(radial@axis);radial/=np.linalg.norm(radial);position=(mesh.vertices-center)@radial;span=np.ptp(position)
    height=(mesh.vertices-center)@axis;basal=int(np.argmin(height));graph=defaultdict(dict);edge_face={}
    for fi,face in enumerate(mesh.faces):
        for aa,bb in zip(face,np.roll(face,-1)):
            aa,bb=int(aa),int(bb);length=np.linalg.norm(mesh.vertices[aa]-mesh.vertices[bb]);penalty=1+4*max(0.,float(position.max()-(position[aa]+position[bb])*.5))/span
            graph[aa][bb]=graph[bb][aa]=length*penalty;edge_face[tuple(sorted([aa,bb]))]=fi
    extra=[]
    for k,comm in enumerate(related):
        face=comm['source_faces'][comm['cusps'].index(i)];point=np.array(comm['point_m']);index=len(mesh.vertices)+k;extra.append((index,face,point))
        for vertex in mesh.faces[face]:
            vertex=int(vertex);distance=np.linalg.norm(point-mesh.vertices[vertex]);graph[index][vertex]=graph[vertex][index]=distance
    first=dijkstra(graph,extra[0][0],basal);last=dijkstra(graph,basal,extra[1][0]);path=first+last[1:];points=[];source_faces=[];weights=[]
    for vertex in path:
        if vertex<len(mesh.vertices):
            fi,corner=map(int,np.argwhere(mesh.faces==vertex)[0]);point=mesh.vertices[vertex];bary=np.eye(3)[corner]
        else:
            _,fi,point=extra[vertex-len(mesh.vertices)];bary=barycentric(point[None],mesh.triangles[fi:fi+1])[0][0]
        points.append(point);source_faces.append(fi);weights.append(bary)
    points=np.array(points);source_faces=np.array(source_faces);weights=np.array(weights)
    residual=float(np.linalg.norm(np.einsum('ni,nij->nj',weights,mesh.triangles[source_faces])-points,axis=1).max());assert residual<1e-12
    repeated=len(path)-len(set(path));witness=D/f'cusp-{i}-provisional-hinge.npz';atomic_npz(witness,points=points,source_faces=source_faces,source_barycentric=weights,source_vertex_path=np.asarray(path),frame=np.array(frame));display=D/f'cusp-{i}-diagnostic-curve.npz';curve_display(points,display)
    arcs.append(dict(cusp=ids[i],source_geometry_sha256=mesh.geometry_sha,witness_file=str(witness),source_controls=len(points),
                     basal_source_vertex=basal,basal_height_mm=float(height[basal]*1000),arc_length_mm=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum()*1000),
                     repeated_path_controls=repeated,maximum_barycentric_residual_m=residual,classification='Provisional source-constrained U-shaped basal path; render review required',confidence='provisional geometric proposal'))
    mp=D/f'cusp-{i}.npz';mesh.save(mp);entries.extend([dict(mesh=str(mp),name=ids[i]+' unchanged source',color=colors[i]),dict(mesh=str(display),name='PROVISIONAL basal hinge — actual source path',color=[1,.12,.55,1])])
atomic_json(D/'review-spec.json',dict(frame=frame,resolution=1100,wall_seconds=50,meshes=entries,views=[dict(name='Provisional source-constrained hinge arcs, oblique',direction=[.4,.7,1]),dict(name='Provisional basal arcs, side',direction=[-1,.1,.4]),dict(name='Provisional annular arcs from below',direction=(-axis).tolist())]))
report=dict(schema='male-aortic-provisional-source-hinge-arcs/v1',source_catalogue_sha256=digest_file(source),source_geometry_modified=False,
             root_axis=axis.tolist(),root_center_m=center.tolist(),commissures=commissures,arcs=arcs,
             native_valve_reference_copied=False,convex_hull_used=False,annulus_anatomically_accepted=False,
             limitation='Basal extremum and peripheral superior contact are explicit geometric hypotheses, not source-declared anatomical hinge ownership.',elapsed_seconds=time.monotonic()-T)
atomic_json(D/'evidence.json',report);print(json.dumps(report))
