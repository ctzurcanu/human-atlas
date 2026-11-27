"""ONE corrected valve-local peripheral source-path proposal after v1 review.

Actual meridional source intersections define the outward radial envelope.
Paths advance monotonically through each cusp's commissure sector and are
restricted to its peripheral sheet, excluding the central coaptation floor.
This remains a qualified attachment proposal, not unique topological proof.
"""
from pathlib import Path
import sys,json,time,heapq,importlib.util
from collections import defaultdict
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.distance import barycentric

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-peripheral-source-hinges-v2b';D.mkdir(exist_ok=False)
previous=ROOT/'.local-models/geometry-workbench/heart/male-aortic-provisional-source-hinges-v1';prior=json.loads((previous/'evidence.json').read_text());center=np.array(prior['root_center_m']);axis=np.array(prior['root_axis']);commissures=prior['commissures'];meshes=[Mesh.load(previous/f'cusp-{i}.npz')for i in range(3)];u=np.cross(axis,[1,0,0]);u/=np.linalg.norm(u);v=np.cross(axis,u)
def polar(points):
    displacement=np.asarray(points)-center;height=displacement@axis;xy=np.stack([displacement@u,displacement@v],axis=-1);return np.arctan2(xy[...,1],xy[...,0]),np.linalg.norm(xy,axis=-1),height
def radial_envelope(mesh,angle):
    radial=np.cos(angle)*u+np.sin(angle)*v;normal=np.cross(axis,radial);values=(mesh.vertices-center)@normal;tri=mesh.faces;d=values[tri];edges=np.stack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]],axis=1);mask=values[edges[:,:,0]]*values[edges[:,:,1]]<0;faces=np.broadcast_to(np.arange(len(tri))[:,None],mask.shape)[mask];edges=edges[mask]
    fraction=values[edges[:,0]]/(values[edges[:,0]]-values[edges[:,1]]);points=mesh.vertices[edges[:,0]]+(mesh.vertices[edges[:,1]]-mesh.vertices[edges[:,0]])*fraction[:,None]
    # A meridional plane through an actual control can touch the outline
    # without crossing an edge. Retain that literal source control as an
    # exact vertex witness; the small bound is floating angular roundoff.
    touched=np.flatnonzero(abs(values)<1e-15)
    if len(touched):
        points=np.concatenate([points,mesh.vertices[touched]])
        faces=np.r_[faces,[int(np.argwhere(mesh.faces==k)[0,0])for k in touched]]
    radius=(points-center)@radial;best=int(np.argmax(radius));face=int(faces[best]);point=points[best];bary=barycentric(point[None],mesh.triangles[face:face+1])[0][0]
    return float(radius[best]),point,face,bary
def shortest(graph,start,end):
    best={start:0.};parent={};queue=[(0.,start)]
    while queue:
        cost,at=heapq.heappop(queue)
        if cost!=best.get(at):continue
        if at==end:
            path=[at]
            while at!=start:at=parent[at];path.append(at)
            return path[::-1]
        for nxt,weight in graph[at].items():
            trial=cost+weight
            if trial<best.get(nxt,np.inf):best[nxt]=trial;parent[nxt]=at;heapq.heappush(queue,(trial,nxt))
    raise ValueError('No sector-monotone actual source-edge path remains in the declared peripheral sheet')
def display(points,path,radius=.00013):
    vertices=[];faces=[]
    for aa,bb in zip(points[:-1],points[1:]):
        direction=bb-aa;direction/=np.linalg.norm(direction);s=np.cross(direction,np.eye(3)[np.argmin(abs(direction))]);s/=np.linalg.norm(s);t=np.cross(direction,s);ring=np.cos(np.arange(8)*np.pi/4)[:,None]*s+np.sin(np.arange(8)*np.pi/4)[:,None]*t;base=len(vertices);vertices.extend(aa+radius*ring);vertices.extend(bb+radius*ring)
        for k in range(8):
            j=(k+1)%8;faces.extend([[base+k,base+j,base+8+j],[base+k,base+8+j,base+8+k]])
    Mesh(np.array(vertices),np.array(faces,int).reshape(-1,3),meshes[0].frame).save(path)
records=[];entries=[];colors=[[.15,.76,.35,1],[.2,.48,.93,1],[.94,.68,.2,1]]
try:
 for i,mesh in enumerate(meshes):
    related=[r for r in commissures if i in r['cusps']];start,end=related;aa=polar(np.array(start['point_m']))[0];bb=polar(np.array(end['point_m']))[0];middle=polar(mesh.vertices.mean(0))[0]
    span=float((bb-aa)%(2*np.pi));inside=float((middle-aa)%(2*np.pi))
    if inside>span:start,end=end,start;aa,bb=bb,aa;span=float((bb-aa)%(2*np.pi))
    angle,radius,height=polar(mesh.vertices);phase=(angle-aa)%(2*np.pi);sector=(phase<=span+1e-12);envelope=np.empty(len(mesh.vertices));envelope_point=np.empty_like(mesh.vertices);envelope_face=np.empty(len(mesh.vertices),int);envelope_bary=np.empty_like(mesh.vertices)
    for k,theta in enumerate(angle):
        envelope[k],envelope_point[k],envelope_face[k],envelope_bary[k]=radial_envelope(mesh,float(theta))
    ratios=radius/np.maximum(envelope,1e-30)
    # Explicit geometric classification fraction, not a clinical thickness or
    # tolerance. We retain the discarded points and nearest-cusp distances.
    allowed=sector&(ratios>=.80)
    other_distance=[]
    for j,other in enumerate(meshes):
        if j!=i:other_distance.append(np.sqrt(igl.point_mesh_squared_distance(mesh.vertices,other.vertices,other.faces)[0]))
    near_other=np.min(other_distance,axis=0)
    graph=defaultdict(dict);edge_faces={}
    for fi,face in enumerate(mesh.faces):
        for va,vb in zip(face,np.roll(face,-1)):
            va,vb=int(va),int(vb)
            if not(allowed[va]and allowed[vb]):continue
            if phase[va]>phase[vb]:va,vb=vb,va
            if phase[vb]-phase[va]<1e-12:continue
            length=np.linalg.norm(mesh.vertices[va]-mesh.vertices[vb]);deficit=1-(ratios[va]+ratios[vb])/2;cost=length*(1+60*deficit*deficit)
            graph[va][vb]=cost;edge_faces[tuple(sorted((va,vb)))]=fi
    extras=[]
    for k,comm in enumerate([start,end]):
        face=comm['source_faces'][comm['cusps'].index(i)];point=np.array(comm['point_m']);node=len(mesh.vertices)+k;extras.append((node,face,point))
        for vertex in mesh.faces[face]:
            vertex=int(vertex)
            if not allowed[vertex]:continue
            length=np.linalg.norm(point-mesh.vertices[vertex])
            if k==0:graph[node][vertex]=length
            else:graph[vertex][node]=length
    path=shortest(graph,extras[0][0],extras[1][0]);points=[];source_faces=[];weights=[]
    for vertex in path:
        if vertex<len(mesh.vertices):fi,corner=map(int,np.argwhere(mesh.faces==vertex)[0]);point=mesh.vertices[vertex];bary=np.eye(3)[corner]
        else:_,fi,point=extras[vertex-len(mesh.vertices)];bary=barycentric(point[None],mesh.triangles[fi:fi+1])[0][0]
        points.append(point);source_faces.append(fi);weights.append(bary)
    points=np.array(points);source_faces=np.array(source_faces);weights=np.array(weights);residual=float(np.linalg.norm(np.einsum('ni,nij->nj',weights,mesh.triangles[source_faces])-points,axis=1).max());assert residual<1e-12
    angles,radii,heights=polar(points);phases=(angles-aa)%(2*np.pi);phases[0]=0;phases[-1]=span;assert np.all(np.diff(phases)>=-1e-10)
    file=D/f'cusp-{i}-qualified-arc.npz';atomic_npz(file,points=points,source_faces=source_faces,source_barycentric=weights,source_vertex_path=np.asarray(path),
             source_control_angles=angle,source_control_envelope_radius=envelope,source_control_envelope_points=envelope_point,source_control_envelope_faces=envelope_face,
             source_control_envelope_barycentric=envelope_bary,source_control_nearest_other_cusp_distance=near_other,source_control_allowed=allowed,frame=np.array(mesh.frame))
    visual=D/f'cusp-{i}-curve-display.npz';display(points,visual);mp=D/f'cusp-{i}.npz';mesh.save(mp);entries.extend([dict(mesh=str(mp),name=prior['arcs'][i]['cusp']+' unchanged source',color=colors[i]),dict(mesh=str(visual),name='Qualified peripheral source hinge proposal',color=[1,.12,.55,1])])
    records.append(dict(cusp=prior['arcs'][i]['cusp'],source_geometry_sha256=mesh.geometry_sha,witness=str(file),controls=len(points),
                        commissure_sector_degrees=float(span*180/np.pi),angular_backtracking_segments=0,minimum_radius_mm=float(radii.min()*1000),
                        local_height_range_mm=[float(heights.min()*1000),float(heights.max()*1000)],minimum_source_envelope_fraction=float(ratios[np.array(path[1:-1])].min()),
                        maximum_barycentric_residual_m=residual,attachment_status='Qualified source-constrained anatomical proposal; requires rendered review'))
 atomic_json(D/'review-spec.json',dict(frame=meshes[0].frame,resolution=1100,wall_seconds=50,meshes=entries,views=[dict(name='Qualified peripheral crown attachments — oblique',direction=[.4,.7,1]),dict(name='Qualified semilunar attachments — side',direction=[-1,.1,.4]),dict(name='Sector-monotone attachments — below',direction=(-axis).tolist())]))
 result=dict(schema='male-aortic-peripheral-source-hinge-arcs/v2',source_catalogue_sha256=prior['source_catalogue_sha256'],source_geometry_changed=False,source_reference_copied=False,
             method='Valve-local commissure sectors; actual meridional radial source envelope; one monotone peripheral source-edge path per cusp. No global Y extremum or hull.',
             clinical_tolerance_inferred=False,geometric_peripheral_envelope_fraction=.80,arcs=records,anatomical_acceptance=False,
             references=[dict(url='https://www.ahajournals.org/doi/pdf/10.1161/CIRCINTERVENTIONS.108.780858',scope='Parent-reviewed primary anatomical classification: crown attachments differ from virtual basal ring; no source coordinates copied',tool_access='403 during this subagent fetch'),dict(url='https://pubmed.ncbi.nlm.nih.gov/28160343/',scope='Parent-provided original CT reference; no coordinates copied')],elapsed_seconds=time.monotonic()-T)
 atomic_json(D/'evidence.json',result);print(json.dumps(result))
except Exception as error:
 atomic_json(D/'failure.json',dict(error=str(error),completed_arcs=records,source_geometry_changed=False,elapsed_seconds=time.monotonic()-T));raise
