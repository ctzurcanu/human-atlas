"""One conforming transverse template for cord, canal, grey/white and pia.

The source curves come from actual fitted tissue sections. A single planar
triangulation owns all compartments. Interior centroids ensure that every
calculation prism has an unambiguous interior material control; interfaces
can carry exactly tied fields without duplicating positive tissue volume.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from scipy.interpolate import CubicSpline
import shapely
from shapely.geometry import Polygon,LineString,Point
from shapely.ops import unary_union,polygonize
import triangle
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();started=time.monotonic()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';source=base/'female-cage-rounded-v3'
d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');v=z['positions'];inner=set();outer=set()
for name,uses in d['cells'].items():
    if name.split('@')[0]!='pia':continue
    for use in uses:
        peer=use.get('neighbor')
        if peer and not peer.startswith('pia'):inner.add(use['surface'])
        elif not peer:outer.add(use['surface'])
triangles={kind:v[np.concatenate([z[s]for s in sorted(sids)])]for kind,sids in [('inner',inner),('outer',outer)]}
def cut_segments(t,y):
    yy=t[:,:,1];t=t[(yy.min(1)<=y)&(yy.max(1)>=y)];result=[]
    for q in t:
        delta=q[:,1]-y;hits=[row[[0,2]]for row in q[delta==0]]
        for a,b in [(0,1),(1,2),(2,0)]:
            if delta[a]*delta[b]<0:hits.append((q[a]+(q[b]-q[a])*delta[a]/(delta[a]-delta[b]))[[0,2]])
        hits=np.unique(np.asarray(hits).reshape(-1,2),axis=0)
        if len(hits)==2:result.append(hits)
    if not result:raise ValueError('No native cord contour at declared section')
    return np.asarray(result)

def radial_section(y,angles):
    segments={k:cut_segments(t,y)for k,t in triangles.items()};p=segments['inner'].reshape(-1,2)
    low,high=p.min(0),p.max(0);center=(low+high)/2;radii=(high-low)/2
    curves={};cross=lambda a,b:a[...,0]*b[...,1]-a[...,1]*b[...,0]
    for kind,s in segments.items():
        a=s[:,0]-center;b=s[:,1]-s[:,0];points=[]
        for theta in angles:
            direction=np.array([np.cos(theta),np.sin(theta)]);den=cross(direction,b)
            good=abs(den)>1e-20
            distance=np.divide(cross(a,b),den,out=np.full(len(a),np.nan),where=good)
            fraction=np.divide(cross(a,direction),den,out=np.full(len(a),np.nan),where=good)
            valid=good&(distance>0)&(fraction>=-1e-10)&(fraction<=1+1e-10)
            if not valid.any():raise ValueError('Cord contour is not star-shaped around the declared centre')
            points.append(center+distance[valid].min()*direction)
        curves[kind]=np.asarray(points)
    assert np.all(np.linalg.norm(curves['outer']-center,axis=1)>np.linalg.norm(curves['inner']-center,axis=1))
    return center,radii,curves

angles=np.linspace(0,2*np.pi,96,endpoint=False);y=1.516826123
center,radii,world=radial_section(y,angles)
curves={k:(p-center)/radii for k,p in world.items()}
right=np.array([[0.,.115],[.18,.15],[.25,.35],[.37,.53],[.52,.54],[.64,.41],
    [.61,.24],[.43,.065],[.41,-.13],[.43,-.33],[.50,-.54],[.57,-.77],
    [.40,-.63],[.27,-.39],[.20,-.17],[0.,-.115]])
controls=np.r_[right,right[-2:0:-1]*[-1.,1.],right[:1]];curve=CubicSpline(np.arange(len(controls)),controls,bc_type='periodic')
grey=curve(np.linspace(0,len(controls)-1,64,endpoint=False));grey[:,0]*=.86;grey[:,1]=np.where(grey[:,1]>0,grey[:,1]*.86,grey[:,1])
canal=np.c_[np.cos(np.linspace(0,2*np.pi,16,endpoint=False)),np.sin(np.linspace(0,2*np.pi,16,endpoint=False))]*.00008/radii
curves.update(grey=grey,canal=canal)
# Declare a reference arithmetic grid before forming the planar graph. Its
# half-step is below 3 pm at the reference physical scale. This removes
# numerical pseudo-regions around sin(pi), shared axes and line intersections.
grid=2.**-30;curves={k:np.rint(q/grid)*grid for k,q in curves.items()}
polys={k:Polygon(q)for k,q in curves.items()}
assert all(p.is_valid for p in polys.values())
assert polys['outer'].contains(polys['inner'])and polys['inner'].contains(polys['grey'])and polys['grey'].contains(polys['canal'])
lines=[p.boundary for p in polys.values()]
xs=np.array([-1.5,0.,1.5])
for height in [.15,-.18]:lines.append(LineString([(-1.5,height),(1.5,height)]).intersection(polys['inner']))
for intercept,slope in [(.13,.9),(-.16,-.9)]:lines.append(LineString(np.c_[xs,intercept+slope*abs(xs)]).intersection(polys['inner']))
lines.append(LineString([(0.,-1.5),(0.,1.5)]).intersection(polys['outer']))
regions=[];hole_seeds=[]
for p in polygonize(shapely.union_all([shapely.set_precision(line,grid)for line in lines],grid_size=grid)):
    c=p.representative_point()
    if not polys['outer'].covers(c):continue
    if polys['canal'].covers(c):hole_seeds.append([c.x,c.y]);continue
    x,zz=c.x,c.y;side='left'if x>0 else'right'
    if not polys['inner'].covers(c):name='pia'
    elif polys['grey'].covers(c):name='grey-'+('anterior'if zz>.15 else'posterior'if zz<-.18 else'intermediate')+'-'+side
    else:name='white-'+('anterior'if zz>.13+.9*abs(x)else'posterior'if zz<-.16-.9*abs(x)else'lateral')+'-'+side
    regions.append((name,p))
names=sorted({n for n,p in regions});assert len(names)==13,names
points=[];lookup={};edges=set();seeds=[]
def node(q):
    key=tuple(map(float,q))
    if key not in lookup:lookup[key]=len(points);points.append(key)
    return lookup[key]
for name,p in regions:
    for ring in [p.exterior,*p.interiors]:
        coords=list(ring.coords);ids=[node(q)for q in coords[:-1]]
        for a,b in zip(ids,ids[1:]+ids[:1]):edges.add(tuple(sorted((a,b))))
    c=p.representative_point();seeds.append([c.x,c.y,names.index(name),.012])
assert len(hole_seeds)==2,'The median constraint divides the canal into two void regions'
result=triangle.triangulate(dict(vertices=np.asarray(points),segments=np.array(sorted(edges),int),holes=hole_seeds,regions=np.array(seeds)), 'pqa0.012A')
uv=np.asarray(result['vertices'],float);faces=np.asarray(result['triangles'],int);labels=np.rint(result['triangle_attributes'].ravel()).astype(int)
area=lambda q:(q[:,1,0]-q[:,0,0])*(q[:,2,1]-q[:,0,1])-(q[:,1,1]-q[:,0,1])*(q[:,2,0]-q[:,0,0])
negative=area(uv[faces])<0;faces[negative]=faces[negative][:,[0,2,1]]
assert np.all(area(uv[faces])>0)
# Exact material/physical boundary controls carry tied zero fields. A new
# interior point in each face guarantees positive-volume material ownership.
membership=defaultdict(set)
for f,label in zip(faces,labels):
    for control in f:membership[int(control)].add(int(label))
edges0=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]);edge0,count0=np.unique(np.sort(edges0,axis=1),axis=0,return_counts=True)
boundary=set(map(int,np.unique(edge0[count0==1])))
interior_label=np.array([next(iter(membership[i]))if len(membership[i])==1 and i not in boundary else-1 for i in range(len(uv))])
first=len(uv);centroids=uv[faces].mean(1);expanded=[]
for k,f in enumerate(faces):
    for a,b in zip(f,np.roll(f,-1)):expanded.append([a,b,first+k])
uv=np.r_[uv,centroids];faces=np.asarray(expanded,int);interior_label=np.r_[interior_label,labels];labels=np.repeat(labels,3)
assert np.all(area(uv[faces])>0),'A centroid split must retain strictly positive reference area'
kind=np.full(len(uv),-1,int);edge_index=np.full(len(uv),-1,int);fraction=np.zeros(len(uv));residual=np.full(len(uv),np.inf)
curve_names=list(curves)
for number,name in enumerate(curve_names):
    poly=curves[name];a=poly;b=np.roll(poly,-1,axis=0);e=b-a;delta=uv[:,None]-a
    t=np.clip(np.einsum('nki,ki->nk',delta,e)/np.einsum('ki,ki->k',e,e),0,1)
    distance=np.linalg.norm(delta-t[:,:,None]*e,axis=2);nearest=distance.argmin(1);dist=distance[np.arange(len(uv)),nearest]
    keep=(dist<1e-9)&(dist<residual);kind[keep]=number;edge_index[keep]=nearest[keep];fraction[keep]=t[np.flatnonzero(keep),nearest[keep]];residual[keep]=dist[keep]
fixed=np.flatnonzero(kind>=0);assert set(boundary)<=set(map(int,fixed))
assert all(any(interior_label[j]==label for j in f)for f,label in zip(faces,labels))
atomic_npz(out/'template.npz',positions=uv,faces=faces,material=labels,interior_material=interior_label,angles=angles,
    fixed_controls=fixed,curve_kind=kind,curve_edge=edge_index,curve_fraction=fraction,reference_center=center,reference_radii=radii,
    **{'curve:'+k:p for k,p in curves.items()})
atomic_json(out/'template.json',dict(passed=True,materials=names,curves=curve_names,positions=len(uv),triangles=len(faces),
    boundary_controls=len(fixed),reference_y_m=y,source_tissues_sha256=digest_file(source/'cells.npz'),
    canal_is_void=True,canal_void_seeds=hole_seeds,interfaces_shared=True,reference_arithmetic_grid=grid,
    maximum_curve_reference_rounding_m=float(grid*.5*np.max(radii)),
    profile_source='Existing fitted cord and reference-guided macroscopic grey profile',
    no_microscopic_boundaries_implied=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(passed=True,materials=13,vertices=len(uv),triangles=len(faces),seconds=time.monotonic()-started)),flush=True)
