"""Variable grey/white tissue inside the literal repaired female cord.

The external source surface is constrained, not replaced with a cylinder.
The canal is a void. Grey/white territories share each interface exactly once.
Microscopic nuclei and tract extents are not established by this core stage.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.interpolate import CubicSpline
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology,winding_contains

BASE=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['domain','partition']);ap.add_argument('--output',required=True)
args=ap.parse_args();start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
surface=BASE/'female-surface-v1';src=np.load(surface/'surface.npz');uv,inv=np.unique(src['positions'],axis=0,return_inverse=True);uf=inv[src['faces']]
source=Mesh(uv,uf,'model:local-female-ta98:assembled-metres');assert topology(source)['closed_oriented_manifold']

def dimensions(y):
    t=source.triangles;hits=[];dd=t[:,:,1]-y
    for i,j in [(0,1),(1,2),(2,0)]:
        keep=dd[:,i]*dd[:,j]<0
        w=dd[keep,i]/(dd[keep,i]-dd[keep,j]);hits.extend(t[keep,i]+w[:,None]*(t[keep,j]-t[keep,i]))
    h=np.asarray(hits);lo=h.min(0);hi=h.max(0)
    return np.array([(lo[0]+hi[0])/2,(lo[2]+hi[2])/2,(hi[0]-lo[0])/2,(hi[2]-lo[2])/2])

# An original smooth macroscopic profile, not copied image coordinates.
# Positive local z is anterior. Thin dorsal and broader ventral horns join
# through the intermediate grey/commissural bridge.
right=np.array([[0.,.115],[.18,.15],[.25,.35],[.37,.53],[.52,.54],[.64,.41],
    [.61,.24],[.43,.065],[.41,-.13],[.43,-.33],[.50,-.54],[.57,-.77],
    [.40,-.63],[.27,-.39],[.20,-.17],[0.,-.115]])
# Mirror x only; keep anterior/posterior anatomy unchanged.
controls=np.r_[right,right[-2:0:-1]*[-1.,1.]]
controls=np.r_[controls,controls[:1]];curve=CubicSpline(np.arange(len(controls)),controls,bc_type='periodic')
profile=curve(np.linspace(0,len(controls)-1,80,endpoint=False))

def local_profile(y):
    # Different regional proportions: thoracic ventral horns are smaller;
    # grey fraction increases in the lumbosacral enlargement and conus.
    # These are macroscopic reference-guided fields, not histological borders.
    scale=np.interp(y,[1.148,1.158,1.178,1.198,1.24,1.42,1.47,1.53,1.565],
        [1.16,1.16,1.08,.95,.69,.69,.86,.86,.81])
    p=profile.copy();p[:,0]*=scale
    p[:,1]=np.where(p[:,1]>0,p[:,1]*scale,p[:,1])
    return p

def polygon_signed(points,polygon):
    a=polygon;b=np.roll(polygon,-1,axis=0);edge=b-a
    d=points[:,None,:]-a[None,:,:];den=np.einsum('ij,ij->i',edge,edge)
    tt=np.clip(np.einsum('nki,ki->nk',d,edge)/den,0,1)
    dist=np.linalg.norm(d-tt[:,:,None]*edge,axis=2).min(1)
    inside=np.zeros(len(points),bool)
    for u,w in zip(a,b):
        crosses=(u[1]>points[:,1])!=(w[1]>points[:,1])
        if u[1]!=w[1]:inside^=crosses & (points[:,0]<(w[0]-u[0])*(points[:,1]-u[1])/(w[1]-u[1])+u[0])
    return np.where(inside,-dist,dist)

if args.stage=='domain':
    import meshpy.tet as tet
    # Native source long edges are retained; interior samples improve the
    # compartment boundary without remeshing the external envelope.
    ys=np.linspace(uv[:,1].min()+.001,uv[:,1].max()-.001,64)
    dims=np.array([dimensions(y)for y in ys]);points=[]
    for y,(cx,cz,rx,rz)in zip(ys,dims):
        points.append([cx,y,cz])
        for radius in [.00016]:
            for angle in np.linspace(0,2*np.pi,8,endpoint=False):points.append([cx+radius*np.cos(angle),y,cz+radius*np.sin(angle)])
        for radius in [.28,.58,.82]:
            for angle in np.linspace(0,2*np.pi,16,endpoint=False):points.append([cx+rx*radius*np.cos(angle),y,cz+rz*radius*np.sin(angle)])
        for x,z in local_profile(y)[::4]:points.append([cx+x*rx,y,cz+z*rz])
    points=np.asarray(points);inside=winding_contains(source,points)['inside'];points=points[inside]
    p=np.r_[uv,points];origin=np.array([0.,1.35,-.05]);info=tet.MeshInfo();info.set_points(((p-origin)*1000).tolist());info.set_facets(uf.tolist())
    options=tet.Options('pQYY');options.nomergefacet=1;options.nomergevertex=1;options.nobisect=2
    result=tet.build(info,options=options);v=np.asarray(result.points)/1000+origin;t=np.asarray(result.elements,dtype='int64')
    assert np.max(abs(v[:len(p)]-p))<1e-12;v[:len(p)]=p
    q=v[t];det=np.linalg.det(q[:,1:]-q[:,:1]);assert np.all(det!=0);negative=det<0;t[negative]=t[negative][:,[0,2,1,3]]
    faces=np.concatenate([t[:,ii]for ii in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]]);keys,counts=np.unique(np.sort(faces,axis=1),axis=0,return_counts=True)
    boundary=keys[counts==1];expected=np.sort(uf,axis=1);expected=expected[np.lexsort(expected.T[::-1])];assert np.array_equal(boundary,expected),'External source facets must remain literal'
    actual=float(abs(det).sum()/6);reference=topology(source)['oriented_triangle_integral_m3'];assert abs(actual-reference)<reference*1e-8
    atomic_npz(out/'domain.npz',positions=v,tetrahedra=t,source_faces=uf,source_position_count=np.array(len(uv)),sample_y=ys,sample_dimensions=dims)
    atomic_json(out/'domain.json',dict(passed=True,source_surface_sha256=digest_file(surface/'surface.npz'),
        source_facets_literal=True,source_points_literal=True,points=len(v),tetrahedra=len(t),volume_mL=actual*1e6,
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
    print(json.dumps(dict(points=len(v),tetrahedra=len(t),volume_mL=actual*1e6,source_facets_literal=True,elapsed_seconds=time.monotonic()-start)),flush=True)
else:
    from geometry_tools.implicit_cells import Arrangement
    domain=BASE/'female-domain-v1';z=np.load(domain/'domain.npz');v=z['positions'];ys=z['sample_y'];dims=z['sample_dimensions']
    fields=np.zeros((len(v),3));fields[:,0]=-1.;profiles=[]
    # Evaluate the profile at each actual vertex level; no nearest-segment
    # discontinuity or one repeated thoracic profile through the whole cord.
    ylevels=np.unique(v[:,1])
    for y in ylevels:
        keep=v[:,1]==y;cx,cz,rx,rz=[np.interp(y,ys,dims[:,i])for i in range(4)]
        xy=np.c_[(v[keep,0]-cx)/rx,(v[keep,2]-cz)/rz]
        fields[keep,1]=polygon_signed(xy,local_profile(y))*min(rx,rz)
        # 0.08 mm canal radius: a boundary-owned void, never a tissue cell.
        fields[keep,2]=np.hypot(v[keep,0]-cx,v[keep,2]-cz)-.00008
    def classify(values):return None if values[2]<0 else 'grey'if values[1]<0 else'white'
    # Profile seeds intentionally lie on the same analytical boundary. Binary
    # evaluation can put them 1e-18 m either side, creating zero-thickness leaves.
    # Classify sub-picometre field residuals as on-boundary; never move controls.
    fields[np.abs(fields)<1e-12]=0.
    arrangement=Arrangement(v,fields,classify)
    positions,surfaces,cells=arrangement.build(z['tetrahedra'])
    arrays={'positions':positions}
    for sid,faces in surfaces.items():arrays[sid]=faces
    atomic_npz(out/'cells.npz',**arrays)
    atomic_json(out/'cells.json',dict(cells=cells,surfaces={k:len(v)for k,v in surfaces.items()},
        source_domain_sha256=digest_file(domain/'domain.npz'),source_envelope_preserved=True,
        boundary_classification='Reference-guided macroscopic grey/white fields; detailed laminae and nuclei not established.',
        canal_representation='Void, no solid plug',field_classification_tolerance_m=1e-12,
        positions=len(positions),elapsed_seconds=time.monotonic()-start,
        anatomical_acceptance=False,installed=False))
    print(json.dumps(dict(cells=list(cells),surfaces={k:len(v)for k,v in surfaces.items()},positions=len(positions),elapsed_seconds=time.monotonic()-start)),flush=True)
