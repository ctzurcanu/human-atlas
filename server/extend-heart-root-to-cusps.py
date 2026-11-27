"""Finite root continuation to the measured peripheral cusp envelope.

The source cusps stay in their native assembly. Their peripheral envelope
defines a provisional hinge crown, not a copied textbook fibrous annulus.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from scipy.spatial import ConvexHull
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,finite_wall_cells
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--valve',required=True);ap.add_argument('--trim-levels',required=True,type=int);ap.add_argument('--stage',required=True);args=ap.parse_args()
start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
src=h/args.input;old=json.loads((src/'construction.json').read_text());assert old['passed']
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());assert digest_file(path)==old['source_catalogue_sha256'];parts={p['id']:p for p in a['parts']};rd=AtlasReader()
z=np.load(src/'source-walls.npz');outer,inner,f,level=z['outer_positions'],z['inner_positions'],z['source_faces'],z['station_indices']
loops=boundary_loops(f);prox=min(loops,key=lambda ring:outer[ring].mean(0)[1]);pl=int(level[prox[0]]);direction=1 if pl==0 else -1;anchor=pl+direction*args.trim_levels
keep=level>=anchor if direction==1 else level<=anchor
ff=f[np.all(keep[f],axis=1)];ids=np.unique(ff);lookup=np.full(len(outer),-1,int);lookup[ids]=np.arange(len(ids));ff=lookup[ff];oo,ii=outer[ids],inner[ids]
other=outer[max(loops,key=lambda ring:outer[ring].mean(0)[1])].mean(0);c=outer[prox].mean(0)
cusps=next(g['partIds'] for g in a['finiteCardiacIntegration']['termGroups'] if g['key']==args.valve)
points=np.unique(np.concatenate([rd.read_part(a,parts[p],frame='body',catalogue_path=path).vertices for p in cusps]),axis=0)
# The source tube's cut is not the valve annulus. Use the cusp assembly's
# own plane and peripheral polygon centre, rather than forcing it concentric
# with that distant tube cut.
c=points.mean(0);_,_,vh=np.linalg.svd(points-c);normal=vh[-1]
if normal@(other-c)<0:normal=-normal
axis=np.array([0.,0.,1.]);axis-=normal*(axis@normal)
if np.linalg.norm(axis)<.1:axis=np.cross(normal,[0.,1.,0.])
axis/=np.linalg.norm(axis);side=np.cross(normal,axis);basis=np.stack([axis,side],axis=1)
xy=(points-c)@basis;hull=ConvexHull(xy);polygon=xy[hull.vertices]
cross=polygon[:,0]*np.roll(polygon[:,1],-1)-np.roll(polygon[:,0],-1)*polygon[:,1]
center=((polygon+np.roll(polygon,-1,axis=0))*cross[:,None]).sum(0)/(3*cross.sum())
c+=center@basis.T;xy=(points-c)@basis;heights=(points-c)@normal;hull=ConvexHull(xy);polygon=xy[hull.vertices];hh=heights[hull.vertices]
ring=min(boundary_loops(ff),key=lambda r:np.linalg.norm(oo[r].mean(0)-outer[level==anchor].mean(0)))
ac=oo[ring].mean(0);_,_,avh=np.linalg.svd(oo[ring]-ac);an=avh[-1]
if an@(other-ac)<0:an=-an
ax=oo[ring[0]]-ac;ax-=an*(ax@an);ax/=np.linalg.norm(ax);ay=np.cross(an,ax)
# Parameterise in the anchor's own plane, then transport the angular frame.
# Projecting around an offset crown centre can reverse angular order and fold
# an otherwise valid ring. No radial contour is resampled in that projection.
if normal@an<0:normal=-normal
crossn=np.cross(an,normal);dot=float(an@normal)
skew=np.array([[0.,-crossn[2],crossn[1]],[crossn[2],0.,-crossn[0]],[-crossn[1],crossn[0],0.]])
rotation=np.eye(3)+skew+skew@skew/(1+dot)
axis=rotation@ax;side=np.cross(normal,axis);basis=np.stack([axis,side],axis=1)
xy=(points-c)@basis;heights=(points-c)@normal;hull=ConvexHull(xy);polygon=xy[hull.vertices];hh=heights[hull.vertices]
rawangle=np.arctan2((oo[ring]-ac)@ay,(oo[ring]-ac)@ax)
steps=np.linalg.norm(np.roll(oo[ring],-1,axis=0)-oo[ring],axis=1)
turn=np.arctan2(np.sin(np.roll(rawangle,-1)-rawangle),np.cos(np.roll(rawangle,-1)-rawangle)).sum()
angle=rawangle[0]+np.sign(turn)*2*np.pi*np.r_[0.,np.cumsum(steps[:-1])]/steps.sum()
crown=[];radii=[]
for theta in angle:
 ray=np.array([np.cos(theta),np.sin(theta)]);found=[]
 for k,(p,q) in enumerate(zip(polygon,np.roll(polygon,-1,axis=0))):
  mat=np.stack([ray,p-q],axis=1)
  if abs(np.linalg.det(mat))<1e-20:continue
  radius,t=np.linalg.solve(mat,p)
  if radius>0 and -.00000001<=t<=1.00000001:found.append((radius,(1-t)*hh[k]+t*hh[(k+1)%len(hh)]))
 assert found,'No source cusp peripheral envelope in one angular sector'
 radius,height=min(found);radii.append(radius);crown.append(c+radius*(ray@basis.T)+height*normal)
crown=np.array(crown);radial=np.cos(angle)[:,None]*axis+np.sin(angle)[:,None]*side
co=crown+.00075*radial;ci=crown-.00045*radial
origlevel=level[ids];level_order=sorted(set(origlevel),reverse=direction==1);levels=np.array([level_order.index(v) for v in origlevel],int)
newlevels=[];previous=anchor+direction
oldring=oo[ring];oldinner=ii[ring]
# Correspondence follows literal native carrier indices at the retained ring.
source_ring=ids[ring];previous_ids=[]
for vertex in source_ring:
 neighbors=np.unique(f[np.any(f==vertex,axis=1)]);possible=neighbors[level[neighbors]==previous]
 assert len(possible),'No retained longitudinal source neighbour'
 previous_ids.append(int(min(possible,key=lambda v:np.linalg.norm(outer[v]-outer[vertex]))))
distance=float(np.linalg.norm(oldring.mean(0)-co.mean(0)))
d=oldring.mean(0)-outer[previous_ids].mean(0);t0=d/np.linalg.norm(d)*distance
# The peripheral cusp crown is scalloped. Its best-fit plane alone must not
# reverse the arriving wall course; use the actual centre displacement here.
t1=co.mean(0)-oldring.mean(0)
for k,t in enumerate(np.linspace(0,1,6)[1:]):
 def hermite(start,end):return (2*t**3-3*t*t+1)*start+(t**3-2*t*t+t)*t0+(-2*t**3+3*t*t)*end+(t**3-t*t)*t1
 no=hermite(oldring,co);ni=hermite(oldinner,ci);newring=np.arange(len(oo),len(oo)+len(ring));ff=np.r_[ff,annular_rim(ring,newring)];oo=np.r_[oo,no];ii=np.r_[ii,ni];levels=np.r_[levels,np.full(len(ring),len(level_order)+k)];ring=newring
fields=np.stack([oo,ii]);wall=finite_wall_cells(fields,ff);n=len(oo);full=[ff,ff[:,[0,2,1]]+n]
for loop in boundary_loops(ff):full.append(annular_rim(loop,loop+n))
mesh=Mesh(np.r_[oo,ii],np.concatenate(full),'model:local-ta98:assembled-metres')
checks=audit(mesh,require_closed=True,max_pairs=250000,max_tests=2000000)
passed=bool(wall['finite_wall_cells_pass'] and checks['passed'])
mesh.save(out/'diagnostic-root.npz')
atomic_npz(out/'diagnostic-walls.npz',outer_positions=oo,inner_positions=ii,source_faces=ff,station_indices=levels)
atomic_json(out/'construction.json',dict(passed=passed,source_part_id=old['source_part_id'],source_catalogue_sha256=digest_file(path),
 source_stage=args.input,cusp_part_ids=cusps,native_cusp_assembly_preserved=True,hinge_extent_provisional=True,
 root_continuation=True,illustrative_crown_wall_thickness_m=.0012,trimmed_levels=args.trim_levels,
 crown_radius_range_m=[min(radii),max(radii)],wall=wall,audit=checks,elapsed_seconds=time.monotonic()-start,installed=False))
if passed:
 mesh.save(out/'finite-root.npz');atomic_npz(out/'source-walls.npz',outer_positions=oo,inner_positions=ii,source_faces=ff,station_indices=levels,centers=np.array([oo[levels==j].mean(0) for j in range(levels.max()+1)]))
print(json.dumps(dict(passed=passed,finite_wall=wall['finite_wall_cells_pass'],failures=checks['failures'],crown_radius_mm=[min(radii)*1000,max(radii)*1000],elapsed_seconds=time.monotonic()-start)),flush=True)
raise SystemExit(0 if passed else 2)
