"""Bounded shared RBF field, constrained by actual liver surface witnesses.

One field transports all four wall networks and their reciprocal source views.
Hilar controls and the caval trunk are fixed. This script saves a candidate,
never accepts or installs it merely because sampled containment improves.
"""
from pathlib import Path
import json,sys,time,argparse
import numpy as np
import scipy.sparse as sp
from scipy.spatial.distance import cdist
import clarabel
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.open3d_backend import library

ap=argparse.ArgumentParser();ap.add_argument('--version',default='v161');ap.add_argument('--gradient-limit',type=float,default=.30);ap.add_argument('--surface-spacing',type=float,default=0.);ap.add_argument('--guard-spacing',type=float,default=.006);ap.add_argument('--extra-probes');args=ap.parse_args()
started=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/('hepatic-contained-course-'+args.version);out.mkdir(exist_ok=False)
a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']};o3d=library()
def caster(m):
    s=o3d.t.geometry.RaycastingScene(nthreads=1);s.add_triangles(o3d.core.Tensor(m.vertices.astype(np.float32)),o3d.core.Tensor(m.faces.astype(np.uint32)));return s
body=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');scene=caster(body);cava=r.read_part(a,parts['TA98HEPATIC:male:region:A12.3.09.001'],frame=body.frame);cv=caster(cava)
anchors=json.loads((base/'separated-courses-v105/evidence.json').read_text())['anchors']
def smooth(t):t=np.clip(t,0,1);return t*t*(3-2*t)
def weights(p):
    distance=cv.compute_distance(o3d.core.Tensor(p.astype(np.float32)),nthreads=1).numpy().astype(float)
    w=smooth((distance-.003)/.009)
    for anchor in anchors.values():w*=smooth((np.linalg.norm(p-np.asarray(anchor['center_m']),axis=1)-anchor['radius_m'])/.009)
    return w
old=json.loads((base/'vessel-containment-v157/evidence.json').read_text());sources={};probes=[]
for row in old['measuredIntrahepaticBranches']:
    p=parts[row['id']];m=r.read_part(a,p,frame=body.frame);sources[row['id']]=m
    probes.extend([m.vertices[np.unique(m.faces)],m.triangles.mean(1)])
if args.extra_probes:probes.append(np.load(args.extra_probes)['points'])
p=np.unique(np.concatenate(probes),axis=0);w=weights(p);signed=scene.compute_signed_distance(o3d.core.Tensor(p.astype(np.float32)),nthreads=1,nsamples=7).numpy().astype(float)
hit=scene.compute_closest_points(o3d.core.Tensor(p.astype(np.float32)),nthreads=1);normal=hit['primitive_normals'].numpy().astype(float)
eligible=(w>.15)&(signed>-.018);points=p[eligible];n=normal[eligible];d=signed[eligible];w=w[eligible]
outside=points[d>0];grid=np.floor(outside/.008).astype(int);_,groups=np.unique(grid,axis=0,return_inverse=True);controls=np.array([outside[groups==i].mean(0)for i in range(groups.max()+1)])
# Deterministic farthest controls cap numerical cost rather than growing a
# new per-vertex deformation problem.
if len(controls)>64:
    chosen=[0];nearest=np.linalg.norm(controls-controls[0],axis=1)
    while len(chosen)<64:
        i=int(nearest.argmax());chosen.append(i);nearest=np.minimum(nearest,np.linalg.norm(controls-controls[i],axis=1))
    controls=controls[chosen]
sigma=.018
phi=np.exp(-.5*cdist(points,controls,'sqeuclidean')/sigma**2)*w[:,None]
# Keep every outside witness and a 2 mm spatial sample of nearby inside
# witnesses. The complete bank is independently remeasured after fitting.
grid=np.floor(points/.002).astype(int);_,first=np.unique(grid,axis=0,return_index=True);take=np.unique(np.r_[np.flatnonzero(d>-.0015),first])
if args.surface_spacing:
    grid=np.floor(points/args.surface_spacing).astype(int);_,groups=np.unique(grid,axis=0,return_inverse=True)
    order=np.lexsort((-d,groups));_,unique=np.unique(groups[order],return_index=True);take=order[unique]
phi=phi[take];n=n[take];d=d[take];count=len(controls)
A=(phi[:,:,None]*n[:,None,:]).reshape(len(take),count*3)
target=-d-.00035
grid=np.floor(p/args.guard_spacing).astype(int);_,sample=np.unique(grid,axis=0,return_index=True);guards=p[sample]
def basis(x):return np.exp(-.5*cdist(x,controls,'sqeuclidean')/sigma**2)*weights(x)[:,None]
guard_basis=basis(guards);pieces=[sp.csc_matrix(A)];bounds=[target]
# Bound all three displacement components and the deformation gradient at
# spatial controls, including deep interior points omitted by surface planes.
for component in range(3):
    unit=np.eye(3)[component:component+1];B=sp.kron(sp.csc_matrix(guard_basis),sp.csc_matrix(unit),format='csc');pieces.extend([B,-B]);bounds.extend([np.full(len(guards),.015),np.full(len(guards),.015)])
for axis in range(3):
    aa=guards.copy();bb=guards.copy();aa[:,axis]+=1e-5;bb[:,axis]-=1e-5;G=(basis(aa)-basis(bb))/2e-5
    for component in range(3):
        B=sp.kron(sp.csc_matrix(G),sp.csc_matrix(np.eye(3)[component:component+1]),format='csc');pieces.extend([B,-B]);bounds.extend([np.full(len(guards),args.gradient_limit),np.full(len(guards),args.gradient_limit)])
matrix=sp.vstack(pieces,format='csc');target=np.concatenate(bounds)
# RKHS displacement energy plus positive coefficient regularization.
K=np.exp(-.5*cdist(controls,controls,'sqeuclidean')/sigma**2);P=sp.triu(sp.kron(K+np.eye(count)*.015,sp.eye(3),format='csc'),format='csc')
settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=10.;settings.max_iter=100;settings.max_threads=1;settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-9
solver=clarabel.DefaultSolver(P,np.zeros(count*3),matrix,target,[clarabel.NonnegativeConeT(len(target))],settings);result=solver.solve()
record=dict(status=str(result.status),solverSeconds=result.solve_time,controls=count,constraints=len(target),sigma=sigma,elapsedSeconds=time.monotonic()-started)
atomic_json(out/'solver.json',record)
if str(result.status)!='Solved':print(json.dumps(record),flush=True);sys.exit(1)
coeff=np.asarray(result.x).reshape(count,3);atomic_npz(out/'field.npz',controls=controls,coefficients=coeff,sigma=np.array(sigma))
def mapped(v):
    return v+weights(v)[:,None]*(np.exp(-.5*cdist(v,controls,'sqeuclidean')/sigma**2)@coeff)
rows=[]
for role in ['artery','portal','hepatic_vein','biliary']:
    m=r.read_part(a,parts['TA98HEPATIC:male:'+role],frame=body.frame);v=mapped(m.vertices);Mesh(v,m.faces,m.frame,m.attributes).save(out/(role+'.npz'))
    print(role,'maximum movement',float(np.linalg.norm(v-m.vertices,axis=1).max()),flush=True)
for row in old['measuredIntrahepaticBranches']:
    m=sources[row['id']];used=np.unique(m.faces);v=mapped(m.vertices[used]);signed=scene.compute_signed_distance(o3d.core.Tensor(v.astype(np.float32)),nthreads=1,nsamples=7).numpy().astype(float);weight=weights(m.vertices[used]);bad=(signed>3e-7)&(weight>.15)
    q=dict(term=row['term'],name=row['name'],unprotectedOutsideVertices=int(bad.sum()),maximumUnprotectedOutsideMetres=float(max(0,signed[weight>.15].max(initial=0))),
        protectedHilarOutsideVertices=int(((signed>3e-7)&(weight<=.15)).sum()),originalMaximumOutsideMetres=row['maximum_outside_m'])
    rows.append(q);print(json.dumps(q),flush=True)
atomic_json(out/'probe.json',dict(**record,branches=rows,geometryAccepted=False,installed=False,
    coordinateGuardMetres=3e-7,marginParameterMetres=.00035,allFourNetworksUseOneField=True,existingCavalAndHilarCoordinatesFixed=True))
