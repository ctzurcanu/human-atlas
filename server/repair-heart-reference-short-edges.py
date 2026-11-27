"""Contract only diagnosed, short, new-reference edges with the link condition."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();base=ROOT/'.local-models/geometry-workbench/heart';src=base/'typed-conus-faithful-v663';out=base/args.stage;out.mkdir(exist_ok=False);started=time.monotonic();term='A12.1.02.008';meta=copy.deepcopy(json.loads((src/'construction.json').read_text()));g=meta['geometries'][term];m=Mesh.load(src/g['file']);v=np.array(m.vertices);f=np.array(m.faces);ev=json.loads((src/'float32-audit.json').read_text())['geometries'][term];assert ev['self_intersections']['complete_within_pair_budget'];affected=set(i for pair in ev['self_intersections']['examples']for i in pair);candidates=[]
for face in f[sorted(affected)]:
    for aa,bb in zip(face,np.roll(face,-1)):
        length=float(np.linalg.norm(v[aa]-v[bb]));
        if length<=.0001:candidates.append((length,int(min(aa,bb)),int(max(aa,bb))))
changes=[]
for length,aa,bb in sorted(set(candidates)):
    if any(aa in row['vertices']or bb in row['vertices']for row in changes):continue
    at_a=np.flatnonzero(np.any(f==aa,axis=1));at_b=np.flatnonzero(np.any(f==bb,axis=1));common=set(at_a)&set(at_b)
    if len(common)!=2:continue
    neighbors_a=set(f[at_a].ravel())-{aa,bb};neighbors_b=set(f[at_b].ravel())-{aa,bb};opposite=set(f[list(common)].ravel())-{aa,bb}
    if neighbors_a&neighbors_b!=opposite:continue
    mapped=np.where(f==bb,aa,f);keep=(mapped[:,0]!=mapped[:,1])&(mapped[:,0]!=mapped[:,2])&(mapped[:,1]!=mapped[:,2]);proposed=v.copy();proposed[aa]=((v[aa]+v[bb])/2).astype('float32').astype(float);incident=np.array(sorted((set(at_a)|set(at_b))-common));before=v[f[incident]];after=proposed[mapped[incident]];n=np.cross(before[:,1]-before[:,0],before[:,2]-before[:,0]);nn=np.cross(after[:,1]-after[:,0],after[:,2]-after[:,0])
    if np.any(np.einsum('ij,ij->i',n,nn)<=0):continue
    v=proposed;f=mapped[keep];changes.append(dict(vertices=[aa,bb],lengthMetres=length,linkCondition=True,allSurvivingIncidentNormalsPreserved=True))
assert changes,'No diagnosed short edge satisfies topology and orientation constraints'
m=Mesh(v,f,m.frame);m.save(out/g['file']);g['shortEdgeReferenceRepair']=dict(changes=changes,maximumPermittedEdgeMetres=.0001,originalAnatomyChanged=False);evidence=audit(m,require_closed=True,max_pairs=1000000,max_tests=4000000);atomic_json(out/'construction.json',meta);atomic_json(out/'float32-audit.json',dict(complete=True,geometries={term:evidence},elapsed_seconds=time.monotonic()-started));print('Conus',evidence['passed'],evidence['failures'],'contractions',changes,'elapsed',round(time.monotonic()-started,2),flush=True)
if not evidence['passed']:raise SystemExit(2)
