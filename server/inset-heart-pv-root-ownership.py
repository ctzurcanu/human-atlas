"""Keep the original atrial lining at the root rather than sharing its exterior."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
from geometry_tools.audit import audit
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--inset-m',type=float,default=.000005)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False)
o=np.array([0.,1.375,.03125]);raw=h/'finite-pv-roots-v136';storage=h/'finite-pv-roots-v140'
fluid_id=json.loads((raw/'construction.json').read_text())['fluid_original_id']
fluid=Mesh.load(h/'blood-spaces-v130/fluid-boundary-4.npz');fs=md.Manifold(md.Mesh64(fluid.vertices-o,fluid.faces.astype('uint64')))
rows=[]
for i in range(4):
    original=Mesh.load(raw/f'root-{i}.npz');m=Mesh.load(storage/f'root-{i}.npz')
    z=np.load(raw/f'root-{i}-provenance.npz');faces=original.faces[z['original_id']==fluid_id]
    _,remap=np.unique((original.vertices-o).astype('<f4'),axis=0,return_inverse=True)
    assert remap.max()+1==len(m.vertices)
    q=original.vertices[faces];norm=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);sums=np.zeros_like(m.vertices)
    for k in range(3):np.add.at(sums,remap[faces[:,k]],norm)
    length=np.linalg.norm(sums,axis=1);active=length>0
    p=m.vertices.copy();p[active]-=args.inset_m*sums[active]/length[active,None]
    p=(p-o).astype('<f4').astype(float)+o
    candidate=Mesh(p,m.faces,m.frame);candidate.save(out/f'root-{i}.npz')
    proof=audit(candidate,require_closed=True,max_pairs=1000000,max_tests=7000000)
    atomic_json(out/f'root-{i}-audit.json',proof)
    body=md.Manifold(md.Mesh64(p-o,m.faces.astype('uint64')))
    overlap=(body^fs).volume()*1e6
    row=dict(i=i,passed=proof['passed']and overlap<1e-12,failures=proof['failures'],
        self_crossings=proof['self_intersections']['defect_pairs'],left_blood_occupied_overlap_mL=overlap,
        moved_controls=int(active.sum()),maximum_root_control_move_m=float(np.linalg.norm(p-m.vertices,axis=1).max()),
        atrial_source_moved=False,inset_m=args.inset_m)
    rows.append(row);print(json.dumps(row),flush=True)
atomic_json(out/'construction.json',dict(passed=all(r['passed']for r in rows),roots=rows,
    elapsed_seconds=time.monotonic()-start,installed=False))
