"""Prepare actual finite-wall/context renders without installing a catalogue."""
from pathlib import Path
import argparse, hashlib, json, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz

t=time.monotonic()
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1'
D=B/'localized-no-rejoin-wall-v1'
parser=argparse.ArgumentParser();parser.add_argument('--output-name',default='render-review')
args=parser.parse_args();assert Path(args.output_name).name==args.output_name
O=D/args.output_name;O.mkdir(exist_ok=False)
ap=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
a=json.loads(ap.read_text());reader=AtlasReader()
assert hashlib.sha256(ap.read_bytes()).hexdigest()=='668c0303671b46bb6c7b9d86e6b7461181d3f815b77cc031bb7cefd255642ac3'
p=D/'finite-wall-source-diagnostic.npz';m=Mesh.load(p)
with np.load(p,allow_pickle=False)as z: fields=z['corner_fields'].copy()
normals=fields[:,:,:3]*fields[:,:,3,None]
normal_length=np.linalg.norm(normals,axis=2)
assert np.isfinite(normals).all()and normal_length.min()>0
gn=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);gn/=np.linalg.norm(gn,axis=1)[:,None]
dots=np.einsum('fci,fi->fc',normals/normal_length[:,:,None],gn)
atomic_json(O/'actual-corner-normal-check.json',dict(nonforward_face_ids=np.flatnonzero(np.any(dots<=0,axis=1)).tolist(),minimum_dot=float(dots.min()),geometry_sha256=m.geometry_sha))
lo=np.array([.012,1.579,-.024]);hi=np.array([.065,1.650,.024])
def crop(mesh):
    tr=mesh.triangles
    keep=np.all(tr.max(1)>=lo,axis=1)&np.all(tr.min(1)<=hi,axis=1)
    return keep,Mesh(mesh.vertices,mesh.faces[keep],mesh.frame,mesh.attributes,
                     dict(mesh.source,diagnostic_uncapped_crop=True))
keep,c=crop(m);cp=O/'actual-wall-uncapped-crop.npz';c.save(cp)
nc=O/'actual-wall-corner-normals.npz';atomic_npz(nc,corner_normals=normals[keep])
meshes=[dict(mesh=str(cp),name='Actual communicating wall candidate; uncapped diagnostic crop',color=[.78,.16,.12,1],
             corner_normal_source=dict(path=str(nc),sha256=hashlib.sha256(nc.read_bytes()).hexdigest()))]
context=[]
for part in a['parts']:
    pid=part['id']
    chosen=pid in ['ZA:Glossopharyngeal nerve (IX).l','ZA:Vagus nerve (X).l','ZA:Internal jugular vein.l'] or pid.startswith('ZA:Sternocleidomastoid muscle.l:')
    if not chosen or part.get('suppressed')or part.get('defaultHidden')or pid in a.get('defaultHidden',[]):continue
    raw=reader.read_part(a,part,frame=m.frame,catalogue_path=ap);_,cut=crop(raw)
    if not len(cut.faces):continue
    path=O/('context-'+str(len(context))+'.npz');cut.save(path)
    color=[.82,.73,.28,.55]if part['system']=='nervous'else[.22,.35,.72,.35]if part['system']=='venous'else[.58,.29,.30,.16]
    meshes.append(dict(mesh=str(path),name='Fixed '+part['name']+'; diagnostic crop',color=color))
    context.append(dict(id=pid,source_geometry_sha256=raw.geometry_sha,source_attributes_sha256=raw.attributes_sha,cropped_render_only=True))
base=dict(frame=m.frame,meshes=meshes,resolution=1200,wall_seconds=40)
base['views']=[dict(name='candidate-fork-posterior-oblique',direction=[1,.1,-1]),
               dict(name='candidate-fork-longitudinal-cut',direction=[1,0,0],section=dict(normal=[1,0,0],point=[.033,1.613,0],clear_outer=True)),
               dict(name='candidate-fork-axial-cut',direction=[0,1,0],section=dict(normal=[0,1,0],point=[.033,1.612,0],clear_outer=True))]
atomic_json(O/'blender-review-spec.json',base)
atomic_json(O/'preparation.json',dict(source_catalogue_sha256=hashlib.sha256(ap.read_bytes()).hexdigest(),wall_file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),contexts=context,
    actual_source_corner_normals_used=True,crop_creates_no_anatomical_caps=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-t))
print(json.dumps(dict(spec=str(O/'blender-review-spec.json'),nonforward_faces=int(np.any(dots<=0,axis=1).sum()),contexts=len(context))))
