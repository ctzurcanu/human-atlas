"""Pack the two registered AV cusps without activating their donor references."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.tissue import pack,normal_words
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'anterior-valves-v52';out=b/'av-bank-v100';positions=[];surfaces={};words={};cells={};ids={};offset=0
for donor,name in [('FJ2421','tricuspid-anterior-cusp'),('FJ2420','mitral-anterior-cusp')]:
 assert json.loads((source/(donor+'-audit.json')).read_text())['passed'];m=Mesh.load(source/(donor+'.npz'));positions.append(m.vertices);surfaces[name]=m.faces+offset;offset+=len(m.vertices);n=m.attributes['normals'].astype(float);n/=np.linalg.norm(n,axis=1)[:,None];field=n[m.faces];normal=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];dot=np.einsum('fki,fi->fk',field,normal);field+=np.maximum(.04-dot,0)[:,:,None]*normal[:,None];field/=np.linalg.norm(field,axis=2)[:,:,None];words[name]=np.rint(field*32767).astype('<i2');cells[name]=[dict(surface=name,side=1,neighbor=None)];ids[name]='TA98HEART:male:'+name
r=pack(np.concatenate(positions),surfaces,cells,ids,m.frame,out,origin=[0,1.375,.03125],corner_words=words,max_pairs=1300000,max_tests=14000000,checkpoint_directory=b/'av-checkpoints-v100');assert r['passed'];print('Registered anterior AV cusps packed',r['total_bytes'],'bytes',round(time.monotonic()-start,2),'s',flush=True)
