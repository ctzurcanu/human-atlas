"""Assign continuous finite PV roots without overlapping existing atrial cells."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import subtract_occupied
from geometry_tools.algorithms.source_face_complement import (
    complement_source_faces,conform_surface_boundary_edges,
)

ap=argparse.ArgumentParser()
ap.add_argument('--source',default='finite-pv-roots-v136')
ap.add_argument('--stage',required=True)
ap.add_argument('--defer-conformity',action='store_true')
ap.add_argument('--exact-source-replay',action='store_true')
args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage
out.mkdir(exist_ok=False);source=h/args.source
source_checks=json.loads((source/'construction.json').read_text())
assert all(row['passed']for row in source_checks['roots']), 'Complete root source checks must pass before partition'
old=h/'regional-bank-v102';meta=json.loads((old/'emission.json').read_text())
z=np.load(old/'arrangement.npz');o=np.array(meta['origin'])
p=(z['positions']-o).astype('<f4').astype(float)
surfaces={sid:z[key]for sid,key in meta['surfaces'].items()}
uses=meta['cells'];cells={};parents={};rows={}
for sid in surfaces:
    owners=[name for name,uu in uses.items()if any(u['surface']==sid for u in uu)]
    rows[sid]=dict(id=sid,kind='outer'if len(owners)==1 else'shared-tissue')
for name,uu in uses.items():
    ff=[];ancestry=[]
    for u in uu:
        fs=surfaces[u['surface']];side=u['side']
        ff.append(fs if side==1 else fs[:,[0,2,1]])
        ancestry.extend((u['surface'],i,side)for i in range(len(fs)))
    cells[name]=np.concatenate(ff);parents[name]=ancestry
roots=[Mesh.load(source/f'root-{i}.npz')for i in range(4)]
offsets=np.cumsum([0]+[len(m.vertices)for m in roots[:-1]])
rv=np.concatenate([m.vertices-o for m in roots])
rf=np.concatenate([m.faces+j for m,j in zip(roots,offsets)])
names=['root:left-inferior-pulmonary-vein','root:left-superior-pulmonary-vein',
       'root:right-inferior-pulmonary-vein','root:right-superior-pulmonary-vein']
owners=np.concatenate([np.full(len(m.faces),name,dtype=object)
                       for m,name in zip(roots,names)])
print('Partitioning atrial tissue against four finite roots',flush=True)
r=subtract_occupied(p,surfaces,rows,cells,parents,uses,rv,rf,
    occupied_face_owners=owners,reconcile_shared_subdivisions=True,
    diagnostic_directory=out/'canonical-diagnostics',
    exact_source_replay=args.exact_source_replay,
    partition_free_outer=args.exact_source_replay)
atomic_json(out/'partition-evidence.json',r['evidence'])
atomic_npz(out/'raw-partition.npz',positions=r['positions'])
print('Reciprocal host partition complete',flush=True)
new_surfaces=r['surfaces'];new_uses=r['cells']
new_rows={row['id']:row for row in r['surfaceRows']}
parent_faces={sid:chart['parentFaces']for sid,chart in r['charts'].items()}
contacts=r['vascularContacts'];all_contact=[];contact_parents=[]
if args.exact_source_replay:
    vv=r['positions'];free=np.concatenate([new_surfaces[u['surface']]for u in r['vascularFreeSurfaces']])
    free_parents=np.concatenate([parent_faces[u['surface']]for u in r['vascularFreeSurfaces']])
    proof=dict(method='Exact source-parent Boolean complement',
               source_replay=r['evidence']['vascularOuterReplay'])
    for u in r['vascularFreeSurfaces']:
        del new_surfaces[u['surface']];del new_rows[u['surface']];del parent_faces[u['surface']]
else:
    for row in contacts:
        sid=row['surface'];all_contact.append(new_surfaces[sid][:,[0,2,1]])
        contact_parents.extend(parent_faces[sid])
    vv,free,free_parents,proof=complement_source_faces(r['positions'],rv,rf,
        np.concatenate(all_contact),np.array(contact_parents))
atomic_json(out/'root-free-complement.json',proof)
for name in names:
    take=owners[free_parents]==name;sid='free:'+name
    new_surfaces[sid]=free[take]
    new_rows[sid]=dict(id=sid,kind='root-free',parentSurface='occupied-outer',
                       parentNormalSign=1)
    parent_faces[sid]=free_parents[take]
    new_uses[name]=[dict(surface=sid,side=1,neighbor=None)]
    for row in contacts:
        if row['vascularOwner']==name:
            new_uses[name].append(dict(surface=row['surface'],side=-1,
                                      neighbor=row['neighbor']))
conformity=dict(deferred=args.defer_conformity)
if not args.defer_conformity:
    new_surfaces,conformity=conform_surface_boundary_edges(vv,new_surfaces,new_uses,
                                                         record_face_parents=True)
    for sid,indices in conformity.pop('sourceFaceParents').items():
        parent_faces[sid]=parent_faces[sid][indices]
atomic_json(out/'edge-conformity.json',conformity)
arrays=dict(positions=vv+o,root_positions=rv+o,root_faces=rf)
manifest=dict(schema='heart-root-reciprocal-ownership/v1',origin=o.tolist(),
    cells=new_uses,surfaces=[],source_surfaces=meta['surfaces'],
    source_words=meta['corner_words'],root_names=names,
    root_source_sha256=digest_file(source/'construction.json'),root_source_stage=args.source)
for i,(sid,faces)in enumerate(new_surfaces.items()):
    arrays[f'surface_{i}']=faces
    arrays[f'parent_{i}']=parent_faces[sid]
    manifest['surfaces'].append(dict(new_rows[sid],array=f'surface_{i}',
                                     parents=f'parent_{i}'))
atomic_npz(out/'arrangement.npz',**arrays)
atomic_json(out/'arrangement.json',manifest)
atomic_json(out/'construction.json',dict(prepared=True,
    cells=len(new_uses),interfaces=len(new_surfaces),
    source_bank_sha256=meta['bank_sha256'],
    source_roots=[m.geometry_sha for m in roots],
    elapsed_seconds=time.monotonic()-start,packed=False,installed=False))
print('Prepared',len(new_uses),'reciprocal cells',round(time.monotonic()-start,2),'s',flush=True)
