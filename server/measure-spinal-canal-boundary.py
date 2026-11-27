"""Locate selectable central-canal wall references in an existing tissue bank.

The opening is not converted into a solid tube. Records point to literal
contiguous face runs already owned by the finite grey-matter cells. This
measurement does not add parts or claim additional installed TA98 coverage.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_continuity import boundary_loops

ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True);ap.add_argument('--source',required=True);ap.add_argument('--output',required=True)
ap.add_argument('--catalogue',help='Explicit staged catalogue; default is the installed female catalogue')
args=ap.parse_args();started=time.monotonic();bank=Path(args.bank).resolve();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
atlas_path=Path(args.catalogue).resolve()if args.catalogue else ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
a=json.loads(atlas_path.read_text());m=json.loads((bank/'tissue.json').read_text());meta=json.loads((source/'cells.json').read_text())
assert m['passed']and a['finiteSpinalIntegration']['bankSha256']==m['bank_sha256']==digest_file(bank/'cells.bin')
raw=(bank/'cells.bin').read_bytes();first=next(iter(m['parts'].values()))
v=np.ndarray((first['vertexCount'],3),'<f4',buffer=raw,offset=first['positions']).astype(float)+first['geometryOrigin']
active={p['sourceId']:p for p in a['parts']if not p.get('suppressed')and p.get('sourceId','').startswith('spinal_cord:')}
ys=np.array([s['y']for s in meta['stations']]);centers=np.array([s['center']for s in meta['stations']])
records=[];selected=[];cap_faces=0
for name,p in m['parts'].items():
    if not meta['attributes'][name]['native_region'].startswith('grey-'):continue
    installed=active['spinal_cord:'+name]
    assert all(installed[key]==p[key]for key in ['positions','indices','indexCount','normalSigns'])
    f=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
    for g in p['tissueDrawGroups']:
        if g.get('neighbor'):continue
        start=g['start']//3;faces=f[start:(g['start']+g['count'])//3]
        q=v[faces];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        keep=abs(normal[:,1])<.5;cap_faces+=int((~keep).sum())
        ids=np.flatnonzero(keep)
        if not len(ids):continue
        midpoint=q[ids].mean(1);center=np.c_[np.interp(midpoint[:,1],ys,centers[:,0]),np.interp(midpoint[:,1],ys,centers[:,1])]
        radii=np.linalg.norm(midpoint[:,[0,2]]-center,axis=1)
        assert radii.max()<.0002,'An external grey chart is not the central-canal wall'
        selected.append(faces[ids])
        for run in np.split(ids,np.flatnonzero(np.diff(ids)!=1)+1):
            records.append(dict(part_id=installed['id'],chunk=installed['chunk'],positions=p['positions'],normals=p['normals'],normalSigns=p['normalSigns'],
                indices=p['indices']+4*(g['start']+int(run[0])*3),indexCount=len(run)*3,vertexCount=p['vertexCount'],geometryOrigin=p['geometryOrigin'],
                source_group_start=g['start'],source_face_range=[start+int(run[0]),start+int(run[-1])+1],
                reference_kind='lumen-boundary',physical_owner_retained=True))
assert selected
f=np.concatenate(selected);used,inv=np.unique(f,return_inverse=True);xyz,mapping=np.unique(v[used],axis=0,return_inverse=True);faces=mapping[inv.reshape(-1,3)]
mesh=Mesh(xyz,faces,'model:local-female-ta98:assembled-metres');check=topology(mesh);loops=boundary_loops(faces)
assert len(check['components'])==1 and len(loops)==2 and not any(check[k]for k in ['nonmanifold_edges','inconsistent_winding_edges','invalid_vertex_links','duplicate_faces','zero_area_faces'])
assert check['open_edges']==sum(map(len,loops))
proof=dict(term='A14.1.02.019',representation='Reference to the actual inner wall of an existing void; no independent solid',
    source_catalogue_sha256=digest_file(atlas_path),bank_sha256=m['bank_sha256'],source_metadata_sha256=digest_file(source/'cells.json'),
    face_runs=records,topology=check,physiological_end_mouths=[dict(vertices=len(loop),y_range_m=[float(xyz[loop,1].min()),float(xyz[loop,1].max())])for loop in loops],
    terminal_tissue_faces_excluded=cap_faces,finite_geometry_added=False,positions_or_indices_copied=False,
    installed=False,anatomical_acceptance=False,passed=True,elapsed_seconds=time.monotonic()-started)
atomic_json(out/'boundary.json',proof)
print(json.dumps(dict(passed=True,face_runs=len(records),triangles=len(faces),mouths=len(loops),new_geometry_bytes=0,installed=False,elapsed_seconds=time.monotonic()-started)),flush=True)
