"""Read-only exact installed cortex/pyramid and collecting free-sheet witnesses."""
import json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_npz,atomic_json,digest_file
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'

def key(t):return b''.join(sorted(p.astype('<f8').tobytes()for p in t))

def main():
    start=time.monotonic();ev=json.loads((BASE/'evidence.json').read_text());assert digest_file(ev['catalogue'])==ev['catalogue_sha256']
    bind={r['part_id']:r for r in ev['source_bindings']}
    def load(pid):
        r=bind[pid];assert digest_file(r['mesh'])==r['mesh_sha256'];return Mesh.load(r['mesh'])
    cortex=load('RENALDETAIL:male:left:cortex');ct=cortex.vertices[cortex.faces];lookup={key(t):i for i,t in enumerate(ct)};rows=[]
    for letter in 'abcdefghijk':
        pid='RENALDETAIL:male:left:pyramid:'+letter;m=load(pid);pt=m.vertices[m.faces];pairs=[(i,lookup[k])for i,t in enumerate(pt)if(k:=key(t))in lookup]
        ids=np.array(pairs,dtype=int).reshape(-1,2);pf=ids[:,0];cf=ids[:,1];p=pt[pf];c=ct[cf]
        pn=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);cn=np.cross(c[:,1]-c[:,0],c[:,2]-c[:,0]);pa=np.linalg.norm(pn,axis=1)/2
        cosine=np.einsum('ij,ij->i',pn,cn)/np.maximum(np.linalg.norm(pn,axis=1)*np.linalg.norm(cn,axis=1),1e-300)
        out=BASE/('shared-cortex-pyramid-'+letter+'.npz');atomic_npz(out,pyramid_face_ids=pf,cortex_face_ids=cf,pyramid_triangles=p,cortex_triangles=c,geometric_normal_cosine=cosine,areas_m2=pa)
        rows.append(dict(pyramid=pid,archive=str(out),sha256=digest_file(out),literal_common_faces=len(pairs),area_m2=float(pa.sum()),
                         all_faces_opposed=bool(np.all(cosine<-.999999999)),centroid_m=(np.einsum('i,ij->j',pa,p.mean(1))/pa.sum()).tolist()if len(pairs)else None,
                         whole_pyramid_base_established=False,anatomical_patch_acceptance=False))
    # The current serialized free mucosal group is an exact physical sheet. It
    # needs inner-versus-outer donor replay before declaring all of it urine-side.
    muc=load('RENALDETAIL:male:left:collecting-wall:mucosa');r=muc.source['record'];gs=[g for g in r['tissueDrawGroups']if not g.get('neighbor')]
    ids=np.concatenate([np.arange(g['start']//3,(g['start']+g['count'])//3)for g in gs]);out=BASE/'collecting-mucosa-free-sheet.npz'
    atomic_npz(out,vertices=muc.vertices,faces=muc.faces[ids],actual_source_face_ids=ids)
    report=dict(schema='renal-blood-literal-tissue-anchors/v1',catalogue_sha256=ev['catalogue_sha256'],source_evidence_sha256=digest_file(BASE/'evidence.json'),
                interfaces=rows,total_literal_shared_faces=sum(r['literal_common_faces']for r in rows),
                collecting_free_sheet=dict(archive=str(out),sha256=digest_file(out),faces=len(ids),source_part='RENALDETAIL:male:left:collecting-wall:mucosa',
                  free_tissue_surface_is_not_entire_urine_volume=True,urine_side_not_yet_independently_replayed=True),
                no_geometry_modified=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start)
    atomic_json(BASE/'tissue-anchors.json',report);print(json.dumps(report))
if __name__=='__main__':main()
