"""Independent decoded-cell/original-triangle renal context gates, one bounded cell."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'
CAND=BASE/'new-A-arterial-walls-v2';LINING=BASE/'interlobar-routing-grid-v1/urine-classifier-v1'

def tri_key(t):return b''.join(sorted(p.astype('<f8').tobytes()for p in t))

def solid_winding(mesh,point):
    t=mesh.triangles-point;a,b,c=t[:,0],t[:,1],t[:,2];la=np.linalg.norm(a,axis=1);lb=np.linalg.norm(b,axis=1);lc=np.linalg.norm(c,axis=1)
    numerator=np.einsum('ij,ij->i',a,np.cross(b,c));denom=la*lb*lc+np.einsum('ij,ij->i',a,b)*lc+np.einsum('ij,ij->i',b,c)*la+np.einsum('ij,ij->i',c,a)*lb
    return float(np.sum(2*np.arctan2(numerator,denom))/(4*np.pi))

def main(cell):
    began=time.monotonic();out=CAND/'independent-context-v2'/cell;assert not out.exists();out.mkdir(parents=True)
    ev=json.loads((BASE/'evidence.json').read_text());catalogue=Path(ev['catalogue']);assert digest_file(catalogue)==ev['catalogue_sha256'];bindings={r['part_id']:r for r in ev['source_bindings']}
    tm=json.loads((CAND/'finite-cells/tissue.json').read_text());bank=Path(tm['bank']);assert digest_file(bank)==tm['bank_sha256'];reader=AtlasReader(ROOT);record=dict(tm['parts'][cell],chunk=0)
    synthetic={'chunks':[{'url':'/local-models/'+str(bank.relative_to(ROOT/'.local-models'))}]};m=reader.read_part(synthetic,record,frame=tm['frame'])
    expected=tm['checks']['cells'][cell];assert m.geometry_sha==expected['geometry_sha256']and m.attributes_sha==expected['attributes_sha256'];tree=BVH(m)
    def load(pid):
        r=bindings[pid];assert digest_file(r['mesh'])==r['mesh_sha256'];q=Mesh.load(r['mesh']);assert q.geometry_sha==r['geometry_sha256']and q.attributes_sha==r['attributes_sha256'];return q,r
    muc,mucbind=load('RENALDETAIL:male:left:collecting-wall:mucosa');z=np.load(LINING/'classifier-surface-ownership.npz');v=z['vertices'];f=z['actual_urine_boundary_faces'];ids=z['actual_source32_face_ids'];assert len(f)==7659
    lt=v[f];st=muc.triangles[ids];assert all(tri_key(a)==tri_key(b)for a,b in zip(lt,st));ln=np.cross(lt[:,1]-lt[:,0],lt[:,2]-lt[:,0]);sn=np.cross(st[:,1]-st[:,0],st[:,2]-st[:,0]);cos=np.einsum('ij,ij->i',ln,sn)/np.maximum(np.linalg.norm(ln,axis=1)*np.linalg.norm(sn,axis=1),1e-300)
    lining=Mesh(v,f,m.frame,source=dict(actual_source_face_ids=ids.tolist(),source_part='RENALDETAIL:male:left:collecting-wall:mucosa',numerical_closure_faces_excluded=True))
    source_lining=out/'actual-urine-facing-lining.npz';lining.save(source_lining)
    contexts=[('actual-urine-facing-lining',lining,dict(path=str(source_lining),sha256=digest_file(source_lining),source_part=mucbind['part_id'],source_mesh_sha256=mucbind['mesh_sha256'],ownership_archive=str(LINING/'classifier-surface-ownership.npz'),ownership_archive_sha256=digest_file(LINING/'classifier-surface-ownership.npz'),all7659_actual_source_triangles_exact=True,all_normals_opposed_to_mucosal_tissue=bool(np.all(cos<-.999999999)),numerical_classifier_closure_faces_used=0),False)]
    for name in ['pyramid:'+c for c in 'abcdefghijk']+['papilla:'+c for c in 'abcdefghijk']+['collecting-wall:'+c for c in ['adventitia','muscularis','mucosa']]:
        pid='RENALDETAIL:male:left:'+name;q,b=load(pid);contexts.append((pid,q,b,True))
    results=[];failed=[]
    for i,(pid,q,b,closed)in enumerate(contexts):
        assert time.monotonic()-began<48,'Context stage nearing declared cell deadline'
        qt=BVH(q);r,s=pair_contacts(m,q,0.,first_bvh=tree,second_bvh=qt,max_pairs=250000,max_tests=2000000);p=out/f'pair-{i:02d}.npz';atomic_npz(p,**r)
        occupancy=None
        if closed and s['complete_within_pair_budget']and not(s['crossing_face_pairs']or s['coplanar_overlap_face_pairs']or s['touching_face_pairs']):
            # For disjoint connected closed boundaries, reciprocal source-point
            # winding checks distinguish separated solids from full containment.
            first=m.vertices[m.faces[0,0]];second=q.vertices[q.faces[0,0]];w1=solid_winding(q,first);w2=solid_winding(m,second);d1=qt.closest_point(first);d2=tree.closest_point(second)
            occupancy=dict(candidate_point_inside_context=abs(w1)>.5,context_point_inside_candidate=abs(w2)>.5,candidate_point_winding=w1,context_point_winding=w2,
              candidate_point_m=first.tolist(),context_point_m=second.tolist(),candidate_point_surface_distance_m=d1['distance_m'],context_point_surface_distance_m=d2['distance_m'],
              method='Float64 analytical triangle solid-angle sums plus original-triangle BVH distances; complete disjoint connected boundaries required',
              neither_boundary_is_a_closed_lumen_classifier=True)
        bad=not s['complete_within_pair_budget']or bool(s['crossing_face_pairs']or s['coplanar_overlap_face_pairs']or s['touching_face_pairs'])or(bool(occupancy)and(occupancy['candidate_point_inside_context']or occupancy['context_point_inside_candidate']))
        if bad:failed.append(pid)
        row=dict(context=pid,context_geometry_sha256=q.geometry_sha,context_attributes_sha256=q.attributes_sha,context_source=b,complete_actual_triangle_gate=s,closed_solid_containment=occupancy,pairs=str(p),pairs_sha256=digest_file(p),passed=not bad)
        results.append(row);atomic_json(out/'progress.json',dict(cell=cell,completed=results,complete=False,failed_contexts=failed,elapsed_seconds=time.monotonic()-began),overwrite=True)
    result=dict(schema='independent-new-renal-arterial-context/v1',cell=cell,candidate_geometry_sha256=m.geometry_sha,candidate_attributes_sha256=m.attributes_sha,
      candidate_bank=str(bank),candidate_bank_sha256=digest_file(bank),candidate_tissue_manifest_sha256=digest_file(CAND/'finite-cells/tissue.json'),
      candidate_decoded_independently=True,source_catalogue=str(catalogue),source_catalogue_sha256=ev['catalogue_sha256'],contexts=results,context_count=len(results),
      failed_contexts=failed,passed=not failed,all_context_pairs_complete=all(r['complete_actual_triangle_gate']['complete_within_pair_budget']for r in results),
      exact_urine_facing_lining_replayed=True,closed_whole_urine_occupancy_remains_unknown=True,no_geometry_or_master_or_report_modified=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(catalogue)==ev['catalogue_sha256'];atomic_json(out/'evidence.json',result)
    print(json.dumps(dict(output=str(out/'evidence.json'),cell=cell,seconds=result['elapsed_seconds'],passed=result['passed'],failed_contexts=failed,summary=[(r['context'],r['complete_actual_triangle_gate']['crossing_face_pairs'],r['complete_actual_triangle_gate']['touching_face_pairs'])for r in results])))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cell',required=True,choices=['interlobar','arcuate','cortical-radiate']);main(p.parse_args().cell)
