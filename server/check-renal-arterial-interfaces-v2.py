"""Exact current arterial crossing witnesses and current urine-sheet source replay."""
import json,sys,time
from pathlib import Path
import numpy as np
import point_cloud_utils as pcu
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'
OUT=BASE/'arterial-course-v2';DEST=OUT/'interfaces'

def main():
    began=time.monotonic();assert not DEST.exists();DEST.mkdir();ev=json.loads((BASE/'evidence.json').read_text());assert digest_file(ev['catalogue'])==ev['catalogue_sha256'];bindings={r['part_id']:r for r in ev['source_bindings']}
    def load(pid):
        b=bindings[pid];assert digest_file(b['mesh'])==b['mesh_sha256'];return Mesh.load(b['mesh'])
    network=load('ZA:Intrarenal arteries of left kidney');rows=[];tree=BVH(network)
    for i,pid in enumerate(['RENALDETAIL:male:left:'+n for n in ['pyramid:'+c for c in 'abcdefghijk']+['collecting-wall:'+c for c in ['adventitia','muscularis','mucosa']]]+['ZA:Left renal artery','ZA:Anterior branch of renal artery.l','ZA:Posterior branch of renal artery.l']):
        assert time.monotonic()-began<48,'Exact contact stage nearing deadline'
        m=load(pid);r,s=pair_contacts(network,m,0.,first_bvh=tree,max_pairs=20000,max_tests=600000);p=DEST/f'pair-{i}.npz';atomic_npz(p,**r)
        assert not s['truncated'],'Keep incomplete contact evidence; do not accept it'
        rows.append(dict(first='ZA:Intrarenal arteries of left kidney',second=pid,second_geometry_sha256=m.geometry_sha,pairs=str(p),pairs_sha256=digest_file(p),summary=s))
        atomic_json(DEST/'progress.json',dict(completed_pairs=rows,seconds=time.monotonic()-began,complete=False));print(json.dumps(dict(stage='pair',second=pid,crossing=s['crossing_face_pairs'],seconds=round(time.monotonic()-began,3))),flush=True)
    # Compare the current actual free mucosal facets with the original, closed
    # urine-side classifier. Store the complete deviation rather than assuming
    # an old success report applies to a newer bank.
    origin=json.loads((OUT/'collecting-fluid-source.json').read_text());assert digest_file(origin['archive'])==origin['archive_sha256'];z=np.load(origin['archive']);v=z['positions'];f=z['innerClassifier'];sheet=np.load(BASE/'collecting-mucosa-free-sheet.npz');sv=sheet['vertices'];sf=sheet['faces'];tri=sv[sf];pts=np.r_[tri.reshape(-1,3),tri.mean(1)]
    d,ids,bc=pcu.closest_points_on_mesh(np.ascontiguousarray(pts),np.ascontiguousarray(v),np.ascontiguousarray(f,dtype=np.int32));d=np.asarray(d);center_ids=np.asarray(ids)[len(tri)*3:]
    cn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);st=v[f[center_ids]];on=np.cross(st[:,1]-st[:,0],st[:,2]-st[:,0]);cos=np.einsum('ij,ij->i',cn,on)/np.maximum(np.linalg.norm(cn,axis=1)*np.linalg.norm(on,axis=1),1e-300)
    path=DEST/'current-mucosa-versus-fluid-source.npz';atomic_npz(path,actual_current_free_source_faces=sheet['actual_source_face_ids'],query_points=pts,distance_m=d,source_classifier_face_ids=ids,source_barycentric=bc,opposed_geometric_normal_cosine=cos)
    lining=dict(witnesses=str(path),sha256=digest_file(path),maximum_point_to_source_triangle_m=float(d.max()),median_distance_m=float(np.median(d)),maximum_centroid_distance_m=float(d[len(tri)*3:].max()),
                maximum_corner_distance_m=float(d[:len(tri)*3].max()),normal_cosine_range=[float(cos.min()),float(cos.max())],all_current_faces_coplanar_with_source_to_1e_8m=bool(d.max()<1e-8),
                opposite_orientation_fraction=float(np.mean(cos<-.999)),source_classifier_complete_current_volume_not_assumed=True)
    arrays={};sample=np.load(OUT/'station-tissue-witnesses.npz');N=len(sample['centers']);queries=np.r_[sample['centers'],sample['profiles'].reshape(-1,3)];ud,ufi,ubc=pcu.signed_distance_to_mesh(np.ascontiguousarray(queries),np.ascontiguousarray(z['urinePositions']),np.ascontiguousarray(z['urineFaces'],dtype=np.int32));ud=np.asarray(ud)
    atomic_npz(DEST/'reference-fluid-station-witnesses.npz',signed_distance_to_legacy_fluid_m=ud,closest_legacy_fluid_face_ids=ufi,barycentric=ubc)
    urine=[]
    for i in range(8):
        take=sample['course_ids']==i;cd=ud[N:].reshape(N,12)[take];urine.append(dict(component=i,legacy_reference_inside_profile_controls=int(np.sum(cd<0)),legacy_reference_inside_centers=int(np.sum(ud[:N][take]<0)),minimum_distance_m=float(cd.min()),
           current_urine_volume_gate=False))
    result=dict(schema='renal-current-arterial-interface-audit/v2',source_catalogue_sha256=ev['catalogue_sha256'],contacts=rows,lining_replay=lining,legacy_fluid_samples=urine,
                no_geometry_or_labels_changed=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(DEST/'evidence.json',result);print(json.dumps(dict(output=str(DEST/'evidence.json'),seconds=result['elapsed_seconds'],lining=lining,legacy_fluid=urine)))
if __name__=='__main__':main()
