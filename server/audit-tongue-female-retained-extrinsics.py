"""Actual six transferred carriers against five unchanged native female muscles."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,msg):
    if not x:raise ValueError(msg)
def main(fp,out):
    began=time.monotonic();fit=read(fp);pp=W/'female-transfer-native-preflight-v2/preflight.json';preflight=read(pp)
    require(out.is_relative_to(W)and not out.exists(),'Use a new mixed-assembly audit directory');require(len(fit['parts'])==6,'Audit exactly the six core/IL parts');out.mkdir(parents=True)
    native=preflight['target_parts'][1:6];require(len(native)==5,'Expected combined GG, bilateral HG and SG')
    for row in native:require(digest_file(row['mesh'])==row['sha256'],'Native female muscle changed')
    hosts={r['part_id']:Mesh.load(r['mesh'])for r in native};trees={k:BVH(v)for k,v in hosts.items()};pairs=[]
    for i,row in enumerate(fit['parts']):
        require(digest_file(row['mesh'])==row['sha256'],'Transferred array changed');m=Mesh.load(row['mesh']);tree=BVH(m)
        for j,r in enumerate(native):
            require(time.monotonic()-began<45,'Mixed assembly audit approached finite deadline');h=hosts[r['part_id']];arrays,summary=pair_contacts(m,h,0.,first_bvh=tree,second_bvh=trees[r['part_id']],max_pairs=250000,max_tests=2000000);path=out/(str(i)+'--'+str(j)+'.npz');atomic_npz(path,**arrays)
            crossings=int(summary['crossing_face_pairs']);wall=row['partId'].endswith(':mucosa')or row['partId'].endswith(':aponeurosis');status=('new tissue-covering interface crossing: unresolved'if wall else'intersection between independently modeled muscle extent carriers: shared interdigitation is anatomically plausible but no exclusive tissue partition is established')if crossings else'No proper surface crossing or coplanar area found; this alone establishes no anatomical attachment'
            # Vertex-to-triangle distances are explicit finite samples, not a
            # claim of the global continuous minimum between two surfaces.
            da=np.array([trees[r['part_id']].closest_point(p)['distance_m']for p in np.unique(m.vertices[np.unique(m.faces)],axis=0)]);db=np.array([tree.closest_point(p)['distance_m']for p in np.unique(h.vertices[np.unique(h.faces)],axis=0)])
            cr=np.asarray(arrays.get('relation',[]))=='crossing';indices=np.asarray(arrays.get('pairs',np.empty((0,2),int)))[cr]if len(cr)else np.empty((0,2),int);bounds=None
            if len(indices):bounds=[np.minimum(m.triangles[indices[:,0]].min((0,1)),h.triangles[indices[:,1]].min((0,1))).tolist(),np.maximum(m.triangles[indices[:,0]].max((0,1)),h.triangles[indices[:,1]].max((0,1))).tolist()]
            pairs.append(dict(first=row['partId'],second=r['part_id'],first_geometry_sha256=m.geometry_sha,second_geometry_sha256=h.geometry_sha,second_mesh=r['mesh'],second_mesh_sha256=r['sha256'],constraint=dict(kind='actual mixed-assembly diagnostic; no blanket contact exception',distance_m=0.,max_pairs=250000,max_tests=2000000),summary=summary,witnesses=str(path),witnesses_sha256=digest_file(path),crossing_triangle_bounds_m=bounds,interpretation=status,covering_interface=wall,anatomical_contact_accepted=False,sampled_vertex_triangle_distances_m=dict(first_to_second_min=float(da.min()),first_to_second_median=float(np.median(da)),second_to_first_min=float(db.min()),second_to_first_median=float(np.median(db)),global_surface_minimum_claimed=False)))
        atomic_json(out/('progress-'+str(i)+'.json'),dict(completed_parts=i+1,pairs=pairs,installed=False))
    report=dict(schema='female-six-carriers-native-extrinsics-actual-audit/v1',fit=str(fp),fit_sha256=digest_file(fp),preflight=str(pp),preflight_sha256=digest_file(pp),script_sha256=digest_file(__file__),native_geometry_unchanged=True,physical_meshes=11,all_30_pairs_complete=all(p['summary']['complete_within_pair_budget']for p in pairs),pairs=pairs,covering_crossings=sum(p['summary']['crossing_face_pairs']for p in pairs if p['covering_interface']),intrinsic_envelope_crossings=sum(p['summary']['crossing_face_pairs']for p in pairs if not p['covering_interface']),no_blanket_overlap_waiver=True,qualifying_candidate=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began);atomic_json(out/'evidence.json',report);print(json.dumps({k:v for k,v in report.items()if k!='pairs'}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.fit.resolve(),a.output.resolve())
