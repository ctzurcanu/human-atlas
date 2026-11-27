"""Source-anchored localization of actual decoded female mucosal crossings."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,msg):
    if not x:raise ValueError(msg)
def area(t):return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
def point_edges(p,edges):
    a=edges[:,0];d=edges[:,1]-a;t=np.clip(((p-a)*d).sum(1)/(d*d).sum(1),0,1);return float(np.linalg.norm(p-a-t[:,None]*d,axis=1).min())
def main(ep,out):
    began=time.monotonic();e=read(ep);require(out.is_relative_to(W)and not out.exists(),'Use a new source-bound localization stage');out.mkdir(parents=True);rows=e['actual_world_meshes'];byid={r['part_id']:r for r in rows};new=e['new_part_ids'];native=e['retained_native_part_ids']
    for row in rows:require(digest_file(row['mesh']['path'])==row['mesh']['sha256'],'Actual emitted geometry changed')
    mp=byid['TA98TONGUE:female:mucosa'];ap=byid['TA98TONGUE:female:aponeurosis'];m=Mesh.load(mp['mesh']['path']);apo=Mesh.load(ap['mesh']['path']);sourcep=W/'female-transfer-installed-source-v1/source-00.npz';source=Mesh.load(sourcep);require(len(m.faces)==824 and len(source.faces)==824,'M source face identities changed')
    # These ranges are authoritative from the packed finite-cell surfaces.
    tp=W/'shared-coverings-exposed-apo-sl-v3/packed/tissue.json';tissue=read(tp);outer_count=tissue['parts']['mucosa']['tissueDrawGroups'][0]['count']//3;interface_start=tissue['parts']['mucosa']['tissueDrawGroups'][2]['start']//3;require(outer_count==348 and interface_start==762,'Source surface grouping changed')
    outer_points={tuple(p)for p in m.triangles[:outer_count].reshape(-1,3)};outer=m.triangles[:outer_count];unique,inv=np.unique(outer.reshape(-1,3),axis=0,return_inverse=True);faces=inv.reshape(-1,3);edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1);ue,counts=np.unique(edges,axis=0,return_counts=True);rim_edges=unique[ue[counts==1]];require(len(rim_edges)==64,'Expected literal64-edge source reflection boundary')
    classes=[]
    for i,tri in enumerate(m.triangles):
        classes.append('exposed-outer'if i<outer_count else'reciprocal-M-APO-interface'if i>=interface_start else'reflection-entry-rim'if any(tuple(p)in outer_points for p in tri)else'free-inner')
    source_outer=source.triangles[:outer_count];ancestryp=W/'exposed-mucosa-v1/source-face-provenance.json';ancestry=read(ancestryp);exposedp=W/'exposed-mucosa-v1/candidate.json';exposed=read(exposedp);wholep=Path(exposed['source']);whole=Mesh.load(wholep);parent={int(r['carrierFace']):int(r['sourceFace'])for r in ancestry['faces']+ancestry['cutFaces']};require(set(parent)==set(range(outer_count)),'Missing original outer source-face ancestry')
    wd=np.all(whole.triangles[:,:,1]>=1.633,axis=1);wb=np.all(whole.triangles[:,:,2]>=.072,axis=1);wp=np.all((whole.triangles[:,:,1]>=1.633)|(whole.triangles[:,:,2]>=.072),axis=1);dorsal=np.array([wd[parent[i]]for i in range(outer_count)]);blade=np.array([wb[parent[i]]for i in range(outer_count)]);protected=np.array([wp[parent[i]]for i in range(outer_count)]);literal_faces={int(r['carrierFace'])for r in ancestry['faces']}
    require(int(protected.sum())==162 and all(i in literal_faces for i in np.flatnonzero(protected)),'Legacy literal protected source ancestry changed');all_faces=set();cross_rows=[]
    for which,carrier in [('M',m),('APO',apo)]:
        tree=BVH(carrier)
        for j,pid in enumerate(native):
            require(time.monotonic()-began<45,'Localization approached finite deadline');r=byid[pid];host=Mesh.load(r['mesh']['path']);arrays,summary=pair_contacts(carrier,host,0.,first_bvh=tree,max_pairs=250000,max_tests=2000000);p=out/(which+'-'+str(j)+'.npz');atomic_npz(p,**arrays);mask=arrays['relation']=='crossing';pairs=arrays['pairs'][mask];ids=np.unique(pairs[:,0]);points=arrays['first'][mask];records=[]
            for fi in ids:
                fi=int(fi);outerflag=which=='M'and fi<outer_count;records.append(dict(face=fi,classification=classes[fi]if which=='M'else'APO surface; unchanged source-group identity',original_whole_source_face=parent[fi]if outerflag else None,area_m2=float(area(carrier.triangles[fi:fi+1])[0]),source_literal_dorsal=bool(dorsal[fi])if outerflag else False,source_literal_free_blade=bool(blade[fi])if outerflag else False,source_literal_protected=bool(protected[fi])if outerflag else False,triangle_m=carrier.triangles[fi].tolist()))
            if which=='M':all_faces.update(int(i)for i in ids)
            distances=[point_edges(p,rim_edges)for p in points]if len(points)else[]
            cross_rows.append(dict(first=mp['part_id']if which=='M'else ap['part_id'],second=pid,first_geometry_sha256=carrier.geometry_sha,second_geometry_sha256=host.geometry_sha,summary=summary,witnesses=str(p),witnesses_sha256=digest_file(p),affected_faces=records,class_pair_counts={c:int(sum(1 for fi in pairs[:,0]if (classes[int(fi)]if which=='M'else'APO')==c))for c in (set(classes)if which=='M'else['APO'])},protected_outer_pair_count=int(sum(protected[int(fi)]for fi in pairs[:,0]if which=='M'and fi<outer_count)),representative_actual_crossing_point_distance_to_64_edge_entry_boundary_m=dict(min=min(distances),max=max(distances),median=float(np.median(distances)))if distances else None,representative_crossing_point_bounds_m=[points.min(0).tolist(),points.max(0).tolist()]if len(points)else None))
    all_faces=sorted(all_faces);protected_hit=[i for i in all_faces if i<outer_count and protected[i]];byclass={c:sorted(i for i in all_faces if classes[i]==c)for c in set(classes)};areas=area(m.triangles);report=dict(schema='female-tongue-actual-M-entry-crossing-localization/v1',packet=str(ep),packet_sha256=digest_file(ep),source_surface_manifest=str(tp),source_surface_manifest_sha256=digest_file(tp),source_mucosa=str(sourcep),source_mucosa_sha256=digest_file(sourcep),source_ancestry=dict(path=str(ancestryp),sha256=digest_file(ancestryp)),original_whole=dict(path=str(wholep),sha256=digest_file(wholep)),script_sha256=digest_file(__file__),source_face_order_retained=True,source_surface_ranges=dict(exposed_outer=[0,347],inner_or_rim=[348,761],reciprocal_M_APO=[762,823]),source_protected_dorsal_Y_m=1.633,source_protected_free_blade_Z_m=.072,source_protected_exact_formula='all three original whole-parent corners satisfy (Y>=1.633 OR Z>=0.072)',source_protected_faces=162,reflection_boundary_edges=64,pairs=cross_rows,all_10_pairs_complete=all(r['summary']['complete_within_pair_budget']for r in cross_rows),M_crossings=sum(r['summary']['crossing_face_pairs']for r in cross_rows if r['first'].endswith(':mucosa')),APO_crossings=sum(r['summary']['crossing_face_pairs']for r in cross_rows if r['first'].endswith(':aponeurosis')),affected_M_faces=all_faces,affected_M_surface_area_m2=float(areas[all_faces].sum()),affected_by_class={c:dict(faces=ids,area_m2=float(areas[ids].sum()))for c,ids in byclass.items()},protected_outer_affected_faces=protected_hit,protected_outer_affected_area_m2=float(areas[protected_hit].sum()),entry_rim_only_correction_feasible_without_changing_protected_faces=len(protected_hit)==0,whole_tissue_reconstruction_attempted=False,qualifying_candidate=False,installed=False,elapsed_seconds=time.monotonic()-began);atomic_json(out/'evidence.json',report);print(json.dumps({k:v for k,v in report.items()if k not in ['pairs','affected_M_faces','affected_by_class']}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packet',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.packet.resolve(),a.output.resolve())
