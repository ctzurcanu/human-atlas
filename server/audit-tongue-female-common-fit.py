"""Finite actual female context/attachment witnesses for one common transfer."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def q(tree,p):
    r=tree.closest_point(p);return {k:(v.tolist()if hasattr(v,'tolist')else v)for k,v in r.items()}
def check_fit(path):
    fit=read(path)
    for p,sha in fit['source_bindings'].items():require(digest_file(p)==sha,'A bound source/target changed: '+p)
    for r in fit['parts']:require(digest_file(r['mesh'])==r['sha256'],'Candidate array changed')
    for r in fit['target_context']:require(digest_file(r['mesh'])==r['sha256'],'Fixed target changed')
    return fit
def pair_stage(fit_path,output,kind):
    began=time.monotonic();fit=check_fit(fit_path);require(output.is_relative_to(W)and not output.exists(),'Use a new audit stage');output.mkdir(parents=True)
    def is_bone(pid):return any(tag in pid for tag in [':mandible:',':hyoid:',':skull_cranium:'])
    rows=[r for r in fit['target_context']if is_bone(r['partId'])]if kind=='bones'else[r for r in fit['target_context']if not is_bone(r['partId'])]
    targets={r['partId']:Mesh.load(r['mesh'])for r in rows};trees={name:BVH(m)for name,m in targets.items()};pairs=[]
    for i,row in enumerate(fit['parts']):
        candidate=Mesh.load(row['mesh']);tree=BVH(candidate)
        for j,other in enumerate(rows):
            require(time.monotonic()-began<45,'Context audit approached its finite deadline')
            host=targets[other['partId']];arrays,summary=pair_contacts(candidate,host,0.,first_bvh=tree,second_bvh=trees[other['partId']],max_pairs=250000,max_tests=2000000)
            path=output/(str(i).zfill(2)+'--'+str(j).zfill(2)+'.npz');atomic_npz(path,**arrays)
            pairs.append(dict(partId=row['partId'],targetPartId=other['partId'],first_geometry_sha256=candidate.geometry_sha,second_geometry_sha256=host.geometry_sha,
                constraint=dict(kind='disjoint',parameters=dict(distance_m=0.,scope='complete proper triangle crossings/coplanar area; closed-volume containment separate',max_pairs=250000,max_tests=2000000)),summary=summary,witnesses=str(path),witnesses_sha256=digest_file(path)))
        atomic_json(output/('progress-'+str(i+1).zfill(2)+'.json'),dict(completed_parts=i+1,total_parts=len(fit['parts']),pairs=pairs,qualifying_candidate=False,installed=False))
    failed=[dict(partId=r['partId'],targetPartId=r['targetPartId'],crossings=r['summary']['crossing_face_pairs'],coplanar=r['summary']['coplanar_overlap_face_pairs'])for r in pairs if r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs']or not r['summary']['complete_within_pair_budget']]
    report=dict(schema='female-common-tongue-context-audit/v1',group=kind,fit=str(fit_path),fit_sha256=digest_file(fit_path),pairs=pairs,all_pairs_complete=all(r['summary']['complete_within_pair_budget']for r in pairs),proper_crossing_pairs=sum(r['summary']['crossing_face_pairs']for r in pairs),coplanar_overlap_pairs=sum(r['summary']['coplanar_overlap_face_pairs']for r in pairs),failed_pairs=failed,
        source_geometry_modified=False,skeleton_modified=False,no_blanket_attachment_exception=True,geometry_fit_accepted=False,qualifying_candidate=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json',report);print(json.dumps({k:v for k,v in report.items()if k!='pairs'}))
def attachments(fit_path,output):
    began=time.monotonic();fit=check_fit(fit_path);require(output.is_relative_to(W)and not output.exists(),'Use a new attachment stage');old=read(W/'source-v1/source-audit.json');old_rows={r['partId']:r for r in old['contextParts']};parts={r['partId']:r for r in fit['parts']};target={r['partId']:Mesh.load(r['mesh'])for r in fit['target_context']};target_trees={pid:BVH(m)for pid,m in target.items()if'mandible'in pid or'hyoid'in pid or'skull_cranium'in pid}
    def ctx(pid):return Mesh.load(old_rows[pid]['meshSavePath'])
    source_groups={'GG':{'mandible':BVH(ctx('ZA:Mandible:Bone'))},'HG':{pid:BVH(ctx(pid))for pid in old_rows if pid.startswith('TA98HYOID:')},'SG-left':{'left-temporal':BVH(ctx('ZA:Temporal bone.l:Bone-1'))},'SG-right':{'right-temporal':BVH(ctx('ZA:Temporal bone.r:Bone-1'))}}
    result=[]
    for pid,row in parts.items():
        group='GG'if'Genioglossus'in pid else'HG'if'Hyoglossus'in pid else'SG-'+('left'if pid.endswith(':left')else'right')if'styloglossus'in pid else None
        if group is None:continue
        source=Mesh.load(row['source_mesh']);candidate=Mesh.load(row['mesh']);unique,indices=np.unique(source.vertices,axis=0,return_index=True);treegroup=source_groups[group];queries=[]
        for vertex,p in zip(indices,unique):
            require(time.monotonic()-began<45,'Attachment audit approached its finite deadline')
            witness=min((dict(host=name,**q(tree,p))for name,tree in treegroup.items()),key=lambda r:r['distance_m']);queries.append(dict(source_vertex=int(vertex),source_point=p.tolist(),source_host=witness))
        minimum=min(r['source_host']['distance_m']for r in queries);selected=[r for r in queries if r['source_host']['distance_m']<=minimum+.00005]
        hosts={name:tree for name,tree in target_trees.items()if('mandible'in name if group=='GG'else'hyoid'in name if group=='HG'else'skull_cranium'in name)}
        for r in selected:
            point=candidate.vertices[r['source_vertex']];witness=min((dict(host=name,**q(tree,point))for name,tree in hosts.items()),key=lambda x:x['distance_m']);r.update(target_point=point.tolist(),target_host=witness)
        gaps=np.array([r['target_host']['distance_m']for r in selected]);result.append(dict(partId=pid,source_origin_definition='Actual source cortex nearest-control band: minimum source distance plus50µm; explicit source control identities, not an inferred histological footprint',source_minimum_gap_m=minimum,source_controls=len(selected),target_gap_range_m=[float(gaps.min()),float(gaps.max())],target_gap_median_m=float(np.median(gaps)),witnesses=selected,
            female_interface_class='adjacent-unjoined candidate; exact shared female bone chart is not established',anatomical_attachment_accepted=False))
    # The native female tongue is an open coarse reference, not an occupancy solid.
    native=Mesh.load(fit['native_female_tongue']);native_tree=BVH(native);mu=Mesh.load(next(r['mesh']for r in fit['parts']if r['partId'].endswith(':mucosa')));mu_tree=BVH(mu)
    a=np.array([q(native_tree,p)['distance_m']for p in np.unique(mu.vertices[mu.faces[:348]].reshape(-1,3),axis=0)]);b=np.array([q(mu_tree,p)['distance_m']for p in np.unique(native.vertices,axis=0)])
    coarse=dict(native_reference_is_open=True,solid_occupancy_inferred=False,transformed_outer_to_native_surface_mm=dict(min=float(a.min()*1000),median=float(np.median(a)*1000),maximum=float(a.max()*1000)),native_to_transformed_surface_mm=dict(min=float(b.min()*1000),median=float(np.median(b)*1000),maximum=float(b.max()*1000)))
    output.mkdir(parents=True);report=dict(schema='female-common-tongue-attachment-audit/v1',fit=str(fit_path),fit_sha256=digest_file(fit_path),attachments=result,coarse_native_reference=coarse,full_fixed_bone_crossings_report_required=True,skeleton_modified=False,qualifying_candidate=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'evidence.json',report);print(json.dumps(dict(output=str(output),attachments=[{k:v for k,v in r.items()if k!='witnesses'}for r in result],coarse_native_reference=coarse,elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['bones','soft','attachments'],required=True);a=p.parse_args();(attachments if a.stage=='attachments'else lambda fit,out:pair_stage(fit,out,a.stage))(a.fit.resolve(),a.output.resolve())
