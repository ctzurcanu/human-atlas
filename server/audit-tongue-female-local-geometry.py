"""Actual local tongue self/orientation and covering checks after basal fit."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,message):
    if not x:raise ValueError(message)
def main(fit_path,output):
    began=time.monotonic();fit=read(fit_path);require(output.is_relative_to(W)and not output.exists(),'Use a new local geometry audit stage');output.mkdir(parents=True)
    loaded={r['partId']:Mesh.load(r['mesh'])for r in fit['parts']};self_checks=[]
    for i,row in enumerate(fit['parts']):
        require(time.monotonic()-began<45,'Local geometry audit approached finite deadline');m=loaded[row['partId']];require(digest_file(row['mesh'])==row['sha256'],'Candidate geometry changed');check=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
        n=m.attributes['normals'][m.faces].astype(float)*m.attributes.get('normalSigns',np.ones(len(m.vertices)))[m.faces][:,:,None];g=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);dots=(n*g[:,None,:]).sum(2);backward=int((dots<=0).sum());self_checks.append(dict(partId=row['partId'],geometry_sha256=m.geometry_sha,check=check,backward_corner_normals=backward))
        atomic_json(output/('self-'+str(i).zfill(2)+'.json'),self_checks[-1])
    covers={pid:m for pid,m in loaded.items()if pid.endswith(':mucosa')or pid.endswith(':aponeurosis')};trees={pid:BVH(m)for pid,m in loaded.items()};pairs=[]
    for i,row in enumerate(fit['parts']):
        if row['partId']in covers or row['partId'].endswith(':superior-longitudinal')or row['partId'].endswith(':transverse-vertical-shared'):continue
        for j,(cid,cover)in enumerate(covers.items()):
            require(time.monotonic()-began<45,'Local covering audit approached finite deadline');m=loaded[row['partId']];arrays,summary=pair_contacts(m,cover,0.,first_bvh=trees[row['partId']],second_bvh=trees[cid],max_pairs=250000,max_tests=2000000);p=output/('cover-'+str(i).zfill(2)+'-'+str(j)+'.npz');atomic_npz(p,**arrays);pairs.append(dict(first=row['partId'],second=cid,first_geometry_sha256=m.geometry_sha,second_geometry_sha256=cover.geometry_sha,summary=summary,witnesses=str(p),witnesses_sha256=digest_file(p)))
    passed=all(r['check']['passed']and r['backward_corner_normals']==0 for r in self_checks)and all(r['summary']['complete_within_pair_budget']and r['summary']['crossing_face_pairs']==r['summary']['coplanar_overlap_face_pairs']==0 for r in pairs)
    report=dict(schema='female-local-tongue-actual-geometry-audit/v1',fit=str(fit_path),fit_sha256=digest_file(fit_path),self_checks=self_checks,covering_pairs=pairs,passed=passed,
        self_defects=sum(r['check']['self_intersections']['defect_pairs']for r in self_checks),backward_corner_normals=sum(r['backward_corner_normals']for r in self_checks),cover_crossings=sum(r['summary']['crossing_face_pairs']for r in pairs),
        all_jobs_complete=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began);atomic_json(output/'evidence.json',report);print(json.dumps({k:v for k,v in report.items()if k not in ['self_checks','covering_pairs']}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.fit.resolve(),a.output.resolve())
