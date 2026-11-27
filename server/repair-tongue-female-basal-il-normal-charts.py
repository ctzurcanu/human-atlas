"""Actual-face normal charts for19 warped IL basal corners; no shape edits."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,message):
    if not x:raise ValueError(message)
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new IL chart stage');fp=W/'female-basal-tongue-fit-v1/fit.json';prior=read(fp);contract=read(prior['basal_contract']);high=contract['protected_Y_m'];output.mkdir(parents=True);rows=[];proof=[]
    for row in prior['parts']:
        if'inferior-longitudinal'not in row['partId']:rows.append(row);continue
        old=Mesh.load(row['mesh']);t=old.triangles;g=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);g/=np.linalg.norm(g,axis=1)[:,None];sign=old.attributes.get('normalSigns',np.ones(len(old.vertices),dtype='i1'));eff=old.attributes['normals'][old.faces].astype(float)*sign[old.faces][:,:,None];bad=np.argwhere((eff*g[:,None,:]).sum(2)<=0)
        vertices=old.vertices.tolist();faces=old.faces.copy();attrs={k:list(v)for k,v in old.attributes.items()};changed=[]
        for face,corner in bad:
            face,corner=int(face),int(corner);index=int(old.faces[face,corner]);require(old.vertices[index,1]<high,'A protected upper IL normal corner would change');new=len(vertices);vertices.append(old.vertices[index]);
            for k,value in old.attributes.items():attrs[k].append(value[index].copy())
            word=np.rint(g[face]*32767).astype('<i2')*int(sign[index]);attrs['normals'][-1]=word;faces[face,corner]=new;changed.append(dict(face=face,corner=corner,source_vertex=index,new_chart_vertex=new,point=old.vertices[index].tolist(),source_word=old.attributes['normals'][index].tolist(),new_word=word.tolist()))
        attr={k:np.asarray(v,dtype=old.attributes[k].dtype)for k,v in attrs.items()};m=Mesh(np.asarray(vertices),faces,old.frame,attr,dict(old.source,local_normal_chart_reason='Warped basal IL actual finite-face normal reversal; positions and source-face identities unchanged'));require(np.array_equal(m.triangles,old.triangles),'Normal chart repair changed geometry');check=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000);require(check['passed'],'Geometry changed during normal-only chart repair')
        effective=m.attributes['normals'][m.faces].astype(float)*m.attributes.get('normalSigns',np.ones(len(m.vertices)))[m.faces][:,:,None];require(np.all((effective*g[:,None,:]).sum(2)>0),'Repaired IL actual corners still face backward');path=output/('IL-'+('left'if row['partId'].endswith('-left')else'right')+'.npz');m.save(path);rows.append(dict(row,mesh=str(path),sha256=digest_file(path),geometry_sha256=m.geometry_sha,normal_chart_source=row['mesh'],normal_chart_corners=len(changed)))
        proof.append(dict(partId=row['partId'],source=row['mesh'],source_sha256=row['sha256'],source_geometry_sha256=old.geometry_sha,actual_geometry_sha256=m.geometry_sha,physical_triangles_literal=True,actual_normal_corners_repaired=len(changed),all_changed_corners_below_protected_plane=True,changes=changed,actual_geometry_audit=check,backward_corner_normals=0))
    result=dict(prior,schema='female-protected-basal-IL-local-normal-charts/v1',parts=rows,IL_normal_chart_proof=proof,source_bindings={**prior['source_bindings'],str(fp):digest_file(fp),str(Path(__file__).resolve()):digest_file(__file__)},qualifying_candidate=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'fit.json',result);atomic_json(output/'evidence.json',dict(schema=result['schema'],proof=proof,only_normal_attribute_charts_changed=True,qualifying_candidate=False,installed=False));print(json.dumps(dict(output=str(output),repaired=[dict(part=r['partId'],corners=r['actual_normal_corners_repaired'],geometry_literal=r['physical_triangles_literal'])for r in proof],elapsed_seconds=result['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
