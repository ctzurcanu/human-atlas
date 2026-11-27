"""One native-guided C1 SG origin-side fit, protecting its intralingual course."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
def read(p):return json.loads(Path(p).read_text())
def require(x,msg):
    if not x:raise ValueError(msg)
def main(out):
    began=time.monotonic();fp=W/'female-basal-IL-normal-charts-v1/fit.json';prior=read(fp);require(out.is_relative_to(W)and not out.exists(),'Use a new SG fit directory');out.mkdir(parents=True)
    parts=[];proof=[];low=1.557
    for i,row in enumerate(prior['parts']):
        if ':styloglossus:'not in row['partId']:parts.append(row);continue
        require(time.monotonic()-began<45,'SG endpoint stage deadline approached');require(digest_file(row['mesh'])==row['sha256'],'Bound SG source changed');old=Mesh.load(row['mesh'])
        gp=W/('female-native-course-guides-v1/course-0'+str(i-6)+'.json');guide=read(gp);target=np.array(guide['axis_start']['point_m']);origin_vertex=int(np.argmax(old.vertices[:,1]));origin=old.vertices[origin_vertex];high=float(origin[1]);delta=target-origin;det_bound=1+min(0.,float(delta[1]))*1.5/(high-low);require(det_bound>0,'Endpoint field does not prove a positive Jacobian')
        u=np.clip((old.vertices[:,1]-low)/(high-low),0,1);weight=u*u*(3-2*u);derivative=np.where((u>0)&(u<1),6*u*(1-u)/(high-low),0.);v=old.vertices+weight[:,None]*delta
        require(np.array_equal(v[old.vertices[:,1]<=low],old.vertices[old.vertices[:,1]<=low]),'Protected SG insertion changed');j=np.broadcast_to(np.eye(3),(len(v),3,3)).copy();j[:,:,1]+=derivative[:,None]*delta;require(float(np.linalg.det(j).min())>0,'Actual field Jacobian reversed')
        attrs={k:a.copy()for k,a in old.attributes.items()};norm=attrs['normals'].astype(float);norm=np.einsum('nij,nj->ni',np.linalg.inv(j).transpose(0,2,1),norm);norm/=np.linalg.norm(norm,axis=1)[:,None];attrs['normals']=np.rint(norm*32767).astype('<i2')
        # The field is identity below low; preserve those raw charts literally.
        attrs['normals'][old.vertices[:,1]<=low]=old.attributes['normals'][old.vertices[:,1]<=low]
        m=Mesh(v,old.faces,old.frame,attrs,dict(old.source,female_SG_endpoint_guide=str(gp),shared_female_cortical_chart=False));path=out/('SG-'+('left'if i==6 else'right')+'.npz');m.save(path);check=audit(m,require_closed=True,max_pairs=250000,max_tests=2000000)
        eff=m.attributes['normals'][m.faces].astype(float)*m.attributes.get('normalSigns',np.ones(len(v)))[m.faces][:,:,None];g=np.cross(m.triangles[:,1]-m.triangles[:,0],m.triangles[:,2]-m.triangles[:,0]);bad=int(((eff*g[:,None,:]).sum(2)<=0).sum());targets=[]
        for ctx in prior['target_context']:
            require(time.monotonic()-began<45,'SG context stage deadline approached');host=Mesh.load(ctx['mesh']);arrays,summary=pair_contacts(m,host,0.,max_pairs=250000,max_tests=2000000);targets.append(dict(target=ctx['partId'],geometry_sha256=host.geometry_sha,summary=summary))
        trees={r['partId']:BVH(Mesh.load(r['mesh']))for r in prior['target_context']if ':skull_cranium:'in r['partId']};patch=np.flatnonzero(old.vertices[:,1]>=1.5719);gaps=[]
        for k in patch:
            witness=min((dict(host=n,**t.closest_point(v[k]))for n,t in trees.items()),key=lambda q:q['distance_m']);gaps.append(dict(vertex=int(k),point=v[k].tolist(),host=witness['host'],host_face=int(witness['face_id']),host_point=witness['point'].tolist(),host_barycentric=witness['barycentric'].tolist(),gap_m=float(witness['distance_m'])))
        item=dict(partId=row['partId'],source=row['mesh'],source_sha256=row['sha256'],guide=str(gp),guide_sha256=digest_file(gp),origin_vertex=origin_vertex,original_origin=origin.tolist(),target_origin=target.tolist(),delta_m=delta.tolist(),identity_below_Y_m=low,full_endpoint_Y_m=high,analytic_jacobian_lower_bound=det_bound,minimum_sampled_jacobian=float(np.linalg.det(j).min()),maximum_displacement_m=float(np.linalg.norm(v-old.vertices,axis=1).max()),physical_topology_preserved=True,protected_lower_positions_literal=True,protected_lower_normal_words_literal=True,actual_self=check,effective_corner_normal_reversals=bad,all_11_contexts=targets,origin_side_fitted_witnesses=gaps,origin_gap_range_m=[min(q['gap_m']for q in gaps),max(q['gap_m']for q in gaps)],female_interface_class='adjacent fitted native origin-side reference; no shared cortical chart or histological enthesis established',qualifying_candidate=False)
        proof.append(item);parts.append(dict(row,mesh=str(path),sha256=digest_file(path),geometry_sha256=m.geometry_sha,SG_endpoint_proof=str(out/'evidence.json')));atomic_json(out/('SG-'+str(i)+'.json'),item)
    result=dict(prior,schema='female-tongue-native-guided-SG-endpoint/v1',parts=parts,SG_endpoint_proof=proof,source_bindings={**prior['source_bindings'],str(fp):digest_file(fp),str(Path(__file__).resolve()):digest_file(__file__)},qualifying_candidate=False,installed=False,elapsed_seconds=time.monotonic()-began);atomic_json(out/'fit.json',result);atomic_json(out/'evidence.json',dict(proof=proof,qualification='Actual SG-specific finite fit; unrelated GG/HG remain failed and unqualified',installed=False));print(json.dumps(dict(output=str(out),parts=[dict(partId=p['partId'],self_pass=p['actual_self']['passed'],normal_reversals=p['effective_corner_normal_reversals'],context_crossings=sum(r['summary']['crossing_face_pairs']for r in p['all_11_contexts']),origin_gap_range_m=p['origin_gap_range_m'],jacobian=p['analytic_jacobian_lower_bound'])for p in proof],elapsed_seconds=result['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
