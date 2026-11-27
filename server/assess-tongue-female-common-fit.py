"""One source-bound common similarity candidate for all installed tongue parts.

Actual bone-surface landmarks initialize a proper uniform transform. This is
not ICP or an anatomical acceptance. Fixed female geometry is never changed;
bounded audit batches retain every actual crossing witness before any install.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.fit import landmark_transform,transform_mesh
W=ROOT/'.local-models/geometry-workbench/tongue-completion';FEMALE_FRAME='model:local-female-ta98:assembled-metres'
MALE_SHA='1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39';FEMALE_SHA='770f47f802c79cf6f3e4336bf2fa7c42cf710998069221f533a1ea3f8322e470'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def vertex_landmark(mesh,ids,axis,maximum=True):
    ids=np.asarray(ids,int);require(len(ids)>0,'An anatomical landmark domain is empty');i=int(ids[(np.argmax if maximum else np.argmin)(mesh.vertices[ids,axis])]);return dict(point=mesh.vertices[i].tolist(),vertex=i,geometry_sha256=mesh.geometry_sha)
def closest(tree,point):
    q=tree.closest_point(np.asarray(point));return {k:(v.tolist()if hasattr(v,'tolist')else v)for k,v in q.items()}
def bone_seed_landmark(trees,seed):
    candidates=[dict(host=name,**closest(tree,seed))for name,tree in trees.items()];return min(candidates,key=lambda r:r['distance_m'])
def register(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use one new common female registration stage')
    source_path=W/'female-transfer-installed-source-v1/source.json';source=read(source_path);preflight_path=W/'female-transfer-native-preflight-v2/preflight.json';pre=read(preflight_path);male=Path(source['source_catalogue']);female=Path(pre['source_catalogues']['female']['path'])
    require(digest_file(male)==source['source_catalogue_sha256']==MALE_SHA and digest_file(female)==FEMALE_SHA,'Installed source/target revision changed')
    audit_path=W/'source-v1/source-audit.json';old=read(audit_path);oldrows={r['partId']:r for r in old['contextParts']};targetrows={r['part_id']:r for r in pre['target_parts']};moving={r['partId']:Mesh.load(r['mesh'])for r in source['parts']}
    for r in source['parts']:require(digest_file(r['mesh'])==r['sha256'],'Installed source export changed')
    for r in pre['target_parts']:require(digest_file(r['mesh'])==r['sha256'],'Fixed target export changed')
    def mcontext(pid):return Mesh.load(oldrows[pid]['meshSavePath'])
    sjaw=mcontext('ZA:Mandible:Bone');tjaw=Mesh.load(targetrows['LOCAL:female:mandible:960']['mesh']);shyo=mcontext('TA98HYOID:A02.1.16.002:whole');thyo=Mesh.load(targetrows['LOCAL:female:hyoid:1002']['mesh'])
    st={side:mcontext('ZA:Temporal bone.'+('l'if side=='left'else'r')+':Bone-1')for side in ['left','right']};tt={r['part_id']:Mesh.load(r['mesh'])for r in pre['target_parts']if r['role'].startswith('immutable fixed cranial')}
    sjtree=BVH(sjaw);tjtree=BVH(tjaw);sttrees={side:BVH(m)for side,m in st.items()};tttrees={name:BVH(m)for name,m in tt.items()}
    landmarks=[]
    def add(name,s,t,method,weight=1.):landmarks.append(dict(name=name,source=s,target=t,method=method,weight=weight))
    for side,sign in [('left',1),('right',-1)]:
        sa=np.flatnonzero(sjaw.vertices[:,0]*sign>.01);ta=np.flatnonzero((tjaw.vertices[:,0]-tjaw.bounds.mean(0)[0])*sign>.01)
        add('mandibular-condylar-summit-'+side,vertex_landmark(sjaw,sa,1),vertex_landmark(tjaw,ta,1),'Actual superior vertex of the corresponding fixed mandibular ramus; surface vertex, not box center')
    # The verified/native extrinsic courses identify the inner genial region;
    # every final landmark is an actual closest fixed-bone surface point.
    source_gg=np.r_[moving['ZA:Genioglossus muscle.l'].vertices,moving['ZA:Genioglossus muscle.r'].vertices];target_gg=Mesh.load(targetrows['LOCAL:female:genioglossus:239']['mesh']).vertices
    def genial(vertices,tree):
        region=np.flatnonzero((np.abs(vertices[:,0])<.006)&(vertices[:,2]>=np.quantile(vertices[:,2],.90)))
        queries=[dict(seed=vertices[i].tolist(),seed_vertex=int(i),**closest(tree,vertices[i]))for i in region]
        return min(queries,key=lambda r:r['distance_m'])
    add('inner-mandibular-genial-origin',genial(source_gg,sjtree),genial(target_gg,tjtree),'Closest actual mandibular facet to the anterior midline GG source course; no free-space landmark',2.)
    smid=np.flatnonzero(np.abs(shyo.vertices[:,0])<.006);tmid=np.flatnonzero(np.abs(thyo.vertices[:,0]-thyo.bounds.mean(0)[0])<.006)
    add('hyoid-body-anterior-midline',vertex_landmark(shyo,smid,2),vertex_landmark(thyo,tmid,2),'Anterior actual body-cortex vertex in the central6mm domain')
    for side,sign in [('left',1),('right',-1)]:
        horn=mcontext('TA98HYOID:A02.1.16.004:'+('l'if side=='left'else'r'));ta=np.flatnonzero((thyo.vertices[:,0]-thyo.bounds.mean(0)[0])*sign>.006)
        add('greater-horn-posterior-tip-'+side,vertex_landmark(horn,np.arange(len(horn.vertices)),2,False),vertex_landmark(thyo,ta,2,False),'Posterior actual greater-horn cortex vertex on corresponding lateral half')
        sm=moving['TA98TONGUE:male:styloglossus:'+side];tm=Mesh.load(targetrows['LOCAL:female:'+('l'if side=='left'else'r')+'_styloglossus:'+('243'if side=='left'else'256')]['mesh']);si=int(np.argmax(sm.vertices[:,1]));ti=int(np.argmax(tm.vertices[:,1]));s=bone_seed_landmark({side:sttrees[side]},sm.vertices[si]);t=bone_seed_landmark(tttrees,tm.vertices[ti]);s.update(seed=sm.vertices[si].tolist(),seed_vertex=si);t.update(seed=tm.vertices[ti].tolist(),seed_vertex=ti)
        add('styloid-origin-surface-'+side,s,t,'Actual closest cranial facet to the superior native/reviewed styloglossus endpoint; cranial surface witness retained')
    # BVH uses closest_point as the actual surface witness coordinate.
    def pt(row):return row.get('point',row.get('closest_point'))
    for row in landmarks:require(pt(row['source'])is not None and pt(row['target'])is not None,'Missing actual surface landmark coordinate')
    sx=np.array([pt(r['source'])for r in landmarks]);tx=np.array([pt(r['target'])for r in landmarks]);matrix,metrics=landmark_transform(sx,tx,weights=[r['weight']for r in landmarks],similarity=True,scale_bounds=(.8,1.05))
    predicted=sx@matrix[:3,:3].T+matrix[:3,3];residual=np.linalg.norm(predicted-tx,axis=1)
    for row,e,q in zip(landmarks,residual,predicted):row.update(residual_m=float(e),transformed_source=q.tolist())
    output.mkdir(parents=True);parts=[]
    for row in source['parts']:
        m=moving[row['partId']];candidate=transform_mesh(m,matrix,FEMALE_FRAME);path=output/('candidate-'+str(len(parts)).zfill(2)+'.npz');candidate.save(path);parts.append(dict(partId=row['partId'],source_mesh=row['mesh'],source_sha256=row['sha256'],source_geometry_sha256=m.geometry_sha,mesh=str(path),sha256=digest_file(path),geometry_sha256=candidate.geometry_sha,faces=len(candidate.faces)))
    bindings={str(p):digest_file(p)for p in [source_path,preflight_path,audit_path,Path(__file__).resolve(),ROOT/'geometry_tools/fit.py',male,female]}
    contexts=[dict(partId=r['part_id'],mesh=r['mesh'],sha256=r['sha256'],geometry_sha256=r['geometry_sha256'],role=r['role'])for r in pre['target_parts']if r['fixed']]
    result=dict(schema='installed-tongue-female-common-similarity/v1',source_catalogue=str(male),source_catalogue_sha256=MALE_SHA,target_catalogue=str(female),target_catalogue_sha256=FEMALE_SHA,matrix=matrix.tolist(),metrics=metrics,landmarks=landmarks,
        parts=parts,target_context=contexts,native_female_tongue=targetrows['LOCAL:female:tongue:1089']['mesh'],native_female_tongue_is_not_occupancy=True,
        dependency_closure_meshes=len(parts),one_proper_uniform_transform=True,unconstrained_ICP=False,skeleton_modified=False,raw_shared_geometry_transformed_identically=True,
        source_bindings=bindings,actual_bone_attachment_checks_pending=True,geometry_fit_accepted=False,qualifying_candidate=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'fit.json',result);print(json.dumps(dict(output=str(output),metrics=metrics,matrix=matrix.tolist(),landmark_residuals_mm=[dict(name=r['name'],mm=r['residual_m']*1000)for r in landmarks],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();register(a.output.resolve())
