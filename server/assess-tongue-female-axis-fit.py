"""One anatomical-axis-preserving tongue-body similarity, with no pose search."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.fit import transform_mesh
W=ROOT/'.local-models/geometry-workbench/tongue-completion';FRAME='model:local-female-ta98:assembled-metres'
def read(p):return json.loads(Path(p).read_text())
def require(x,message):
    if not x:raise ValueError(message)
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use a new axis-preserving stage');prior_path=W/'female-common-tongue-fit-v1/fit.json';prior=read(prior_path)
    for path,sha in prior['source_bindings'].items():require(digest_file(path)==sha,'Frozen source/target changed: '+path)
    sm=Mesh.load(next(r['source_mesh']for r in prior['parts']if r['partId'].endswith(':mucosa')));tm=Mesh.load(prior['native_female_tongue']);su=np.unique(sm.faces[:348]);tu=np.unique(tm.faces)
    def witness(m,used,axis,maximum):
        i=int(used[(np.argmax if maximum else np.argmin)(m.vertices[used,axis])]);return dict(vertex=i,point=m.vertices[i].tolist(),geometry_sha256=m.geometry_sha)
    sa=witness(sm,su,2,True);ta=witness(tm,tu,2,True);sp=witness(sm,su,2,False);tp=witness(tm,tu,2,False);sr=witness(sm,su,1,True);tr=witness(tm,tu,1,True)
    sl=sa['point'][2]-sp['point'][2];tl=ta['point'][2]-tp['point'][2];scale=tl/sl;require(.8<scale<1.,'Actual anterior/posterior body length does not support the expected bounded reduction')
    s_mid=(sm.vertices[su,0].max()+sm.vertices[su,0].min())/2;t_mid=(tm.vertices[tu,0].max()+tm.vertices[tu,0].min())/2
    translation=np.array([t_mid-scale*s_mid,tr['point'][1]-scale*sr['point'][1],ta['point'][2]-scale*sa['point'][2]]);matrix=np.eye(4);matrix[:3,:3]*=scale;matrix[:3,3]=translation
    output.mkdir(parents=True);parts=[]
    for i,row in enumerate(prior['parts']):
        source=Mesh.load(row['source_mesh']);m=transform_mesh(source,matrix,FRAME);path=output/('candidate-'+str(i).zfill(2)+'.npz');m.save(path);require(np.array_equal(source.faces,m.faces),'Common transform changed source topology');require(np.array_equal(source.attributes['normals'],m.attributes['normals']),'Axis-preserving common transform changed raw normal words')
        parts.append(dict(partId=row['partId'],source_mesh=row['source_mesh'],source_sha256=row['source_sha256'],source_geometry_sha256=source.geometry_sha,mesh=str(path),sha256=digest_file(path),geometry_sha256=m.geometry_sha,faces=len(m.faces)))
    anchors=dict(source_apex=sa,target_apex=ta,source_posterior=sp,target_posterior=tp,source_roof=sr,target_roof=tr,source_mid_sagittal_x=float(s_mid),target_mid_sagittal_x=float(t_mid))
    report=dict(prior,schema='installed-tongue-female-axis-body-similarity/v1',matrix=matrix.tolist(),metrics=dict(scale=float(scale),rotation_determinant=1.,rotation_degrees=0.,source_AP_length_m=sl,target_AP_length_m=tl),
        landmarks=anchors,parts=parts,registration_decision='Uniform length ratio; exact roofY/apexZ alignment; bilateral surface mid-sagittalX alignment; native anatomical axes retained with zero rotation',
        no_native_open_height_box_fit=True,no_condylar_width_fit=True,no_asymmetric_extreme_pose_fit=True,body_placement_review_pending=True,extrinsic_attachment_correction_pending=True,
        source_bindings={**prior['source_bindings'],str(prior_path):digest_file(prior_path),str(Path(__file__).resolve()):digest_file(__file__)},elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'fit.json',report)
    spec=dict(frame=FRAME,resolution=850,wall_seconds=60,anatomical_acceptance=False,assembled_tissues=True,meshes=[dict(name=r['partId'],part_id=r['partId'],mesh=r['mesh'],color=[.76,.38,.39,1.])for r in parts]+[dict(name=r['partId'],part_id=r['partId'],mesh=r['mesh'],color=[.81,.79,.68,1.])for r in prior['target_context']if any(k in r['partId']for k in [':mandible:',':hyoid:',':skull_cranium:'])],
        views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[1,0,0]),dict(name='dorsal',direction=[0,1,0]),dict(name='sagittal',direction=[1,0,0],section=dict(clear_outer=True,normal=[1,0,0],point=[float(t_mid),1.545,.05]))])
    atomic_json(output/'review-spec.json',spec);print(json.dumps(dict(output=str(output),matrix=matrix.tolist(),metrics=report['metrics'],elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
