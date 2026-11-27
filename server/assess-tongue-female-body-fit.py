"""One body-guided uniform tongue transfer; fixed anatomy and origins audited.

Open female tongue is a surface guide only, never a volume or height-box target.
Distant condyles/styloids do not drive the tongue-body placement. No ICP.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.fit import landmark_transform,transform_mesh
W=ROOT/'.local-models/geometry-workbench/tongue-completion';FRAME='model:local-female-ta98:assembled-metres'
def require(x,message):
    if not x:raise ValueError(message)
def read(p):return json.loads(Path(p).read_text())
def landmark(m,used,kind):
    v=m.vertices;used=np.asarray(used,int);center=(v[used,0].min()+v[used,0].max())/2;low,high=np.quantile(v[used,2],[.25,.75]);middle=used[(v[used,2]>=low)&(v[used,2]<=high)]
    if kind=='apex':domain=used;i=int(domain[np.argmax(v[domain,2])])
    elif kind=='dorsum':domain=used[np.abs(v[used,0]-center)<.006];i=int(domain[np.argmax(v[domain,1])])
    elif kind=='left-margin':domain=middle;i=int(domain[np.argmax(v[domain,0])])
    elif kind=='right-margin':domain=middle;i=int(domain[np.argmin(v[domain,0])])
    else:
        domain=used[(np.abs(v[used,0]-center)<.006)&(v[used,1]>np.quantile(v[used,1],.25))];i=int(domain[np.argmin(v[domain,2])])
    return dict(point=v[i].tolist(),vertex=i,source_geometry_sha256=m.geometry_sha,domain_definition=kind+' actual source-surface vertex; margins restricted to middle longitudinal half, dorsal/posterior sites to central6mm')
def main(output):
    began=time.monotonic();require(output.is_relative_to(W)and not output.exists(),'Use one new body registration stage');prior_path=W/'female-common-tongue-fit-v1/fit.json';prior=read(prior_path)
    for path,sha in prior['source_bindings'].items():require(digest_file(path)==sha,'Source/target revision changed: '+path)
    mu_row=next(r for r in prior['parts']if r['partId'].endswith(':mucosa'));male=Mesh.load(mu_row['source_mesh']);female=Mesh.load(prior['native_female_tongue']);sm=np.unique(male.faces[:348]);tm=np.unique(female.faces);rows=[]
    for kind in ['apex','dorsum','left-margin','right-margin','posterior-dorsum']:rows.append(dict(name=kind,source=landmark(male,sm,kind),target=landmark(female,tm,kind),weight=1.,representation='explicit homologous source-surface landmark; regional correspondence remains approximate'))
    origin=next(r for r in prior['landmarks']if r['name']=='inner-mandibular-genial-origin');rows.append(dict(name=origin['name'],source=origin['source'],target=origin['target'],weight=2.,representation='actual fixed mandibular facet identified from native GG origin course'))
    def pt(r):return r.get('point',r.get('closest_point'))
    x=np.array([pt(r['source'])for r in rows]);y=np.array([pt(r['target'])for r in rows]);matrix,metrics=landmark_transform(x,y,similarity=True,weights=[r['weight']for r in rows],scale_bounds=(.8,1.05));predicted=x@matrix[:3,:3].T+matrix[:3,3]
    for r,q in zip(rows,predicted):r.update(transformed_source=q.tolist(),residual_m=float(np.linalg.norm(q-np.asarray(pt(r['target'])))))
    output.mkdir(parents=True);parts=[]
    for i,row in enumerate(prior['parts']):
        source=Mesh.load(row['source_mesh']);m=transform_mesh(source,matrix,FRAME);path=output/('candidate-'+str(i).zfill(2)+'.npz');m.save(path);parts.append(dict(partId=row['partId'],source_mesh=row['source_mesh'],source_sha256=row['source_sha256'],source_geometry_sha256=source.geometry_sha,mesh=str(path),sha256=digest_file(path),geometry_sha256=m.geometry_sha,faces=len(m.faces)))
    report=dict(prior,schema='installed-tongue-female-body-similarity/v1',matrix=matrix.tolist(),metrics=metrics,landmarks=rows,parts=parts,
        registration_decision='One body-surface+mandibular-origin similarity; condyles, hyoid and styloids are clearance/attachment targets, not forced exact global correspondences',
        no_native_open_height_box_fit=True,source_bindings={**prior['source_bindings'],str(prior_path):digest_file(prior_path),str(Path(__file__).resolve()):digest_file(__file__)},elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'fit.json',report);print(json.dumps(dict(output=str(output),matrix=matrix.tolist(),metrics=metrics,residuals_mm=[dict(name=r['name'],mm=r['residual_m']*1000)for r in rows],elapsed_seconds=report['elapsed_seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
