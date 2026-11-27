"""Measured uniform heart similarity from native valve/apical landmarks."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'female-fit-v68';out.mkdir(exist_ok=False);reader=AtlasReader();catalogues={sex:ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json'for sex in ['male','female']};atlases={sex:json.loads(path.read_text())for sex,path in catalogues.items()}
def mesh(sex,ids):
 rows=[reader.read_part(atlases[sex],next(p for p in atlases[sex]['parts']if p['id']==pid),frame='model:'+sex,catalogue_path=catalogues[sex])for pid in ids];return rows
def centroid(rows):
 q=np.concatenate([m.triangles for m in rows]);area=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1);return np.average(q.mean(1),axis=0,weights=area)
male_ao=centroid(mesh('male',['ZA:Right coronary leaflet','ZA:Left coronary leaflet','ZA:Non-coronary leaflet']));male_pv=centroid(mesh('male',['ZA:Right semilunar leaflet of pulmonary valve','ZA:Left semilunar leaflet of pulmonary valve','ZA:Anterior semilunar leaflet of pulmonary valve']));female_ao=centroid(mesh('female',['LOCAL:female:aortic_semilunar_valve:1504']));female_pv=centroid(mesh('female',['LOCAL:female:pulmonary_semilunar_valve:1505']));male=Mesh.load(b/'whole-wall-v25/whole-wall.npz');female=mesh('female',['LOCAL:female:ventricles:1506'])[0]
def apex(m):
 q=m.vertices;return q[q[:,1]<q[:,1].min()+.001].mean(0)
x=np.array([male_ao,male_pv,apex(male)]);y=np.array([female_ao,female_pv,apex(female)]);weights=np.array([2.,2.,1.]);weights/=weights.sum();cx=np.average(x,axis=0,weights=weights);cy=np.average(y,axis=0,weights=weights);xc=x-cx;yc=y-cy;u,s,vt=np.linalg.svd(xc.T@(weights[:,None]*yc));sign=np.linalg.det(vt.T@u.T);r=vt.T@np.diag([1,1,sign])@u.T;scale=float(np.sum(weights*np.einsum('ij,ij->i',xc@r.T,yc))/np.sum(weights*np.einsum('ij,ij->i',xc,xc)));shift=cy-scale*r@cx;residual=np.linalg.norm(scale*x@r.T+shift-y,axis=1);rotation_angle=float(np.degrees(np.arccos(np.clip((np.trace(r)-1)/2,-1,1))));assert .7<scale<1.1 and rotation_angle<35
report=dict(method='Proper rigid rotation plus one uniform scale; weighted area-centroid aortic/pulmonary valve and native apical landmarks',scale=scale,rotation=r.tolist(),translation_m=shift.tolist(),rotation_angle_degrees=rotation_angle,source_landmarks_m=x.tolist(),target_landmarks_m=y.tolist(),landmark_names=['aortic-valve-area-centroid','pulmonary-valve-area-centroid','native-ventricular-apex'],residuals_m=residual.tolist(),source_catalogue_sha256={sex:digest_file(path)for sex,path in catalogues.items()},approximate_native_landmarks=True,clinical_size_established=False,skeleton_moved=False,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'fit.json',report);print(json.dumps(report),flush=True)
