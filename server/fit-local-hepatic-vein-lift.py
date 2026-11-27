"""Fit only the left hepatic-vein end using an exactly volume-preserving shear."""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.optimize import minimize,LinearConstraint
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.open3d_backend import library
from geometry_tools.hepatic_course_field import HepaticCourseField,HepaticVeinShear
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'hepatic-vein-lift-v190';out.mkdir(exist_ok=False);start=time.monotonic();a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());parts={p['id']:p for p in a['parts']};r=AtlasReader();frame='model:local-ta98:assembled-metres'
protect=HepaticCourseField(a,base/'hepatic-contained-course-v163/field.npz',base/'separated-courses-v105/evidence.json');cava=r.read_part(a,parts['TA98HEPATIC:male:region:A12.3.09.001'],frame=frame);caval_max=float(cava.vertices[np.unique(cava.faces),0].max());record=dict(schema='left-hepatic-vein-shear/v1',startX=caval_max+.001,endX=.035,shiftYZMetres=[.0035,.0008]);field=HepaticVeinShear(record)
probes=[]
for term in ['A12.3.09.006','A12.3.09.007','A12.3.09.008']:
    m=r.read_part(a,parts['TA98HEPATIC:male:region:'+term],frame=frame);f=m.faces;edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0);probes.extend([m.vertices[np.unique(f)],m.triangles.mean(1),m.vertices[edges].mean(1)])
p=np.unique(np.concatenate(probes),axis=0);eligible=protect.weight(p)>.15;p=p[eligible];weight=field.weight(p);moving=weight>1e-7;p=p[moving];weight=weight[moving]
body=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');origin=body.vertices.mean(0);lib=library();scene=lib.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(lib.core.Tensor((body.vertices-origin).astype(np.float32)),lib.core.Tensor(body.faces.astype(np.uint32)))
history=[]
for step in range(5):
    mapped=field(p);query=lib.core.Tensor((mapped-origin).astype(np.float32));d=scene.compute_signed_distance(query,nthreads=1,nsamples=7).numpy().astype(float);hit=scene.compute_closest_points(query,nthreads=1);normal=hit['primitive_normals'].numpy().astype(float);margin=.00025*weight
    history.append(dict(maximumOutsideMetres=float(d.max(initial=0)),outsideSamples=int((d>3e-7).sum())))
    if np.all(d<=3e-7):break
    matrix=normal[:,1:3]*weight[:,None];current=np.array(record['shiftYZMetres']);bound=(matrix@current-d-margin)/.001
    fit=minimize(lambda x:.5*np.dot(x,x),current/.001,jac=lambda x:x,constraints=[LinearConstraint(matrix,-np.inf,bound)],bounds=[(0,6),(-3,3)],method='SLSQP',options={'maxiter':50,'ftol':1e-11})
    if not fit.success:raise ValueError('Fixed-cava shear constraints infeasible: '+fit.message)
    record['shiftYZMetres']=(fit.x*.001).tolist();field=HepaticVeinShear(record)
else:raise ValueError('Local vein-lift iteration budget reached')
assert np.array_equal(field(cava.vertices),cava.vertices),'Caval coordinates moved'
record.update(passedContainment=True,history=history,cavalMaximumX=caval_max,cavalSurfaceLiteral=True,determinantIdenticallyOne=True,sourceCatalogueSha256=digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'),seconds=time.monotonic()-start)
atomic_json(out/'field.json',record);print(json.dumps(record),flush=True)
