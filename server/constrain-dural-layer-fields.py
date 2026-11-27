"""Solve the few failed coating patches with exact finite-prism constraints.

The carrier and every unrelated control stay literal. Only the generated
outer/inner coat controls in the failed one-ring patches can change.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import minimize
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.tissue import layer_arrangement
from geometry_tools.algorithms.tube_continuity import finite_wall_cells

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
z=np.load(source/'fields.npz');fields=z['fields'].copy();faces=z['faces'];record=json.loads((source/'field-check.json').read_text())
initial=fields.copy();bad=np.unique(np.concatenate([np.asarray(x,int)for x in record['prisms']['invalid_faces_by_interval']]))
active=np.unique(faces[bad]);patch_face_ids=np.flatnonzero(np.isin(faces,active).any(1));patch_used=np.unique(faces[patch_face_ids])
mapping=np.full(fields.shape[1],-1,int);mapping[patch_used]=np.arange(len(patch_used));local=fields[:,patch_used].copy();ff=mapping[faces[patch_face_ids]];variable=mapping[active]
middle=local[1];q=middle[ff];area=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)
normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])/area[:,None]
scale=np.maximum(np.linalg.norm(local[[0,2]][:,variable]-local[1,variable],axis=2),.00004)
targets=np.stack([local[0,variable],local[2,variable]])
x0=np.zeros((2,len(variable),3)).ravel()
height=np.array([np.linalg.norm(local[0,ff]-local[1,ff],axis=2).mean(1),np.linalg.norm(local[1,ff]-local[2,ff],axis=2).mean(1)])
normalizer=np.maximum(area[None,:,None]*height[:,:,None],1e-24)
def updated(x):
    result=local.copy();delta=x.reshape(2,len(variable),3)*scale[:,:,None]
    result[0,variable]=targets[0]+delta[0];result[2,variable]=targets[1]+delta[1]
    return result
def jacobian_constraints(p):
    result=[]
    for outer,inner in zip(p[:-1],p[1:]):
        a,b=outer[ff],inner[ff];delta=b-a;u,v=a[:,1]-a[:,0],a[:,2]-a[:,0];du,dv=delta[:,1]-delta[:,0],delta[:,2]-delta[:,0]
        c0=np.einsum('fi,fvi->fv',np.cross(u,v),delta)
        c1=np.einsum('fi,fvi->fv',np.cross(du,v)+np.cross(u,dv),delta);c2=np.einsum('fi,fvi->fv',np.cross(du,dv),delta)
        maximum=np.maximum(c0,c0+c1+c2);stationary=np.divide(-c1,2*c2,out=np.zeros_like(c0),where=c2!=0)
        valid=(c2!=0)&(stationary>0)&(stationary<1);extrema=c0+c1*stationary+c2*stationary**2
        maximum[valid]=np.maximum(maximum[valid],extrema[valid]);result.append(maximum)
    return (-np.asarray(result)/normalizer-1e-5).ravel()
def constraints(x):
    p=updated(x);jac=jacobian_constraints(p)
    # Normal-depth bounds prevent a mathematical rescue by collapsing either
    # tissue. Root sleeves allow the published 10–60 um arachnoid range.
    dura=np.einsum('fki,fi->fk',p[0,ff]-p[1,ff],normal)
    ara=np.einsum('fki,fi->fk',p[1,ff]-p[2,ff],normal)
    return np.r_[jac,(dura.ravel()-.00005)/.0001,(.00040-dura.ravel())/.0001,
        (ara.ravel()-.000012)/.00004,(.000060-ara.ravel())/.00004]
iteration=0
def callback(x):
    global iteration
    iteration+=1
    if iteration%10==0:print(json.dumps(dict(iteration=iteration,minimum_constraint=float(constraints(x).min()),elapsed_seconds=time.monotonic()-start)),flush=True)
solution=minimize(lambda x:float(x@x),x0,method='SLSQP',bounds=[(-2,2)]*len(x0),constraints=[dict(type='ineq',fun=constraints)],
    options=dict(maxiter=120,ftol=1e-10,eps=1e-6),callback=callback)
solved=updated(solution.x);fields[:,patch_used]=solved;prisms=finite_wall_cells(fields,faces)
atomic_npz(out/'fields.npz',fields=fields,faces=faces,frame=np.array('model:local-female-ta98:assembled-metres'))
proof=dict(source_sha256=digest_file(source/'fields.npz'),optimizer_success=bool(solution.success),message=str(solution.message),iterations=solution.nit,
    minimum_constraint=float(constraints(solution.x).min()),prisms=prisms,active_vertices=active.tolist(),patch_faces=patch_face_ids.tolist(),
    maximum_generated_control_change_m=float(np.linalg.norm(fields-initial,axis=2).max()),source_carrier_literal=bool(np.array_equal(fields[1],initial[1])),
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'field-check.json',proof)
assert solution.success and prisms['finite_wall_cells_pass']and proof['minimum_constraint']>=-1e-7,proof
positions,surfaces,cells,_=layer_arrangement(fields,faces,['dura','arachnoid'])
atomic_npz(out/'cells.npz',positions=positions,**surfaces)
atomic_json(out/'cells.json',dict(cells=cells,surfaces={k:len(v)for k,v in surfaces.items()},attributes={n:dict(native_region=n,axial_region='whole-covering')for n in cells},
    construction=dict(finite_wall_cells_pass=True,source_carrier_literal=True,field_constraint_proof=str(out/'field-check.json'),
        source_carrier_role='Dura/arachnoid interface; the native source did not establish an outer histological boundary',
        cranial_rim_is_annular=True,dura_arachnoid_interface_shared=True,subarachnoid_cavity_not_plugged=True),
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(passed=True,variables=len(x0),patch_faces=len(ff),iterations=solution.nit,
    maximum_control_change_um=proof['maximum_generated_control_change_m']*1e6,elapsed_seconds=time.monotonic()-start)),flush=True)
