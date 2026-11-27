"""One transform per explicit dependency closure; candidates never install here."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize,LinearConstraint
from .common import atomic_json,digest_file,identity
from .model import Mesh
from .registry import load_registered
from .constraints import closure,check
from .audit import audit


def landmark_transform(source,target,*,weights=None,similarity=False,scale_bounds=(.5,2.)):
    x,y=np.asarray(source,float),np.asarray(target,float)
    if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=3 or len(x)<3 or not np.isfinite(np.r_[x,y]).all():raise ValueError('At least three finite corresponding 3D landmarks are required')
    w=np.ones(len(x))if weights is None else np.asarray(weights,float)
    if w.shape!=(len(x),)or np.any(w<=0)or not np.isfinite(w).all():raise ValueError('Every landmark needs a positive finite weight')
    w=w/w.sum();cx=w@x;cy=w@y;xx=x-cx;yy=y-cy
    if np.linalg.matrix_rank(xx)<2 or np.linalg.matrix_rank(yy)<2:raise ValueError('Collinear landmarks do not determine a proper 3D rotation')
    u,s,vt=np.linalg.svd((xx*w[:,None]).T@yy);sign=np.ones(3);sign[-1]=np.sign(np.linalg.det(vt.T@u.T));rotation=vt.T@np.diag(sign)@u.T
    scale=float(np.dot(s,sign)/np.sum(w*np.sum(xx*xx,axis=1)))if similarity else 1.
    if not scale_bounds[0]<=scale<=scale_bounds[1]or scale<=0:raise ValueError('Requested landmark scale is outside the declared positive scale bounds')
    translation=cy-scale*rotation@cx;matrix=np.eye(4);matrix[:3,:3]=scale*rotation;matrix[:3,3]=translation
    residual=scale*x@rotation.T+translation-y
    return matrix,dict(scale=scale,rotation_determinant=float(np.linalg.det(rotation)),landmark_rms_m=float(np.sqrt(w@np.sum(residual**2,axis=1))),maximum_landmark_error_m=float(np.linalg.norm(residual,axis=1).max()))


def constrained_translation(desired,halfspaces,*,bounds=None):
    desired=np.asarray(desired,float)
    if desired.shape!=(3,)or not np.isfinite(desired).all():raise ValueError('Finite desired metre translation required')
    normals=np.asarray([x['normal']for x in halfspaces],float).reshape(-1,3);lower=np.asarray([x['minimum_dot_m']for x in halfspaces],float)
    if not np.isfinite(normals).all()or not np.isfinite(lower).all()or np.any(np.linalg.norm(normals,axis=1)==0):raise ValueError('Finite nonzero half-space normals required')
    constraints=[LinearConstraint(normals,lower,np.inf)]if len(normals)else[]
    result=minimize(lambda v:.5*np.sum((v-desired)**2),desired,jac=lambda v:v-desired,method='SLSQP',bounds=bounds,constraints=constraints,options={'ftol':1e-14,'maxiter':200})
    slack=normals@result.x-lower
    if not result.success or len(slack)and slack.min()<-1e-10:raise ValueError('No feasible bounded translation: '+str(result.message))
    matrix=np.eye(4);matrix[:3,3]=result.x
    return matrix,dict(translation_m=result.x.tolist(),halfspace_slack_m=slack.tolist(),iterations=result.nit,linearized_constraints_only=True,nonlinear_mesh_checks_required=True)


def transform_mesh(mesh,matrix,target_frame=None):
    matrix=np.asarray(matrix,float)
    if matrix.shape!=(4,4)or not np.isfinite(matrix).all()or not np.array_equal(matrix[3],[0,0,0,1]):raise ValueError('Finite affine 4x4 matrix required')
    linear=matrix[:3,:3];gram=linear.T@linear;scale=float(np.sqrt(np.trace(gram)/3))
    if np.linalg.det(linear)<=0 or not np.allclose(gram,np.eye(3)*scale**2,rtol=1e-10,atol=1e-14):raise ValueError('Use an explicit cage for nonuniform deformation; this operation allows only proper rigid or uniform similarity transforms')
    attrs=dict(mesh.attributes)
    if 'normals'in attrs:
        raw=attrs['normals'];rotation=linear/scale
        if not np.array_equal(rotation,np.eye(3)):
            q=raw.astype(float)@rotation.T;length=np.linalg.norm(q,axis=1);q=np.divide(q,length[:,None],out=np.zeros_like(q),where=length[:,None]>0)
            attrs['normals']=np.rint(np.clip(q,-1,1)*32767).astype('<i2')
    source=dict(mesh.source,derived_from_geometry=mesh.geometry_sha,transform=matrix.tolist())
    return Mesh(mesh.vertices@linear.T+matrix[:3,3],mesh.faces,target_frame or mesh.frame,attrs,source)


def cage_mesh(mesh,spec):
    from geometry_tools.algorithms.longitudinal_organ_cage import LongitudinalOrganCage
    cage=LongitudinalOrganCage(**spec);q=cage.map(mesh.vertices);jac=cage.jacobians(mesh.vertices);attrs=dict(mesh.attributes)
    if 'normals'in attrs:
        words=np.linalg.solve(jac.transpose(0,2,1),attrs['normals'].astype(float)[...,None])[...,0];length=np.linalg.norm(words,axis=1);words=np.divide(words,length[:,None],out=np.zeros_like(words),where=length[:,None]>0);attrs['normals']=np.rint(words*32767).astype('<i2')
    metrics=dict(minimum_sampled_jacobian=float(np.linalg.det(jac).min()),continuous_cage_invertible=True,meshed_self_intersection_proof=False)
    if cage.protected_basal_lift is not None:metrics['protected_basal_lift_jacobian_lower_bound']=cage.protected_basal_lift_jacobian_lower_bound
    return Mesh(q,mesh.faces,mesh.frame,attrs,dict(mesh.source,cage=spec,derived_from_geometry=mesh.geometry_sha)),metrics


def fit(db,spec,output):
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    if not spec.get('moving'):raise ValueError('Explicit moving mesh references required')
    group=closure(db,spec['moving'],spec.get('model'))
    if group['fixed']:raise ValueError('Moving closure contains fixed structures; choose another assembly, never move the skeleton implicitly')
    protected={db.mesh(p,spec.get('model'))['id']for p in spec.get('protected',[])}
    if protected.intersection(x['id']for x in group['members']):raise ValueError('Moving dependency closure contains an explicitly protected item')
    loaded={r['id']:load_registered(db,r['id'])[1]for r in group['members']if r['geometry_sha256']}
    if len(loaded)!=len(group['members']):raise ValueError('A dependency member lacks usable geometry; resolve it explicitly before moving the assembly')
    frames={m.frame for m in loaded.values()}
    if len(frames)!=1:raise ValueError('Dependency closure contains mixed frames')
    mode=spec['mode'];matrix=None
    if mode in ['rigid','similarity']:matrix,metrics=landmark_transform(spec['source_landmarks'],spec['target_landmarks'],weights=spec.get('weights'),similarity=mode=='similarity',scale_bounds=spec.get('scale_bounds',[.5,2.]))
    elif mode=='translation':matrix,metrics=constrained_translation(spec['translation_m'],spec.get('halfspaces',[]),bounds=spec.get('translation_bounds_m'))
    elif mode=='matrix':
        matrix=np.asarray(spec['matrix'],float);metrics={}
    elif mode=='icp':
        from .open3d_backend import icp_transform
        source_id=db.mesh(spec['registration_source'],spec.get('model'))['id']
        if source_id not in loaded:raise ValueError('ICP source must belong to the explicit moving dependency closure')
        target_ref=spec['registration_target']
        if Path(target_ref).suffix=='.npz'and Path(target_ref).is_file():target=Mesh.load(target_ref)
        else:target=load_registered(db,target_ref,spec.get('target_model',spec.get('model')))[1]
        if 'initial_matrix'in spec:initial=np.asarray(spec['initial_matrix'],float)
        else:initial,_=landmark_transform(spec['source_landmarks'],spec['target_landmarks'],weights=spec.get('weights'),similarity=spec.get('initial_similarity',False),scale_bounds=spec.get('scale_bounds',[.5,2.]))
        matrix,metrics=icp_transform(loaded[source_id],target,initial,target_frame=spec['target_frame'],**spec['icp'])
    elif mode=='cage':metrics={}
    else:raise ValueError('Mode must be rigid, similarity, translation, matrix, icp or cage')
    # Shared equal coordinates receive the same transformation for every part.
    candidates={};checks={};output.mkdir(parents=True)
    rows=[]
    for row in group['members']:
        old=loaded[row['id']]
        if mode=='cage':candidate,metrics=cage_mesh(old,spec['cage'])
        else:candidate=transform_mesh(old,matrix,spec.get('target_frame'))
        check_result=audit(candidate,require_closed=spec.get('require_closed',False),max_pairs=spec.get('max_pairs',250000));checks[row['id']]=check_result;candidates[row['id']]=candidate
        path=output/(row['id']+'.npz');candidate.save(path);rows.append(dict(id=row['id'],part_id=row['part_id'],path=str(path),sha256=digest_file(path),source_geometry_sha256=old.geometry_sha,geometry_sha256=candidate.geometry_sha))
    constraints=check(db,candidates=candidates)
    passed=all(x['passed']for x in checks.values())and constraints['passed']
    result=dict(schema='geometry-fit/v1',passed=passed,model=spec.get('model'),mode=mode,transform=matrix.tolist()if matrix is not None else None,metrics=metrics,dependency_closure=group,
                candidates=rows,geometry_checks=checks,constraint_checks=constraints,atlas_modified=False,anatomical_acceptance=False)
    atomic_json(output/'fit.json',result);return result


def add_parser(sub):
    p=sub.add_parser('fit');p.add_argument('spec');p.add_argument('--output',required=True);p.set_defaults(handler=lambda args,db:fit(db,json.loads(Path(args.spec).read_text()),args.output))
