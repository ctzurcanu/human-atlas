"""Saved finite stages of a shared, positive piecewise-affine cord fit.

Every source tissue chart is split on the same cage planes before transport.
No triangle is allowed to bridge two different affine transformations.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.tetrahedral_cage import RegularTetrahedralCage
from geometry_tools.algorithms.selective_plane_refinement import refine_selected_plane

ap=argparse.ArgumentParser();ap.add_argument('--course',required=True);ap.add_argument('--stage',choices=['prepare','grid','diagxy','diagxz','diagyz','map'],required=True);ap.add_argument('--source');ap.add_argument('--output',required=True);ap.add_argument('--begin',type=int,default=0);ap.add_argument('--count',type=int,default=64);ap.add_argument('--coalesce-affine',action='store_true');args=ap.parse_args();start=time.monotonic()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);course_path=Path(args.course).resolve();course=np.load(course_path/'course.npz')
minimum=np.array([-.03,1.14,-.11]);steps=np.array([.015,.002,.02]);shape=[5,221,7];cage=RegularTetrahedralCage([0.,0.,0.],1.,shape)
world=cage.nodes*steps+minimum;y=np.clip(world[:,1],course['y'][0],course['y'][-1]);values={k:PchipInterpolator(course['y'],course[k])(y)for k in ['center_x','native_center_z','section_ap_scale','cord_shift_z']}
target=world.copy();target[:,0]=values['center_x']+(world[:,0]-values['center_x'])/values['section_ap_scale'];target[:,2]=values['native_center_z']+(world[:,2]-values['native_center_z'])*values['section_ap_scale']+values['cord_shift_z']
displacement=(target-world)/steps;displacement[cage.boundary]=0.;proof=cage.check(displacement,minimum_jacobian=.2);assert proof['passed'],proof
levels=minimum[1]+np.arange(shape[1])*steps[1];aq=PchipInterpolator(course['y'],course['section_ap_scale'])(np.clip(levels,course['y'][0],course['y'][-1]));varying=(np.abs(np.diff(aq))>1e-12)|(np.abs(np.diff(1/aq))>1e-12)
def complex_faces(p,f):
    local=(p-minimum)/steps;q=local[f];used=np.unique(f)
    # Only the interior has constant transverse affine coefficients in a
    # shear-only slab; the fixed outer cage boundary cannot be coalesced.
    assert np.all(local[used]>=1-2e-9)and np.all(local[used]<=np.asarray(shape)-2+2e-9),'Coalescing is restricted to literal interior cage cells'
    lo=q[:,:,1].min(1);hi=q[:,:,1].max(1)
    assert not np.any(np.floor(hi-2e-9)>np.floor(lo+2e-9)),'Conform to every axial grid plane before coalescing'
    region=np.clip(np.floor(q[:,:,1].mean(1)).astype(int),0,len(varying)-1)
    return varying[region]
if args.stage=='prepare':
    bank=base/'female-regional-pia-bank-v1';m=json.loads((bank/'tissue.json').read_text());assert m['passed'];z=np.load(bank/'arrangement.npz');p=z['positions'];names=list(m['surfaces']);f=np.concatenate([z[m['surfaces'][n]]for n in names]);labels=np.concatenate([np.full(len(z[m['surfaces'][n]]),i,int)for i,n in enumerate(names)])
    original=json.loads((base/'female-pia-v4/cells.json').read_text());meta=dict(cells=m['cells'],surfaces=names,attributes=original['attributes'],bands=original.get('bands'),source_sha256=digest_file(bank/'arrangement.npz'),course_sha256=digest_file(course_path/'course.npz'),stages=[])
else:
    source=Path(args.source).resolve();z=np.load(source/'surface-complex.npz');p=z['positions'];f=z['faces'];labels=z['labels'];meta=json.loads((source/'surface-complex.json').read_text());assert meta['course_sha256']==digest_file(course_path/'course.npz')
if args.stage in ['grid','diagxy','diagxz','diagyz']:
    planes=cage.conformance_planes(args.stage);stop=min(len(planes),args.begin+args.count);records=[]
    assert 0<=args.begin<stop
    for index in range(args.begin,stop):
        if time.monotonic()-start>40:raise TimeoutError('Spinal conformance stage reached its native 40 second deadline')
        normal,coordinate,scope=planes[index]
        local=(p-minimum)/steps;active=cage.active_faces(local,f,scope,tolerance_m=1e-9)
        if args.coalesce_affine and args.stage.startswith('diag'):active&=complex_faces(p,f)
        r=refine_selected_plane(p,f,normal=normal/steps,coordinate=coordinate,origin=minimum,active_faces=active,plane_tolerance_m=2e-12)
        labels=labels[r['parent_faces']];p,f=r['positions'],r['faces'];records.append(dict(plane=index,new_controls=r['evidence']['new_controls'],faces=len(f)))
        if(index-args.begin)%16==15:print(json.dumps(dict(stage=args.stage,plane=index,total_planes=len(planes),faces=len(f),elapsed_seconds=time.monotonic()-start)),flush=True)
    meta['stages'].append(dict(stage=args.stage,begin=args.begin,end=stop,total_planes=len(planes),records=records))
if args.stage=='map':
    complex_mask=complex_faces(p,f)if args.coalesce_affine else np.ones(len(f),bool)
    conform=cage.triangle_conformance((p-minimum)/steps,f[complex_mask],tolerance_m=2e-9);assert conform['passed'],conform
    conform.update(coalesced_affine_faces=int((~complex_mask).sum()),proof='Constant transverse scale and one affine shear within a conformed interior axial slab'if args.coalesce_affine else'Every actual tetrahedral facet')
    mapped=cage.map((p-minimum)/steps,displacement,minimum_jacobian=.2)*steps+minimum
    atomic_npz(out/'cells.npz',positions=mapped,**{name:f[labels==i]for i,name in enumerate(meta['surfaces'])})
    atomic_json(out/'cells.json',dict(cells=meta['cells'],surfaces={name:int((labels==i).sum())for i,name in enumerate(meta['surfaces'])},attributes=meta['attributes'],bands=meta['bands'],
        construction=dict(shared_positive_cage=proof,source_triangles_conform=True,source_conformance=conform,source_cell_bank_sha256=meta['source_sha256'],
            source_course_sha256=meta['course_sha256'],coordinate_units='Anisotropic unit cage; actual geometry remains metres',analytic_map_is_piecewise_affine=True,
            central_canal_is_void=True,all_tissue_interfaces_transported_together=True,source_envelope_preserved=False),
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
atomic_npz(out/'surface-complex.npz',positions=p,faces=f,labels=labels)
meta.update(cage_proof=proof,cage_minimum_m=minimum.tolist(),cage_steps_m=steps.tolist(),cage_shape=shape,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'surface-complex.json',meta)
print(json.dumps(dict(stage=args.stage,positions=len(p),faces=len(f),minimum_cage_jacobian=proof['minimum_jacobian'],elapsed_seconds=time.monotonic()-start)),flush=True)
