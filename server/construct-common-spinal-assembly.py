"""Transport every cord tissue together and construct its reviewed sac carrier.

The shared map preserves Y, transverse area, tissue ownership and the canal.
No skeleton or unrelated organ geometry is modified by this saved candidate.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.refine_mapped_surface import refine_mapped_surface

ap=argparse.ArgumentParser();ap.add_argument('--course',required=True);ap.add_argument('--output',required=True);ap.add_argument('--refine-field',action='store_true');ap.add_argument('--refinement-from');ap.add_argument('--encoded-source',action='store_true');args=ap.parse_args();start=time.monotonic()
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';course_path=Path(args.course).resolve();c=np.load(course_path/'course.npz')
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);core=out/'cord';core.mkdir();carrier=out/'carrier';carrier.mkdir()
source=base/'female-pia-v4';d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');v=z['positions'].copy()
source_binding=source/'cells.npz';source_surfaces={k:z[k]for k in d['surfaces']}
if args.encoded_source:
    installed=base/'female-regional-pia-bank-v1';manifest=json.loads((installed/'tissue.json').read_text());assert manifest['passed']
    arr=np.load(installed/'arrangement.npz');v=arr['positions'].copy();source_surfaces={k:arr[array]for k,array in manifest['surfaces'].items()};d['cells']=manifest['cells'];source_binding=installed/'arrangement.npz'
interpolators={key:PchipInterpolator(c['y'],c[key])for key in ['center_x','native_center_z','section_ap_scale','cord_shift_z']}
def map_points(points):
    y=np.clip(points[:,1],c['y'][0],c['y'][-1]);cx,cz,aa,dz=[interpolators[k](y)for k in interpolators]
    result=points.copy();result[:,0]=cx+(points[:,0]-cx)/aa;result[:,2]=cz+(points[:,2]-cz)*aa+dz;return result
vv=map_points(v);surfaces=source_surfaces;refinements=[]
if args.refine_field:
    names=list(surfaces);ff=np.concatenate([surfaces[k]for k in names]);labels=np.concatenate([np.full(len(surfaces[k]),i,int)for i,k in enumerate(names)])
    intervals=[]
    if args.refinement_from:
        prior=Path(args.refinement_from).resolve();provenance=json.loads((prior/'cord/cells.json').read_text());assert provenance['shared_map_limits']['y_positions_literal']
        prior_data=np.load(prior/'cord/cells.npz');prior_pia=np.concatenate([prior_data[u['surface']]if u['side']==1 else prior_data[u['surface']][:,[0,2,1]]for u in provenance['cells']['pia']]);screen=json.loads((prior/'fixed-context-screen.json').read_text());bad_faces=[]
        for check in screen['checks']:
            if not check['crossings']:continue
            bad_faces.extend(int(pair[0 if check['tissue']=='pia'else 1])for pair in check['examples'])
        assert bad_faces,'Witness-guided refinement needs actual prior triangle crossings'
        yy=prior_data['positions'][prior_pia[bad_faces],1];intervals=[[float(yy.min())-.003,float(yy.max())+.003]]
    for iteration in range(7):
        if time.monotonic()-start>45:raise TimeoutError('Shared field refinement native stage deadline')
        edges=np.unique(np.sort(np.concatenate([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]),axis=1),axis=0);error=np.zeros(len(edges))
        for t in [.25,.5,.75]:
            src=v[edges[:,0]]*(1-t)+v[edges[:,1]]*t;expected=vv[edges[:,0]]*(1-t)+vv[edges[:,1]]*t
            error=np.maximum(error,np.linalg.norm(map_points(src)-expected,axis=1))
        mask=error>(.00004 if intervals else .000012)
        if intervals:
            ev=v[edges,1];mask&=np.any([(ev.min(1)<=hi)&(ev.max(1)>=lo)for lo,hi in intervals],axis=0)
        selected=edges[mask]
        if not len(selected):break
        refined=refine_mapped_surface(v,vv,ff,selected,map_points,protect_open_rims=False,centroid_fan=False);labels=labels[refined['parent_faces']]
        v,vv,ff=refined['source_positions'],refined['positions'],refined['faces'];refinements.append(dict(iteration=iteration,edges=len(selected),faces=len(ff),maximum_edge_map_error_m=float(error.max())))
        print(json.dumps(refinements[-1]),flush=True)
        assert len(ff)<500000,'Shared field tessellation exceeds the declared face budget'
    else:raise ValueError('Shared field interpolation has not reached its explicit local error bound')
    surfaces={k:ff[labels==i]for i,k in enumerate(names)}
physical=np.unique(np.concatenate([f.ravel()for f in surfaces.values()]));aa=interpolators['section_ap_scale'](np.clip(v[physical,1],c['y'][0],c['y'][-1]))
atomic_npz(core/'cells.npz',positions=vv,**surfaces)
meta=dict(d,surfaces={k:len(f)for k,f in surfaces.items()},construction=dict(source_cell_bank_sha256=digest_file(source_binding),shared_course_sha256=digest_file(course_path/'course.npz'),
    continuous_map_positive_unit_jacobian=True,all_tissue_interfaces_transported_together=True,source_envelope_preserved=False,
    central_canal_is_void=True,regional_fields_provisional=True),
    shared_map_limits=dict(refinements=refinements,refinement_y_intervals_m=intervals if args.refine_field else[],maximum_control_displacement_m=float(np.linalg.norm(vv[physical]-v[physical],axis=1).max()),
        minimum_section_ap_scale=float(aa.min()),y_positions_literal=bool(np.array_equal(v[:,1],vv[:,1])),analytic_section_area_ratio=1.),
    anatomical_acceptance=False,installed=False)
atomic_json(core/'cells.json',meta)
# A clean complete meningeal carrier replaces the unverified native cuffs.
# Its sole cranial continuation remains open. Root ports and sleeves are an
# explicit subsequent interface stage, never silently asserted by this sac.
verts=[];faces=[];previous=None;theta=np.arange(32)*2*np.pi/32
for i,yy in enumerate(c['y']):
    xx,zz=float(c['center_x'][i]),float(c['center_z'][i]);rx,rz=float(c['rx'][i]),float(c['rz'][i])
    fade=np.clip((yy-c['y'][0])/.0045,0,1);dome=np.sqrt(max(0.,1-(1-fade)**2))
    if i==0:verts.append([xx,float(yy),zz]);previous=np.array([0]);continue
    row=np.c_[xx+rx*dome*np.cos(theta),np.full(32,yy),zz+rz*dome*np.sin(theta)]
    current=np.arange(len(verts),len(verts)+32);verts.extend(row)
    if len(previous)==1:
        for a,b in zip(current,np.roll(current,-1)):faces.append([int(previous[0]),int(a),int(b)])
    else:
        for a,b,cc,dd in zip(previous,np.roll(previous,-1),current,np.roll(current,-1)):
            faces.extend([[int(a),int(cc),int(dd)],[int(a),int(dd),int(b)]])
    previous=current
pv=np.asarray(verts);pf=np.asarray(faces,int);mesh=Mesh(pv,pf,'model:local-female-ta98:assembled-metres')
proof=audit(mesh,require_closed=False,max_pairs=250000,max_tests=2000000)
atomic_npz(carrier/'surface.npz',vertices=pv,faces=pf,frame=np.array(mesh.frame))
atomic_json(carrier/'surface.json',dict(audit=proof,source_course_sha256=digest_file(course_path/'course.npz'),carrier_role='outer dura',
    complete_central_sac=True,cranial_continuation_open=True,caudal_dome_y_m=float(c['y'][0]),
    root_ports_and_sleeves_pending=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
assert proof['passed'],proof['failures']
print(json.dumps(dict(passed=True,cord_cells=len(d['cells']),carrier_triangles=len(pf),
    maximum_tissue_displacement_mm=meta['shared_map_limits']['maximum_control_displacement_m']*1000,elapsed_seconds=time.monotonic()-start)),flush=True)
