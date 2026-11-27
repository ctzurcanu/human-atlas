"""Fit one shared cord/meningeal course in a literal fixed native spine.

Only an axial Z displacement field is proposed for the cord: Y and X stay
literal and its continuous Jacobian is one. Meningeal sections are separate
reference-guided envelopes around it; unsupported native cuffs are not copied.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import PchipInterpolator
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.open3d_backend import SurfaceIndex,library

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--area-preserving-shapes',action='store_true');ap.add_argument('--adaptive-clearance',action='store_true');args=ap.parse_args();start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';source=base/'female-surface-v1/surface.npz';zz=np.load(source);tri=zz['positions'][zz['faces']]
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(path.read_text());reader=AtlasReader();bones=[];pp=[];ff=[];offset=0
for part in a['parts']:
    pid=part['id']
    if part.get('suppressed')or not(pid.startswith('LOCAL:female:')and part['system']=='skeletal'and('vertebra'in pid or ':sacrum:'in pid)):continue
    mesh=reader.read_part(a,part,frame='model:local-female-ta98:assembled-metres',catalogue_path=path)
    bones.append(dict(part_id=pid,geometry_sha256=mesh.geometry_sha));pp.append(mesh.vertices);ff.append(mesh.faces+offset);offset+=len(mesh.vertices)
index=SurfaceIndex(Mesh(np.concatenate(pp),np.concatenate(ff),'model:local-female-ta98:assembled-metres'))
ys=np.linspace(.987,float(tri[:,:,1].max())-.00001,400);old=np.load(base/'female-canal-course-v1/course.npz');native_min=float(tri[:,:,1].min());native_max=float(tri[:,:,1].max());centers=[];dims=[]
for y in ys:
    if y<native_min:
        z=float(np.interp(y,old['y'],old['z']));r=float(np.clip(np.interp(y,old['y'],old['radius'])-.00045,.0013,.0055));centers.append([0.,z]);dims.append([r,r]);continue
    cuts=[]
    for i,j in [(0,1),(1,2),(2,0)]:
        selected=(tri[:,i,1]-y)*(tri[:,j,1]-y)<0;pa,pb=tri[selected,i],tri[selected,j]
        if len(pa):cuts.append(pa+(pb-pa)*(y-pa[:,1])[:,None]/(pb[:,1]-pa[:,1])[:,None])
    p=np.concatenate(cuts)[:,[0,2]]
    # A sloping artificial cranial truncation must not be mistaken for a
    # tapering anatomical cord or a reason to close the meningeal opening.
    if y>=1.559:p=np.concatenate([p,zz['positions'][zz['positions'][:,1]>=1.558][:,[0,2]]])
    lo,hi=p.min(0),p.max(0);center=(lo+hi)/2
    centers.append(center);dims.append((hi-lo)/2+.00115)
centers=np.asarray(centers);dims=np.asarray(dims)
# Smooth the cavity around the conus without treating the terminal cord point
# as the end of the dural sac. Subsequent literal-triangle checks are required.
for k in range(2):
    centers[:,k]=gaussian_filter1d(centers[:,k],1.2);dims[:,k]=gaussian_filter1d(dims[:,k],1.2)
native_centers=centers.copy();deltas=np.linspace(-.006,.004,51);shapes=np.array([1.,.95,1.05,.90,1.10,.85,1.15])if args.area_preserving_shapes else np.array([1.])
base_grid=np.tile(deltas,len(shapes));base_shape=np.repeat(shapes,len(deltas));theta=np.arange(32)*2*np.pi/32
margins=np.array([.00115,.00085,.00060,.00045])if args.adaptive_clearance else np.array([.00115])
clearances=[];enclosures=[];eccentricities=[]
pia=np.load(base/'female-pia-v4/pial-fields.npz');pial_tri=pia['outer'][pia['faces']];profiles=[]
for y in ys:
    if y<native_min or y>float(pial_tri[:,:,1].max()):profiles.append(None);continue
    cuts=[]
    for i,j in [(0,1),(1,2),(2,0)]:
        selected=(pial_tri[:,i,1]-y)*(pial_tri[:,j,1]-y)<0;pa,pb=pial_tri[selected,i],pial_tri[selected,j]
        if len(pa):cuts.append(pa+(pb-pa)*(y-pa[:,1])[:,None]/(pb[:,1]-pa[:,1])[:,None])
    p=np.concatenate(cuts)[:,[0,2]]if cuts else None
    if y>=1.559:p=pia['outer'][pia['outer'][:,1]>=1.558][:,[0,2]]
    profiles.append(p)
for margin_value in margins:
    ss=np.where(ys[:,None]>=native_min,base_shape[None,:],1.)
    rx0=np.where(ys[:,None]>=native_min,(dims[:,None,0]-.00115)/ss+margin_value,dims[:,None,0])
    rz0=np.where(ys[:,None]>=native_min,(dims[:,None,1]-.00115)*ss+margin_value,dims[:,None,1])
    points=np.empty((len(ys),len(base_grid),len(theta),3));points[:,:,:,0]=centers[:,None,None,0]+rx0[:,:,None]*np.cos(theta)[None,None,:]
    points[:,:,:,1]=ys[:,None,None];points[:,:,:,2]=centers[:,None,None,1]+base_grid[None,:,None]+rz0[:,:,None]*np.sin(theta)[None,None,:]
    tensor=library().core.Tensor(((points.reshape(-1,3)-index.origin)/index.scale).astype(np.float32))
    sd=index.scene.compute_signed_distance(tensor,nthreads=2,nsamples=3).numpy().reshape(points.shape[:-1]).astype(float)*index.scale
    contains=np.ones(rx0.shape,bool);eccentric=np.zeros(rx0.shape)
    for i,p in enumerate(profiles):
        if p is None:continue
        dx=(p[:,None,0]-centers[i,0])/base_shape[None,:];dz=(p[:,None,1]-centers[i,1])*base_shape[None,:]
        contains[i]=False
        for beta in [0.,-.0002,.0002,-.0004,.0004,-.0006,.0006]:
            rho=(dx/(rx0[i]-.00018)[None,:])**2+((dz-beta)/(rz0[i]-.00018)[None,:])**2
            take=(rho.max(0)<.98)&~contains[i];eccentric[i,take]=beta;contains[i,take]=True
    clearances.append(sd.min(2));enclosures.append(contains);eccentricities.append(eccentric);print(json.dumps(dict(margin_mm=float(margin_value*1000),feasible_stations=int(((clearances[-1]>=.00015)&contains).any(1).sum()),elapsed_seconds=time.monotonic()-start)),flush=True)
grid=np.tile(base_grid,len(margins));shape=np.tile(base_shape,len(margins));margin=np.repeat(margins,len(base_grid));state_shape=np.where(ys[:,None]>=native_min,shape[None,:],1.)
clearance=np.concatenate(clearances,axis=1);enclosure=np.concatenate(enclosures,axis=1);eccentric=np.concatenate(eccentricities,axis=1);valid=(clearance>=.00015)&enclosure;bad=np.flatnonzero(~valid.any(1))
atomic_npz(out/'search.npz',y=ys,centers=centers,radii=dims,grid=grid,shape=shape,margin=margin,clearance=clearance,valid=valid)
assert not len(bad),dict(no_feasible_section_y_m=ys[bad].tolist(),best_clearance_m=clearance[bad].max(1).tolist())
node_cost=((grid[None,:]-eccentric)/.001)**2*.025+np.maximum(0,(.0008-clearance)/.001)**2+((1-state_shape)/.05)**2*.08+((.00115-margin[None,:])/.0003)**2*.12+(eccentric/.0004)**2*.02
cost=np.where(valid[0],node_cost[0],1e20);parents=np.zeros((len(ys),len(grid)),int)
penalty=((grid[:,None]-grid[None,:])/.0006)**2+((shape[:,None]-shape[None,:])/.025)**2+((margin[:,None]-margin[None,:])/.0002)**2
for i in range(1,len(ys)):
    values=cost[:,None]+penalty+((eccentric[i-1,:,None]-eccentric[i,None,:])/.0004)**2;parents[i]=values.argmin(0);cost=np.where(valid[i],values.min(0)+node_cost[i],1e20)
selected=np.empty(len(ys),int);selected[-1]=int(cost.argmin())
for i in range(len(ys)-1,0,-1):selected[i-1]=parents[i,selected[i]]
shift=grid[selected];smooth=gaussian_filter1d(shift,1.3);chosen_shape=gaussian_filter1d(state_shape[np.arange(len(ys)),selected],1.3);chosen_margin=gaussian_filter1d(margin[selected],1.3);beta=gaussian_filter1d(eccentric[np.arange(len(ys)),selected],1.3)
final_rx=np.where(ys>=native_min,(dims[:,0]-.00115)/chosen_shape+chosen_margin,dims[:,0]);final_rz=np.where(ys>=native_min,(dims[:,1]-.00115)*chosen_shape+chosen_margin,dims[:,1])
# Save both the discrete feasible solution and its smooth interpolation.
# Neither sampled sign nor smoothing is an anatomical or crossing certificate.
atomic_npz(out/'course.npz',y=ys,center_x=centers[:,0],center_z=centers[:,1]+smooth,rx=final_rx,rz=final_rz,
    cord_shift_z=smooth-beta,cord_eccentricity_z=beta,section_ap_scale=chosen_shape,outer_envelope_margin=chosen_margin,discrete_cord_shift_z=shift,native_center_z=centers[:,1])
atomic_json(out/'course.json',dict(source_surface_sha256=digest_file(source),source_catalogue_sha256=digest_file(path),fixed_bones=bones,
    map='Y remains literal; X scales reciprocally to Z about one shared centre; Z also follows one smooth displacement. Exact positive unit Jacobian and section-area preservation.',shared_assembly=['All cord grey/white cells','Pia','New dural/arachnoid course','Root attachment region'],
    minimum_ap_scale=float(chosen_shape.min()),maximum_ap_scale=float(chosen_shape.max()),
    maximum_shift_m=float(abs(smooth).max()),sampled_minimum_clearance_m=float(clearance[np.arange(len(ys)),selected].min()),
    unsupported_native_root_cuffs_not_copied=True,root_ports_and_sleeves_pending=True,search_signed_distance_not_acceptance=True,
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(stations=len(ys),maximum_shift_mm=float(abs(smooth).max()*1000),minimum_discrete_clearance_mm=float(clearance[np.arange(len(ys)),selected].min()*1000),elapsed_seconds=time.monotonic()-start)),flush=True)
