"""Extend the native dural carrier within the fixed lumbar/sacral canal.

This stage constructs only the complete carrier, before finite dura/arachnoid
layers. Its one cranial opening is retained. S2 and the caudal course remain
reference-guided until the actual bone/context checks and root review pass.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--course');args=ap.parse_args()
start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';source=base/'female-covering-source-v1/spinal_cord_dura_mater-959.npz'
z=np.load(source);original=z['vertices'];source_faces=z['faces'];v,inv=np.unique(original,axis=0,return_inverse=True);f=inv[source_faces]
cut=1.178;end_y=.987
atomic_json(out/'plan.json',dict(schema='caudal-dural-carrier-plan/v1',source_path=str(source),source_sha256=digest_file(source),
    whole_assembly='Female cord, coverings and native root interfaces',fixed=['Entire skeleton','Every unrelated organ','Native upper dural coordinates and cranial rim'],
    interfaces=['Native upper dura retained above the shared caudal rim','Dura/arachnoid share a finite interface in the next stage',
        'Caudal sac closes at the proposed S2 level; the subarachnoid cavity is not a solid plug','Root cuffs need separate nerve-surface mouth review'],
    cut_y_m=cut,proposed_s2_y_m=end_y,skeleton_witnesses=str(base/'female-caudal-canal-v1/canal-witnesses.json'),
    references=['https://nba.uth.tmc.edu/neuroscience/m/s2/chapter03.html','https://pubmed.ncbi.nlm.nih.gov/8886963/'],
    anatomical_acceptance=False,installed=False))
points=list(v);edge_cache={};faces=[]
def cross(a,b):
    key=tuple(sorted((a,b)))
    if key not in edge_cache:
        p,q=v[list(key)];t=(cut-p[1])/(q[1]-p[1]);point=p+t*(q-p);point[1]=cut
        edge_cache[key]=len(points);points.append(point)
    return edge_cache[key]
for tri in f:
    polygon=[]
    for a,b in zip(tri,np.roll(tri,-1)):
        a,b=int(a),int(b);da=v[a,1]-cut;db=v[b,1]-cut
        if da>=0:polygon.append(a)
        if (da>=0)!=(db>=0):polygon.append(cross(a,b))
    polygon=list(dict.fromkeys(polygon))
    for j in range(1,len(polygon)-1):faces.append([polygon[0],polygon[j],polygon[j+1]])
faces=np.asarray(faces,int);points=np.asarray(points);loops=boundary_loops(faces)
lower=[r for r in loops if np.max(abs(points[r,1]-cut))<1e-10]
assert len(lower)==1,'A single actual caudal carrier rim is required'
ring=lower[0];first=points[ring];center=first.mean(0);xz=first[:,[0,2]]-center[[0,2]]
scale=np.ptp(xz,axis=0)/2;shape=xz/scale
# Use medians of the actual fixed-bone ray witnesses at a coarse set of
# anatomical stations. This suppresses individual facet/cancellous outliers;
# full bone surface crossings are still tested on the emitted carrier.
witness=json.loads((base/'female-caudal-canal-v1/canal-witnesses.json').read_text())['rows']
anchors=[]
for y in [.987,1.005,1.025,1.045,1.065,1.085,1.105,1.125,1.145,1.16]:
    near=[r for r in witness if abs(r['y_m']-y)<.0041]
    if not near:continue
    mid=float(np.median([r['canal_midpoint_z_m']for r in near]));half=float(np.quantile([r['ap_clearance_m']/2 for r in near],.25))
    anchors.append([y,mid,max(.0015,min(.0065,half-.0008))])
anchors.append([cut,float(center[2]),float(scale[1])]);anchors=np.asarray(anchors)
cz=PchipInterpolator(anchors[:,0],anchors[:,1]);rz=PchipInterpolator(anchors[:,0],anchors[:,2])
if args.course:
    course_path=Path(args.course).resolve();course=np.load(course_path)
    cz=PchipInterpolator(course['y'],course['z']);rz=PchipInterpolator(course['y'],course['radius'])
    # Resample the new portion only. The native caudal rim remains literal,
    # while 24 circumferential controls keep the long sac inexpensive.
    perimeter=np.r_[0.,np.cumsum(np.linalg.norm(xz-np.roll(xz,-1,axis=0),axis=1))]
    phase=perimeter/perimeter[-1];sample=np.linspace(0,1,24,endpoint=False)
    sampled=np.c_[[np.interp(sample,phase,np.r_[xz[:,i],xz[0,i]])for i in range(2)]].T
    native_radius=float(np.linalg.norm(xz,axis=1).max());shape=sampled/native_radius
    ys=np.linspace(cut,end_y,int(np.ceil((cut-end_y)/.0006))+1)[1:]
else:ys=np.linspace(cut,end_y,100)[1:]
rows=[];previous=ring;all_points=list(points);all_faces=list(faces)
for y in ys:
    # Preserve the literal native rim shape at the join. Distally its smooth
    # affine transverse profile follows the independently measured canal.
    ap=float(rz(y));width=np.interp(y,[end_y,1.005,1.04,1.12,cut],[.0035,.006,.008,.0085,scale[0]])
    z_center=float(cz(y))
    if args.course:
        blend=float(np.clip((cut-y)/.0025,0,1));blend=blend*blend*(3-2*blend)
        width=ap=native_radius*(1-blend)+float(rz(y))*blend
        z_center=float(center[2]*(1-blend)+cz(y)*blend)
    fade=np.clip((y-end_y)/.005,0,1)
    if fade<1:
        # A smooth, finite terminal dome closes the dural sac around its CSF
        # space. It does not cap the spinal cord, filum or a nerve.
        factor=np.sqrt(max(0,1-(1-fade)**2));width*=factor;ap*=factor
    row=np.c_[center[0]+shape[:,0]*width,np.full(len(shape),y),z_center+shape[:,1]*ap]
    if y==ys[-1]:
        tip=len(all_points);all_points.append([center[0],float(y),z_center])
        for a,b in zip(previous,np.roll(previous,-1)):all_faces.append([int(b),int(a),tip])
    else:
        current=np.arange(len(all_points),len(all_points)+len(row));all_points.extend(row)
        if len(previous)==len(current):
            for a,b,c,d in zip(previous,np.roll(previous,-1),current,np.roll(current,-1)):
                all_faces.extend([[int(b),int(a),int(c)],[int(b),int(c),int(d)]])
        else:
            # Ordered unequal-rim zipper. Only the new connecting strip is
            # triangulated; no disk closes the dural cavity at this join.
            i=j=0;n=len(previous);m=len(current)
            ta=phase if args.course else np.linspace(0,1,n+1);tb=np.linspace(0,1,m+1)
            while i<n or j<m:
                aa=int(previous[i%n]);bb=int(current[j%m])
                if j==m or i<n and ta[i+1]<=tb[j+1]:
                    all_faces.append([int(previous[(i+1)%n]),aa,bb]);i+=1
                else:all_faces.append([aa,bb,int(current[(j+1)%m])]);j+=1
        previous=current
    rows.append(dict(y_m=float(y),center_z_m=z_center,half_width_m=float(width),ap_radius_m=float(ap)))
vv=np.asarray(all_points);ff=np.asarray(all_faces,int)
mesh=Mesh(vv,ff,'model:local-female-ta98:assembled-metres')
proof=audit(mesh,require_closed=False,max_pairs=250000,max_tests=2000000)
atomic_npz(out/'surface.npz',vertices=vv,faces=ff,frame=np.array(mesh.frame))
remaining=boundary_loops(ff);assert len(remaining)==1 and np.min(vv[remaining[0],1])>1.56
atomic_json(out/'surface.json',dict(source_sha256=digest_file(source),audit=proof,stations=rows,anchors=anchors.tolist(),
    original_upper_positions_literal=True,cranial_rim_preserved=True,caudal_sac_end_y_m=end_y,
    course_sha256=digest_file(args.course)if args.course else None,
    remaining_limits=['S2 termination and canal clearance require actual surface/context review','Root cuff openings and cauda courses remain uncorrected'],
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start))
assert proof['passed'],proof['failures']
print(json.dumps(dict(triangles=len(ff),components=len(proof['topology']['components']),cranial_open_edges=proof['topology']['open_edges'],
    self_crossings=proof['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start)),flush=True)
