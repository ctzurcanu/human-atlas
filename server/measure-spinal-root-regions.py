"""Measure native root entry patches before assigning cord segment regions.

The source nerves and skeleton are immutable. Distances are replayed against
literal cord triangles; proximity alone does not certify a root junction.
"""
from pathlib import Path
import argparse,json,re,time
import numpy as np
import nibabel as nib
from scipy.interpolate import PchipInterpolator
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.open3d_backend import SurfaceIndex

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
atlas=json.loads(catalogue.read_text());reader=AtlasReader()
surface=np.load(base/'female-surface-v1/surface.npz');cord=Mesh(surface['positions'],surface['faces'],'model:local-female-ta98:assembled-metres')
index=SurfaceIndex(cord);rows=[];centres={}
for p in atlas['parts']:
    if p.get('suppressed') or p['id'].startswith('HRA:'):continue
    c=re.search(r'cervical_spinal_nerve_c(\d+)',p['id'])
    t=re.search(r'intercostal_nerve(\d+)',p['id'])
    if not(c or t):continue
    level=('C' if c else 'T')+str(int((c or t)[1]))
    mesh=reader.read_part(atlas,p,frame=cord.frame,catalogue_path=catalogue)
    points=np.unique(mesh.vertices,axis=0);points=points[np.abs(points[:,0])<.013]
    q,evidence=index.closest_points(points);near=q['distance_m']<=.00025
    # Keep the whole source attachment patch, rather than choosing the one
    # closest vertex and silently treating it as the centre of a segment.
    patch=q['points'][near];dist=q['distance_m'][near]
    if len(patch):
        y=np.quantile(patch[:,1],[0,.1,.5,.9,1]);centres.setdefault(level,[]).append(float(y[2]))
    else:y=None
    rows.append(dict(part_id=p['id'],level=level,source_geometry_sha256=mesh.geometry_sha,
        patch_points=len(patch),patch_y_quantiles_m=None if y is None else y.tolist(),
        maximum_patch_distance_m=float(dist.max())if len(dist)else None,
        closest_distance_m=float(q['distance_m'].min()),
        certified_interface=False,native_triangle_witnesses=True))
ref=base/'pam50-reference-v1/template/PAM50_spinal_levels.nii.gz'
im=nib.load(ref);data=np.asanyarray(im.dataobj);bounds={}
for label in range(1,31):
    zz=np.flatnonzero(np.any(data==label,axis=(0,1)));bounds[label]=[float(zz[0]-.5),float(zz[-1]+.5)]
# Root-defined anchors C2..T12. The template supplies only relative spacing
# below the last named native root, not a claim of donor-specific segments.
anchors=[]
for label in range(2,21):
    level=('C'+str(label))if label<=8 else('T'+str(label-8))
    anchors.append([float(np.mean(bounds[label])),float(np.median(centres[level]))])
co1_lower=bounds[30][0]-(bounds[30][1]-bounds[30][0])
anchors.extend([[bounds[1][1],float(cord.bounds[1,1])],
                [co1_lower,float(cord.bounds[0,1])]])
anchors.sort();x,y=np.asarray(anchors).T
assert np.all(np.diff(y)>0),'Native source roots are not monotonically ordered'
# Do not flatten the last segments with a shape-preserving cubic endpoint.
# The relative-length estimate below T12 is affine; a cubic endpoint made S5
# implausibly thin even though every anchor was monotonic.
def mapping(value):return np.interp(value,x,y)
segments=[]
for label in range(1,31):
    low,high=map(float,mapping(bounds[label]));level=('C'+str(label))if label<=8 else('T'+str(label-8))if label<=20 else('L'+str(label-20))if label<=25 else('S'+str(label-25))
    segments.append(dict(name=level,y_range_m=[low,high],native_root_patch=label<=20,
        boundary_status='Root-registered approximate segment boundary; rootlets overlap and do not define an exact transverse plane'))
segments.append(dict(name='Co1',y_range_m=[float(cord.bounds[0,1]),float(mapping(bounds[30][0]))],
    native_root_patch=False,boundary_status='Distal estimated segment, assigned the same axial extent as the S5 reference; not supplied by PAM50 labels'))
atomic_json(out/'root-regions.json',dict(schema='native-cord-root-regions/v1',
    source_catalogue_sha256=digest_file(catalogue),source_surface_sha256=digest_file(base/'female-surface-v1/surface.npz'),
    reference_sha256=digest_file(ref),reference_commit=json.loads((base/'pam50-reference-v1/manifest.json').read_text())['commit'],
    root_patches=rows,registration_anchors=anchors,segments=segments,
    skeleton_changed=False,geometry_changed=False,anatomical_acceptance=False,
    limitations=['C1 has no direct cord contact in the native source.',
        'Lumbar, sacral and coccygeal boundaries need native cauda/root verification.',
        'Transverse segment planes are display partitions of overlapping root territories.'],
    elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(root_sources=len(rows),segments=len(segments),
    regions={r['name']:r['y_range_m']for r in segments if r['name']in['C3','C8','T1','T12','L1','L5','S1','S2','S5','Co1']},
    elapsed_seconds=time.monotonic()-start)),flush=True)
