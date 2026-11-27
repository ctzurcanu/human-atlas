"""Diagnostic blood-space flood fill with native, computation-only root gates.

No gate is a tissue solid, display triangle, picking target or installed cap.
The voxel result supplements actual root sections and does not certify them.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.ndimage import label
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.tube_continuity import boundary_loops
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';out=base/'blood-space-v12';out.mkdir(exist_ok=False);z=np.load(base/'whole-wall-v7/field.npz');d=z['distance']+.00001;o=z['origin'];h=float(z['spacing']);shape=np.array(d.shape);free=d>0;path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();gates=[];atrial_center=np.array([.021,1.40,-.004]);whole_center=np.array([.02,1.38,.03])
ids=['ZA:Ascending aorta','ZA:Pulmonary trunk','ZA:Superior vena cava','ZA:Inferior vena cava (thoracic part)','ZA:Left inferior pulmonary vein','ZA:Left superior pulmonary vein','ZA:Right inferior pulmonary vein','ZA:Right superior pulmonary vein']
for pid in ids:
 m=reader.read_part(a,parts[pid],frame='model:local-ta98:assembled-metres',catalogue_path=path);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);rings=boundary_loops(inv[m.faces]);ref=atrial_center if'pulmonary vein'in pid else whole_center;ring=min(rings,key=lambda q:np.linalg.norm(v[q].mean(0)-ref));p=v[ring];c=p.mean(0);_,s,vh=np.linalg.svd(p-c,full_matrices=False);n=vh[-1];radius=np.linalg.norm(p-c,axis=1).max()+.0015;extent=radius+2*h;lo=np.maximum(0,np.floor((c-extent-o)/h).astype(int));hi=np.minimum(shape,np.ceil((c+extent-o)/h).astype(int)+1);indices=np.indices(tuple(hi-lo)).reshape(3,-1).T+lo;points=o+indices*h;delta=points-c;axial=delta@n;radial=np.linalg.norm(delta-axial[:,None]*n,axis=1);mask=(np.abs(axial)<1.5*h)&(radial<=radius);chosen=indices[mask];free[tuple(chosen.T)]=False;gates.append(dict(part_id=pid,center=c.tolist(),normal=n.tolist(),radius_m=float(radius),half_thickness_m=1.5*h,source_controls=p.tolist(),rendered=False))
regions,n=label(free);counts=np.bincount(regions.ravel());outside=int(regions[0,0,0]);seeds={'right-atrium':[-.022,1.39,.022],'right-ventricle':[.025,1.36,.06],'left-atrium':[.018,1.405,-.012],'left-ventricle':[.05,1.36,.025]};observed={}
for key,point in seeds.items():
 ijk=np.rint((np.array(point)-o)/h).astype(int);region=int(regions[tuple(ijk)]);observed[key]=dict(point=point,region=region,outside=region==outside,occupied_tissue=region==0,voxel_volume_mL=float(counts[region]*h**3*1e6))
result=dict(computed=True,gates=gates,seeds=observed,outside_region=outside,positive_region_count=n,compartment_volumes_mL=sorted((counts[1:]*h**3*1e6).tolist(),reverse=True)[:20],elapsed_seconds=time.monotonic()-start,diagnostic_only=True,anatomical_acceptance=False)
atomic_npz(out/'regions.npz',regions=regions,origin=o,spacing=np.array(h));atomic_json(out/'compartments.json',result);print(json.dumps(result),flush=True)
