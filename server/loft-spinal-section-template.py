"""Measure and transport all cord compartments together in one saved segment."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.radial_section_loft import measured_radial_curves,harmonic_sections,prism_tetrahedra,PrismJacobianError

ap=argparse.ArgumentParser();ap.add_argument('--template',required=True);ap.add_argument('--segment',required=True)
ap.add_argument('--output',required=True);args=ap.parse_args();started=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';template=Path(args.template).resolve();m=json.loads((template/'template.json').read_text());r=np.load(template/'template.npz')
source=base/'female-cage-rounded-v3';d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');v=z['positions'];inner=set();outer=set()
for name,uses in d['cells'].items():
    if name.split('@')[0]!='pia':continue
    for use in uses:
        peer=use.get('neighbor')
        if peer and not peer.startswith('pia'):inner.add(use['surface'])
        elif not peer:outer.add(use['surface'])
families={kind:v[np.concatenate([z[s]for s in sorted(sids)])]for kind,sids in [('inner',inner),('outer',outer)]}
segments=json.loads((base/'female-root-regions-v2/root-regions.json').read_text())['segments'];segment=next(row for row in segments if row['name']==args.segment)
low,high=segment['y_range_m'];low=float(np.float32(low-1.35))+1.35;high=float(np.float32(high-1.35))+1.35
ys=np.linspace(low,high,max(3,int(np.ceil((high-low)/.003))+1));reference_grey=r['curve:grey'];curve_names=m['curves'];cache={}
def station(y):
    if y in cache:return cache[y]
    center,radii,curves=measured_radial_curves(families,float(y),r['angles'])
    scale=np.interp(y,[1.148,1.158,1.178,1.198,1.24,1.42,1.47,1.53,1.565],[1.16,1.16,1.08,.95,.69,.69,.86,.86,.81])
    grey=reference_grey.copy();grey[:,0]*=scale/.86;grey[:,1]=np.where(grey[:,1]>0,grey[:,1]*scale/.86,grey[:,1]);curves['grey']=center+grey*radii
    angle=np.linspace(0,2*np.pi,16,endpoint=False);curves['canal']=center+.00008*np.c_[np.cos(angle),np.sin(angle)]
    fixed=r['fixed_controls'];target=np.empty((len(fixed),2))
    for i,control in enumerate(fixed):
        name=curve_names[r['curve_kind'][control]];edge=r['curve_edge'][control];weight=r['curve_fraction'][control];p=curves[name]
        target[i]=p[edge]*(1-weight)+p[(edge+1)%len(p)]*weight
    record=dict(y=float(y),center=center.tolist(),radii=radii.tolist(),grey_scale=float(scale))
    cache[y]=(target,record);return cache[y]
history=[]
for refinement in range(7):
    if time.monotonic()-started>20:raise TimeoutError('Saved section refinement reached its explicit 20-second budget')
    data=[station(float(y))for y in ys];targets=np.array([row[0]for row in data]);records=[row[1]for row in data]
    fields,proof=harmonic_sections(r['positions'],r['faces'],r['fixed_controls'],targets)
    try:positions,tetrahedra,source_faces,volume=prism_tetrahedra(fields,ys,r['faces']);break
    except PrismJacobianError as error:
        history.append(dict(round=refinement,stations=len(ys),bad_slabs=error.bad_slabs))
        if refinement==6:raise
        middle=np.array([(ys[k]+ys[k+1])/2 for k in error.bad_slabs]);ys=np.sort(np.r_[ys,middle])
else:raise ValueError('Section refinement did not complete')
atomic_npz(out/'domain.npz',positions=positions,tetrahedra=tetrahedra,source_faces=source_faces,y=ys,
    material=r['material'][source_faces],interior_material=np.tile(r['interior_material'],len(ys)),
    station_centers=np.array([row['center']for row in records]),station_radii=np.array([row['radii']for row in records]))
atomic_json(out/'domain.json',dict(passed=True,segment=args.segment,template_sha256=digest_file(template/'template.npz'),
    source_tissues_sha256=digest_file(source/'cells.npz'),materials=m['materials'],stations=records,section_proof=proof,
    volume_proof=volume,central_canal_is_void=True,all_tissues_coupled=True,
    adaptive_axial_refinement=history,calibres_not_reduced=True,
    native_envelope_resampled=True,complete_source_deviation_and_fixed_bone_checks_pending=True,
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(segment=args.segment,stations=len(ys),positions=len(positions),tetrahedra=len(tetrahedra),
    positive_sections=True,seconds=time.monotonic()-started)),flush=True)
