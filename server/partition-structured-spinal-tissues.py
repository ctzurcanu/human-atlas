"""One disjoint field partition of a conforming native section loft.

All compartments use one common order of field cuts. Their macro identity
comes from the supplied volume mesh. No Boolean copy of an existing parent
tissue is superimposed, and no artificial all-label tie is imposed on a wall.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.score_cells import ScoreArrangement
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('domain');ap.add_argument('--output',required=True);ap.add_argument('--seconds',type=float,default=25.)
ap.add_argument('--count',type=int,default=0);ap.add_argument('--resume')
ap.add_argument('--source-tie-tolerance',type=float,default=1e-12)
args=ap.parse_args();started=time.monotonic();domain=Path(args.domain).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
m=json.loads((domain/'domain.json').read_text());z=np.load(domain/'domain.npz');p=z['positions'];t=z['tetrahedra'];tet_material=z['material'];macros=m['materials']
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';reference=base/'female-tract-fields-v1';field=np.load(reference/'fields.npz')
fm=json.loads((reference/'fields.json').read_text());plans=json.loads((reference/'plan.json').read_text())['parts']
segments=json.loads((base/'female-root-regions-v2/root-regions.json').read_text())['segments'];levels=[s['name']for s in segments];index=levels.index(m['segment'])
names=macros+['tract-'+r['name']for r in plans];probabilities=np.zeros((len(p),len(plans)),float)
middle=field['segment_y'][:30].mean(1);order=np.argsort(middle);interp=interp1d(middle[order],field['scores'][order],axis=0,bounds_error=False,
    fill_value=(field['scores'][order[0]],field['scores'][order[-1]]))
eligible={};parent_names={};n=len(field['axis']);ny=len(z['y']);nc=len(p)//ny
for k,plan in enumerate(plans):
    eligible[k]=levels.index(plan['first_level'])<=index<=levels.index(plan['last_level'])
    tract=plan['name'].rsplit('-',1)[0];side=plan['name'].rsplit('-',1)[1]
    funiculus=('posterior'if tract in ['gracile','cuneate']else'lateral'if tract in ['lateral-corticospinal','rubrospinal','anterior-spinocerebellar']else'anterior')
    parent_names[k]='white-'+funiculus+'-'+side
for layer,y in enumerate(z['y']):
    ids=np.arange(layer*nc,(layer+1)*nc);uv=(p[ids][:,[0,2]]-z['station_centers'][layer])/z['station_radii'][layer]
    query=((uv-field['axis'][0])/(field['axis'][-1]-field['axis'][0])*(n-1)).T
    values=interp(float(y)).reshape(n,n,len(plans))
    for k,plan in enumerate(plans):
        if not eligible[k]:continue
        prob=map_coordinates(values[:,:,k],query,order=1,mode='constant',cval=0.)
        probabilities[ids,k]=np.clip(prob,0,1)
allowed=[]
for i,name in enumerate(macros):allowed.append([k for k in range(len(plans))if eligible[k]and parent_names[k]==name])
background=len(plans);comparisons=[]
for candidates in allowed:
    options=[background,*candidates]
    for i,a in enumerate(options):
        for b in options[i+1:]:
            if(a,b)not in comparisons:comparisons.append((a,b))
def classify(value):
    material=int(tet_material[begin+arrangement.current_tet_index]);winner=background
    for candidate in allowed[material]:
        field,sign=arrangement.comparison_index[(winner,candidate)]
        if sign*value[field]<0:winner=candidate
    name=macros[material]if winner==background else'tract-'+plans[winner]['name']
    return name+'@'+m['segment']
arrangement=ScoreArrangement(p,np.c_[probabilities,np.full(len(p),fm['probability_floor'])],comparisons,classify,
    background_index=background,source_tie_tolerance=args.source_tie_tolerance)
inputs=dict(domain_sha256=digest_file(domain/'domain.npz'),fields_sha256=digest_file(reference/'fields.npz'),
    script_sha256=digest_file(__file__),score_algorithm_sha256=digest_file(ROOT/'geometry_tools/score_cells.py'),
    arrangement_algorithm_sha256=digest_file(ROOT/'geometry_tools/implicit_cells.py'),segment=m['segment'],source_tie_tolerance=args.source_tie_tolerance)
begin=0;previous=[]
if args.resume:
    saved=arrangement.restore_checkpoint(args.resume,inputs);begin=saved['next_tetrahedron'];previous=saved['stages']
stop=min(len(t),begin+(args.count or len(t)));assert 0<=begin<stop
def progress(number,total,points,leaves):
    if time.monotonic()-started>args.seconds:raise TimeoutError('Shared field arrangement reached its explicit saved-stage budget')
    print(json.dumps(dict(stage='shared-field-arrangement',tetrahedra=number,total=total,controls=points,leaves=leaves)),flush=True)
positions,surfaces,cells=arrangement.build(t[begin:stop],progress=progress)
stage=dict(begin=begin,end=stop,elapsed_seconds=time.monotonic()-started)
arrangement.save_checkpoint(out,dict(inputs=inputs,next_tetrahedron=stop,total_tetrahedra=len(t),stages=[*previous,stage],complete=stop==len(t)))
if stop<len(t):
    print(json.dumps(dict(segment=m['segment'],completed_tetrahedra=stop,total_tetrahedra=len(t),
        saved_stage=True,installable=False,seconds=time.monotonic()-started)),flush=True)
    raise SystemExit(0)
result=dict(positions=positions,cells=cells,evidence=dict(calculation_tetrahedra=len(t),
    comparison_fields=len(comparisons),common_field_order=True,shared_position_bank=True,
    probability_arithmetic_tolerance=1e-12,calculation_leaves=arrangement.leaves,
    reference_score_tie_regularization=args.source_tie_tolerance,
    zero_volume_contacts_not_tissue=arrangement.zero_volume_leaves,
    maximum_source_score_projection=arrangement.maximum_source_projection,
    maximum_derived_score_projection=arrangement.maximum_derived_projection,
    symbolic_contact_leaves=arrangement.symbolic_contacts,
    all_material_boundaries_from_supplied_domain=True,anatomical_acceptance=False))
attributes={};checks={}
for name,uses in cells.items():
    native=name.split('@')[0];tract=None
    if native.startswith('tract-'):
        k=next(i for i,r in enumerate(plans)if r['name']==native.removeprefix('tract-'));tract=plans[k];native=parent_names[k]
    attributes[name]=dict(native_region=native,axial_region=m['segment'],y_range_m=[float(z['y'][0]),float(z['y'][-1])])
    if tract:attributes[name].update(tract=tract['name'],ta98_term=tract['term'],probabilistic_reference=True)
    f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
    q=topology(Mesh(result['positions'],f,'model:local-female-ta98:assembled-metres'))
    checks[name]={k:q[k]for k in ['closed_oriented_manifold','open_edges','nonmanifold_edges','zero_area_faces','components','oriented_triangle_integral_m3']}
bad=[n for n,r in checks.items()if not r['closed_oriented_manifold']]
atomic_npz(out/'cells.npz',positions=result['positions'],**surfaces)
atomic_json(out/'cells.json',dict(cells=cells,surfaces={k:len(f)for k,f in surfaces.items()},attributes=attributes,segment=m['segment'],topology=checks,passed=not bad,
    construction=dict(shared_structured_volume=True,source_domain_sha256=digest_file(domain/'domain.npz'),field_sha256=digest_file(reference/'fields.npz'),
        calculation=result['evidence'],source_volume_proof=m['volume_proof'],central_canal_is_void=True,
        microscopic_borders_are_reference_estimates=True,complete_emitted_checks_pending=True),
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
assert not bad,bad
print(json.dumps(dict(segment=m['segment'],cells=len(cells),surfaces=len(surfaces),triangles=sum(len(f)for f in surfaces.values()),
    passed=True,seconds=time.monotonic()-started)),flush=True)
