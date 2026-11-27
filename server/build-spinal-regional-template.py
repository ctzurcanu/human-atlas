"""Declare one common cross-sectional graph for regional spinal tissues.

Nuclear and tract borders are anatomical reference estimates, not histological
measurements of the atlas donor. Longitudinal eligibility is retained separately
from the common mesh: an absent structure returns its tissue to its parent.
The source envelope, canal, grey/white and pial boundaries remain shared.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
import triangle,shapely
from shapely.geometry import Polygon,Point,LineString,box
from shapely.ops import polygonize
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True)
ap.add_argument('--radial-carrier',action='store_true',help='Overlay a positive radial carrier so every fine tissue triangle has one affine owner')
ap.add_argument('--field-directory',help='Source-bound probabilistic fields; defaults to the installed nine-family construction reference')
ap.add_argument('--regional-details',action='store_true',help='Add reference-guided motor pools, dorsal sensory tracts and actual white commissural tissue')
args=ap.parse_args();started=time.monotonic()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
source=base/'female-structured-template-v2';m=json.loads((source/'template.json').read_text());z=np.load(source/'template.npz')
grid=1e-4
snap=lambda p:shapely.set_precision(p,grid)
curves={n:z['curve:'+n]for n in m['curves']}
polys={n:snap(Polygon(p))for n,p in curves.items()}
macros={'pia':snap(polys['outer'].difference(polys['inner']))}
white=snap(polys['inner'].difference(polys['grey']))
occupied_grey=snap(polys['grey'].difference(polys['canal']))
for side,sign in [('left',1),('right',-1)]:
    half=box(0,-2,2,2)if sign>0 else box(-2,-2,0,2)
    for region,strip in [('anterior',box(-2,.15,2,2)),('posterior',box(-2,-2,2,-.18)),('intermediate',box(-2,-.18,2,.15))]:
        macros['grey-'+region+'-'+side]=snap(occupied_grey.intersection(half).intersection(strip))
    anterior=Polygon([(0,.13),(2,1.93),(2,2),(0,2)])
    posterior=Polygon([(0,-.16),(2,-1.96),(2,-2),(0,-2)])
    if sign<0:
        anterior=Polygon([(-x,y)for x,y in anterior.exterior.coords]);posterior=Polygon([(-x,y)for x,y in posterior.exterior.coords])
    macros['white-anterior-'+side]=snap(white.intersection(half).intersection(anterior))
    macros['white-posterior-'+side]=snap(white.intersection(half).intersection(posterior))
    macros['white-lateral-'+side]=snap(white.intersection(half).difference(shapely.union_all([anterior,posterior])))
levels=[*[f'C{i}'for i in range(1,9)],*[f'T{i}'for i in range(1,13)],
        *[f'L{i}'for i in range(1,6)],*[f'S{i}'for i in range(1,6)],'Co1']
features=[]
def feature(name,term,shape,*,owner,first='C1',last='S5',parents=(),priority=0,source_ref='UTHealth',level_ranges=None):
    assert levels.index(first)<=levels.index(last)
    shape=snap(shape.intersection(polys['inner']))
    if not shape.is_empty:
        features.append(dict(name=name,term=term,shape=shape,first=first,last=last,
            parents=list(parents),priority=priority,source_ref=source_ref,owner=owner,
            **(dict(level_ranges=level_ranges)if level_ranges else{})))
def ellipse(x,z,rx,rz):
    a=np.linspace(0,2*np.pi,24,endpoint=False)
    return Polygon(np.c_[x+rx*np.cos(a),z+rz*np.sin(a)])
grey=polys['grey']
# Transverse laminae follow the dorsal horn from apex to base; laminae VII/X
# and motor nuclei have their own shared territories, rather than point tags.
for side,sign in [('left',1),('right',-1)]:
    posterior=macros['grey-posterior-'+side]
    for name,term,lo,hi,parents in [
        ('marginal','117',-2,-.72,['116']),
        ('gelatinous','119',-.72,-.58,['118']),
        ('proprius','121',-.58,-.39,['120']),
        ('lamina-V','122',-.39,-.25,['120']),
        ('lamina-VI','124',-.25,2,['123'])]:
        feature(name+'-'+side,'A14.1.02.'+term,posterior.intersection(box(-2,lo,2,hi)),
            owner='grey-posterior-'+side,parents=['A14.1.02.'+p for p in parents],priority=20)
    intermediate=macros['grey-intermediate-'+side]
    central=intermediate.intersection(box(-.145,-.115,.145,.115))
    feature('grey-commissure-anterior-'+side,'A14.1.02.303',central.intersection(box(-2,0,2,2)),
        owner='grey-intermediate-'+side,parents=['A14.1.02.025','A14.1.02.302','A14.1.02.134'],priority=50)
    feature('grey-commissure-posterior-'+side,'A14.1.02.304',central.intersection(box(-2,-2,2,0)),
        owner='grey-intermediate-'+side,parents=['A14.1.02.025','A14.1.02.302','A14.1.02.134'],priority=50)
    # Small bilateral nuclear territories lie in the intermediate grey, with
    # eligibility independent of their reference cross-sectional footprint.
    for name,term,x,zz,rx,rz,first,last in [
        ('clarke','135',.235,-.137,.055,.042,'C8','L3'),
        ('intermediomedial','137',.175,.025,.035,.040,'C1','S5'),
        ('intermediolateral','133',.373,.005,.045,.062,'T1','L2'),
        ('sacral-parasympathetic','138',.343,.015,.055,.050,'S2','S4')]:
        feature(name+'-'+side,'A14.1.02.'+term,
            ellipse(sign*x,zz,rx,rz).intersection(intermediate),owner='grey-intermediate-'+side,first=first,last=last,
            parents=['A14.1.02.131'],priority=60,
            source_ref='UTHealth; human thoracic histology (Bergman et al.)')
    # This is a lamina VII territory, not a new nucleus covering the whole horn.
    feature('lateral-horn-'+side,'A14.1.02.132',
        intermediate.intersection(box(.30 if sign>0 else -2,-.16,2 if sign>0 else -.30,.15)),
        owner='grey-intermediate-'+side,first='T1',last='L2',parents=['A14.1.02.022','A14.1.02.131'],priority=25)
    feature('lateral-intermediate-'+side,'A14.1.02.136',intermediate.difference(box(-.23,-2,.23,2)),
        owner='grey-intermediate-'+side,parents=['A14.1.02.131'],priority=5)
    # Somatotopic motor pools are deliberately not invented here. The phrenic
    # and accessory territories have independently established level limits;
    # the remaining individual motor nuclei need their own construction data.
    anterior=macros['grey-anterior-'+side]
    feature('accessory-nucleus-'+side,'A14.1.02.112',
        ellipse(sign*.495,.225,.052,.045).intersection(anterior),
        owner='grey-anterior-'+side,first='C1',last='C5',parents=['A14.1.02.104'],priority=80)
    feature('phrenic-nucleus-'+side,'A14.1.02.113',
        ellipse(sign*.315,.292,.045,.054).intersection(anterior),
        owner='grey-anterior-'+side,first='C3',last='C5',parents=['A14.1.02.104'],priority=80)
    # Perigray association-fibre territories are finite tissue, not outlines.
    ring=snap(grey.buffer(.045,quad_segs=8).difference(grey))
    for fun,term in [('anterior','203'),('lateral','216'),('posterior','240')]:
        feature(fun+'-proprius-'+side,'A14.1.02.'+term,
            ring.intersection(macros['white-'+fun+'-'+side]),owner='white-'+fun+'-'+side,priority=100)

    if args.regional_details:
        # Independent, deliberately simplified territories follow the human
        # teaching topography, rather than tracing or copying the source image.
        # The two medial pools are continuous/axial or thoracic respectively;
        # lateral pools have separate cervical and lumbosacral eligibility.
        # These small volumes represent reference territories, never measured
        # donor cell columns or their exact cytoarchitectural borders.
        motor_reference='https://www.nan.upol.cz/neuro/cd43_en.html; https://open.umich.edu/sites/default/files/downloads/030909.m1-cns.spinalcordpti-handout.pdf'
        for name,term,x,zz,rx,rz,ranges in [
            ('anteromedial-motor','107',.332,.409,.035,.038,[['C1','S5']]),
            ('posteromedial-motor','110',.264,.237,.029,.034,[['T1','L2']]),
            ('anterolateral-motor','105',.497,.322,.035,.035,[['C5','T1'],['L2','S2']]),
            ('posterolateral-motor','108',.474,.221,.030,.026,[['C6','T1'],['L3','S3']]),
            ('retroposterolateral-motor','109',.424,.171,.026,.018,[['C8','T1'],['S1','S2']]),
            ('central-motor','111',.411,.344,.034,.032,[['C5','T1'],['L2','S2']])]:
            feature(name+'-'+side,'A14.1.02.'+term,
                ellipse(sign*x,zz,rx,rz).intersection(anterior),owner='grey-anterior-'+side,
                first=ranges[0][0],last=ranges[-1][1],level_ranges=ranges,priority=75,
                parents=['A14.1.02.104'],source_ref=motor_reference)

        # The dorsal spinocerebellar territory is superficial to the lateral
        # corticospinal field. Lissauer is immediately dorsal to the horn apex;
        # neither is obtained by relabelling the combined PAM spinal lemniscus.
        tract_reference='https://nba.uth.tmc.edu/neuroscience/m/s2/chapter03.html; https://dspace.uzhnu.edu.ua/jspui/bitstream/lib/49081/3/El-Variant-Poslidovno2018C-M-Y-K.pdf'
        lateral=macros['white-lateral-'+side]
        rim=polys['inner'].difference(polys['inner'].buffer(-.105,quad_segs=8))
        posterior_rim=rim.intersection(box(.38 if sign>0 else -2,-2,2 if sign>0 else -.38,-.19))
        feature('posterior-spinocerebellar-'+side,'A14.1.02.227',posterior_rim.intersection(lateral),
            owner='white-lateral-'+side,first='C1',last='L3',priority=70,source_ref=tract_reference)
        for fun in ['posterior','lateral']:
            feature('lissauer-'+fun+'-'+side,'A14.1.02.228',
                ellipse(sign*.466,-.789,.078,.072).intersection(macros['white-'+fun+'-'+side]),
                owner='white-'+fun+'-'+side,priority=110,source_ref=tract_reference)

        # A white commissure is occupied white matter immediately adjacent to
        # the grey bridge, not the canal or a marker placed at the centre.
        for name,term,shape in [
            ('anterior','305',box(-.145,.107,.145,.194)),
            ('posterior','306',box(-.115,-.197,.115,-.125))]:
            for fun in ['anterior','lateral','posterior']:
                owner='white-'+fun+'-'+side
                feature('white-commissure-'+name+'-'+fun+'-'+side,'A14.1.02.'+term,
                    shape.intersection(macros[owner]),owner=owner,priority=120,
                    parents=['A14.1.02.301'],source_ref='https://nba.uth.tmc.edu/neuroscience/m/s2/chapter03.html')

# Source-bound PAM50 probabilities are clipped as affine fields in a common
# triangular grid. All competing fields use the same underlying samples.
field=Path(args.field_directory).resolve()if args.field_directory else base/'female-tract-fields-v1';fm=json.loads((field/'fields.json').read_text());fz=np.load(field/'fields.npz')
plans=json.loads((field/'plan.json').read_text())['parts'];axis=fz['axis'];n=len(axis)
xy=fz['normalized_controls'];scores=fz['scores'][3];floor=float(fm['probability_floor'])
grid_faces=[]
for i in range(n-1):
    for j in range(n-1):
        a=i*n+j;b=a+n;c=a+1;d=b+1;grid_faces.extend([[a,b,d],[a,d,c]])
def clip(p,v,k):
    result=[];values=[]
    for i,j in zip(range(len(p)),np.roll(np.arange(len(p)),-1)):
        a,b=v[i,k],v[j,k]
        if a>=0:result.append(p[i]);values.append(v[i])
        if(a>=0)!=(b>=0):
            w=a/(a-b);result.append(p[i]*(1-w)+p[j]*w);values.append(v[i]*(1-w)+v[j]*w)
    return np.asarray(result),np.asarray(values)
for macro,parent in macros.items():
    if not macro.startswith('white-'):continue
    fun,side=macro.split('-')[1:];eligible=[]
    for k,plan in enumerate(plans):
        tract,s=plan['name'].rsplit('-',1)
        f=plan.get('funiculus')or('posterior'if tract in ['gracile','cuneate']else'lateral'if tract in ['lateral-corticospinal','rubrospinal','anterior-spinocerebellar']else'anterior')
        if s==side and f==fun:eligible.append(k)
    for k in eligible:
        pieces=[]
        for face in grid_faces:
            vv=scores[face,k,None]-np.c_[scores[np.ix_(face,[q for q in eligible if q!=k])],np.full(3,floor)]
            p=xy[face].copy()
            for column in range(vv.shape[1]):
                if not len(p):break
                p,vv=clip(p,vv,column)
            if len(p)>=3:
                poly=snap(Polygon(p))
                if poly.is_valid and poly.area>0:pieces.append(poly)
        shape=snap(shapely.union_all(pieces,grid_size=grid).intersection(parent))
        plan=plans[k]
        feature('tract-'+plan['name'],plan['term'],shape,owner=macro,first=plan['first_level'],last=plan['last_level'],
            parents=plan.get('parent_terms',[]),priority=40,source_ref='PAM50 commit '+fm['reference_commit']+'; independent anatomical extent')

lines=[p.boundary for p in macros.values()]
lines.extend(f['shape'].boundary for f in features)
graph=shapely.union_all([snap(line)for line in lines],grid_size=grid)
carrier_points=[];carrier_faces=[];carrier_polys=[]
radial_fractions=np.array([.04,.1,.25,.45,.65,.82,1.])
if args.radial_carrier:
    # The measured cord outline is star shaped; the H-shaped grey matter need
    # not be. Mapping a common carrier moves every embedded tissue together.
    # Its boundaries are overlaid once, so fine triangles cannot straddle two
    # different affine maps. Do not round the resulting line intersections.
    inner=np.rint(curves['inner']/grid)*grid;outer=np.rint(curves['outer']/grid)*grid;n=len(inner)
    assert n==len(outer)
    carrier_points=np.vstack([[[0.,0.]],*[t*inner for t in radial_fractions],outer])
    for j in range(n):carrier_faces.append([0,1+j,1+(j+1)%n])
    for k in range(len(radial_fractions)):
        a=1+k*n;b=a+n
        for j in range(n):
            jj=(j+1)%n;carrier_faces.extend([[a+j,b+j,b+jj],[a+j,b+jj,a+jj]])
    carrier_faces=np.array(carrier_faces,int);q=carrier_points[carrier_faces]
    cross=(q[:,1,0]-q[:,0,0])*(q[:,2,1]-q[:,0,1])-(q[:,1,1]-q[:,0,1])*(q[:,2,0]-q[:,0,0])
    assert np.all(cross>0),'The reference carrier must be positive before partitioning'
    carrier_polys=[Polygon(p)for p in q]
    graph=shapely.union_all([graph,*[p.boundary for p in carrier_polys]])
    carrier_tree=shapely.STRtree(carrier_polys)
regions=[];holes=[]
for poly in polygonize(graph):
    p=poly.representative_point()
    if not polys['outer'].covers(p):continue
    if polys['canal'].covers(p):holes.append([p.x,p.y]);continue
    # The common graph owns the rounded borders. Independently rounded unions
    # of thousands of source triangles may leave sub-grid strips; do not use
    # those separate unions as a second, conflicting ownership authority.
    xx,zz=p.x,p.y;side='left'if xx>0 else'right'
    if not polys['inner'].covers(p):macro='pia'
    elif polys['grey'].covers(p):
        macro='grey-'+('anterior'if zz>.15 else'posterior'if zz<-.18 else'intermediate')+'-'+side
    else:
        macro='white-'+('anterior'if zz>.13+.9*abs(xx)else'posterior'if zz<-.16-.9*abs(xx)else'lateral')+'-'+side
    candidates=[i for i,f in enumerate(features)if f['owner']==macro and f['shape'].covers(p)]
    carrier=None
    if args.radial_carrier:
        owners=carrier_tree.query(p,predicate='within');assert len(owners)==1,'A common region needs one affine carrier owner'
        carrier=int(owners[0])
    regions.append((macro,sorted(candidates,key=lambda i:(-features[i]['priority'],i)),poly,carrier))
points=[];lookup={};edges=set();seeds=[]
def node(q):
    key=tuple(map(float,q))
    if key not in lookup:lookup[key]=len(points);points.append(key)
    return lookup[key]
for i,(macro,candidates,p,carrier)in enumerate(regions):
    for ring in[p.exterior,*p.interiors]:
        ids=[node(q)for q in list(ring.coords)[:-1]]
        for a,b in zip(ids,ids[1:]+ids[:1]):edges.add(tuple(sorted((a,b))))
    p=p.representative_point();seeds.append([p.x,p.y,i,.025])
assert holes and regions, 'The canal needs an explicit void seed inside its shared boundary'
result=triangle.triangulate(dict(vertices=np.asarray(points),segments=np.array(sorted(edges)),
    holes=holes,regions=np.array(seeds)), 'pa0.025A')
uv=np.asarray(result['vertices'],float);faces=np.asarray(result['triangles'],int)
assert len(uv)<=(12000 if args.radial_carrier else 6000), 'The common section exceeds its declared control memory budget'
labels=np.rint(result['triangle_attributes'].ravel()).astype(int)
q=uv[faces];a=q[:,1]-q[:,0];b=q[:,2]-q[:,0];area=a[:,0]*b[:,1]-a[:,1]*b[:,0]
negative=area<0;faces[negative]=faces[negative][:,[0,2,1]]
assert np.all(abs(area)>0)
names=m['curves'];kind=np.full(len(uv),-1,int);edge=np.full(len(uv),-1,int);fraction=np.zeros(len(uv));residual=np.full(len(uv),np.inf)
# Only controls on a physical interface are fixed. A nearby nuclear or tract
# control must not be projected onto that interface merely by proximity.
incident=defaultdict(list)
for i,f in enumerate(faces):
    for a,b in zip(f,np.roll(f,-1)):incident[tuple(sorted((int(a),int(b))))].append(i)
physical=defaultdict(set)
for(a,b),owners in incident.items():
    macros0=[regions[int(labels[i])][0]for i in owners]
    if len(owners)==1:
        point=Point(uv[[a,b]].mean(0));name=min(['outer','canal'],key=lambda k:polys[k].boundary.distance(point))
    elif ('pia'in macros0)and len(set(macros0))>1:name='inner'
    elif {n.split('-')[0]for n in macros0}=={'grey','white'}:name='grey'
    else:continue
    physical[name].update([a,b])
for number,name in enumerate(names):
    poly=curves[name];a=poly;b=np.roll(poly,-1,axis=0);e=b-a;delta=uv[:,None]-a
    t=np.clip(np.einsum('nki,ki->nk',delta,e)/np.einsum('ki,ki->k',e,e),0,1)
    distance=np.linalg.norm(delta-t[:,:,None]*e,axis=2);nearest=distance.argmin(1);dist=distance[np.arange(len(uv)),nearest]
    possible=np.zeros(len(uv),bool);possible[list(physical[name])]=True
    assert np.all(dist[possible]<grid*2), 'A physical interface lost its declared source correspondence'
    keep=possible&(dist<residual)
    kind[keep]=number;edge[keep]=nearest[keep];fraction[keep]=t[np.flatnonzero(keep),nearest[keep]];residual[keep]=dist[keep]
fixed=np.flatnonzero(kind>=0)
extra={}
if args.radial_carrier:
    owners=np.array([regions[i][3]for i in labels],int)
    # Choose one incident carrier per vertex and prove agreement in every
    # incident triangle. Shared edges use the identical position word.
    vertex_owner=np.full(len(uv),-1,int)
    for tri,owner in zip(faces,owners):vertex_owner[tri]=owner
    assert np.all(vertex_owner>=0)
    q=carrier_points[carrier_faces[vertex_owner]]
    matrix=np.stack([q[:,1]-q[:,0],q[:,2]-q[:,0]],axis=2)
    vw=np.linalg.solve(matrix,(uv-q[:,0])[...,None])[...,0]
    bary=np.c_[1-vw.sum(1),vw]
    assert bary.min()>-1e-8,'A fine control left its reference carrier'
    extra=dict(carrier_positions=carrier_points,carrier_faces=carrier_faces,
        carrier_owner=vertex_owner,carrier_barycentric=bary,triangle_carrier=owners,radial_fractions=radial_fractions)
atomic_npz(out/'template.npz',positions=uv,faces=faces,region=labels,
    material=np.array([m['materials'].index(regions[i][0])for i in labels]),
    angles=z['angles'],fixed_controls=fixed,curve_kind=kind,curve_edge=edge,curve_fraction=fraction,
    reference_center=z['reference_center'],reference_radii=z['reference_radii'],
    **{'curve:'+n:p for n,p in curves.items()},**extra)
atomic_json(out/'template.json',dict(passed=True,materials=m['materials'],curves=names,
    features=[{k:v for k,v in f.items()if k!='shape'}for f in features],
    regions=[dict(macro=macro,candidates=c,carrier=carrier)for macro,c,p,carrier in regions],
    source_template_sha256=digest_file(source/'template.npz'),
    field_sha256=digest_file(field/'fields.npz'),positions=len(uv),triangles=len(faces),
    reference_arithmetic_grid=grid,maximum_profile_rounding_m=float(grid*np.max(z['reference_radii'])),
    macro_ownership_from_one_declared_common_graph=True,
    finite_regions_not_points=True,nuclear_borders_are_anatomical_estimates=True,
    common_positive_radial_carrier=bool(args.radial_carrier),
    unestablished_motor_pool_territories_not_invented=not args.regional_details,
    motor_pool_borders_and_longitudinal_extents_are_reference_estimates=bool(args.regional_details),
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(passed=True,features=len(features),regions=len(regions),positions=len(uv),
    triangles=len(faces),seconds=time.monotonic()-started)),flush=True)
