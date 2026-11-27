"""Transport a common regional graph along the fitted female cord.

The prism orientation gate checks every physical interpolation element, while
only the exterior/reciprocal tissue faces are retained. Axial interfaces use
the identical graph at each end; no independent segment-cap triangulation is
permitted. Oblique terminal sections are measured in their own source planes.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time,inspect,hashlib
import numpy as np
from shapely.geometry import Polygon,Point
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.radial_section_loft import measured_radial_curves,harmonic_sections,bound_section_deformations,transport_section_sequence,radial_carrier_section,carrier_prism_edge_grid
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('--template',required=True)
ap.add_argument('--output',required=True);ap.add_argument('--seconds',type=float,default=45.)
ap.add_argument('--segment',help='Bounded regional probe; omit for the continuous cord')
ap.add_argument('--section-fit-limit',type=float,default=0.,help='Explicit maximum deviation from measured target controls, in metres; zero stops on any section inversion')
ap.add_argument('--regional-profile',action='store_true',help='Keep reference grey inside the measured cord with a 10-micrometre white rim and taper the conus canal inside pericanal grey')
ap.add_argument('--successive-transport',action='store_true')
ap.add_argument('--interior-section-repair-bound',type=float,default=0.)
ap.add_argument('--station-cache',help='Reuse a completed, source/hash-bound profile stage; not an installable assembly')
ap.add_argument('--emit-stage',action='store_true',help='Save completed construction before separately supervised topology batches; never marks the stage passed')
ap.add_argument('--measured-conus-profile',action='store_true',help='Construct shorter distal posterior grey horns inside the measured envelope; retain at least 65% of their reference length. Does not change anterior-grey width, relax any outer-cord fit or waive stored-geometry checks')
args=ap.parse_args();started=time.monotonic()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
template=Path(args.template).resolve();m=json.loads((template/'template.json').read_text())
with np.load(template/'template.npz')as archive:r={key:archive[key]for key in archive.files}
source=base/'female-cage-rounded-v3';d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');v=z['positions']
ids={'inner':set(),'outer':set()};exterior=set()
for name,uses in d['cells'].items():
    for u in uses:
        if not u.get('neighbor'):exterior.add((u['surface'],u['side']))
        if name.split('@')[0]!='pia':continue
        peer=u.get('neighbor')
        if peer and not peer.startswith('pia'):ids['inner'].add(u['surface'])
        elif not peer:ids['outer'].add(u['surface'])
families={kind:v[np.concatenate([z[s]for s in sorted(sids)])]for kind,sids in ids.items()}
ef=np.concatenate([z[s]if sign==1 else z[s][:,[0,2,1]]for s,sign in exterior]);eq=v[ef]
en=np.cross(eq[:,1]-eq[:,0],eq[:,2]-eq[:,0]);en/=np.linalg.norm(en,axis=1)[:,None];ec=eq.mean(1)
segments=json.loads((base/'female-root-regions-v2/root-regions.json').read_text())['segments']
levels=[s['name']for s in segments];end_sections={}
if not args.segment:
    # The true caps are oblique and slightly nonplanar. A nearby complete
    # source section is used, with the discarded terminal distance explicit.
    for label,mask,sign in [('upper',(ec[:,1]>1.56)&(en[:,1]>.5),-1),
                            ('lower',(ec[:,1]<1.148)&(en[:,1]<-.5),1)]:
        p=v[np.unique(ef[mask])];coef=np.linalg.lstsq(np.c_[p[:,0],p[:,2],np.ones(len(p))],p[:,1],rcond=None)[0]
        transformed={k:q.copy()for k,q in families.items()}
        for q in transformed.values():q[:,:,1]-=coef[0]*q[:,:,0]+coef[1]*q[:,:,2]
        cap_y=p[:,1]-coef[0]*p[:,0]-coef[1]*p[:,2]
        plane=float(cap_y.min()-.00002 if sign<0 else cap_y.max()+.00002)
        center,radii,curves=measured_radial_curves(transformed,plane,r['angles'])
        y=float(plane+coef[0]*center[0]+coef[1]*center[1])
        end_sections[label]=dict(y=y,center=center,radii=radii,curves=curves,
            slope=coef[:2],plane=plane,maximum_source_cap_offset_m=float(abs(cap_y-plane).max()))
    low,high=end_sections['lower']['y'],end_sections['upper']['y']
else:
    s=next(s for s in segments if s['name']==args.segment);low,high=s['y_range_m']
    low=float(np.float32(low-1.35))+1.35;high=float(np.float32(high-1.35))+1.35
bounds=sorted({low,high,*[float(np.float32(y-1.35))+1.35 for s in segments for y in s['y_range_m']if low<y<high]})
ys=np.unique(np.concatenate([np.linspace(a,b,max(2,int(np.ceil((b-a)/.003))+1))for a,b in zip(bounds,bounds[1:])]))
cache={};reference_grey=r['curve:grey'];fixed=r['fixed_controls']
radial_carrier=bool(m.get('common_positive_radial_carrier'))
if radial_carrier:
    grey_faces=np.array([m['materials'][int(i)].startswith('grey-')for i in r['material']])
    q=r['positions'][r['faces']];a=q[:,1]-q[:,0];b=q[:,2]-q[:,0]
    area=a[:,0]*b[:,1]-a[:,1]*b[:,0]
    reference_grey_fraction=float(.5*area[grey_faces].sum()/Polygon(r['curve:inner']).area)
def station(y):
    if y in cache:return cache[y]
    end=next((e for e in end_sections.values()if e['y']==y),None)
    if end:center,radii,curves=end['center'],end['radii'],dict(end['curves'])
    else:center,radii,curves=measured_radial_curves(families,y,r['angles'])
    scale=float(np.interp(y,[1.148,1.158,1.178,1.198,1.24,1.42,1.47,1.53,1.565],
        [1.16,1.16,1.08,.95,.69,.69,.86,.86,.81]))
    if radial_carrier:
        # The previous explicit H-profile can pinch white matter against the
        # conus. Regional grey proportion is carried through the whole graph
        # instead; neither donor surface is pushed to make room for it.
        cache[y]=(np.empty((0,2)),dict(y=y,center=center.tolist(),radii=radii.tolist(),grey_scale=scale,
            inner_profile=curves['inner'].tolist(),outer_profile=curves['outer'].tolist(),
            desired_grey_fraction=reference_grey_fraction*(scale/.86)**2,
            outer_inner_pial_profiles_literal=True,grey_borders_are_reference_estimates=True))
        return cache[y]
    grey=reference_grey.copy();grey[:,0]*=scale/.86;grey[:,1]=np.where(grey[:,1]>0,grey[:,1]*scale/.86,grey[:,1])
    curves['grey']=center+grey*radii
    grey_correction=0.;grey_posterior_fraction=1.;canal_radius=.00008
    if args.regional_profile:
        occupied=Polygon(curves['inner']).buffer(-.000010)
        assert occupied.is_valid and not occupied.is_empty
        before=curves['grey'].copy()
        if not occupied.covers(Polygon(before)):
            delta=before-center;lower,upper=.1,1.;posterior_base=-.18*radii[1]
            def shortened(fraction):
                probe=before.copy();negative=delta[:,1]<posterior_base
                probe[negative,1]=center[1]+posterior_base+(delta[negative,1]-posterior_base)*fraction
                return probe
            probe=shortened(lower)
            if not occupied.covers(Polygon(probe)):raise ValueError('The regional grey profile needs a different shape, not further posterior shortening')
            for _ in range(42):
                middle=(lower+upper)/2;probe=shortened(middle)
                if occupied.covers(Polygon(probe)):lower=middle
                else:upper=middle
            curves['grey']=shortened(lower)
            grey_posterior_fraction=lower
        grey_correction=float(np.linalg.norm(curves['grey']-before,axis=1).max())
        if args.measured_conus_profile:
            if grey_posterior_fraction<.65:
                atomic_npz(out/'rejected-profile.npz',inner=curves['inner'],outer=curves['outer'],
                    reference_grey=before,fitted_grey=curves['grey'],center=center,radii=radii)
                atomic_json(out/'rejected-profile.json',dict(y_m=y,grey_reference_posterior_length_fraction=grey_posterior_fraction,
                    maximum_grey_reference_change_m=grey_correction,minimum_permitted_posterior_length_fraction=.65,
                    outer_profiles_unchanged=True,anatomical_acceptance=False,installed=False))
                raise ValueError('The measured grey profile needs a different shape at y='+str(y)+
                    '; retained posterior length fraction='+str(grey_posterior_fraction))
        elif grey_correction>.0001:raise ValueError('Regional grey correction exceeds its declared 0.1 mm bound')
        gp=Polygon(curves['grey']);assert gp.is_valid and Polygon(curves['inner']).covers(gp)
        canal_radius=min(canal_radius,float(gp.boundary.distance(Point(center)))*.5)
        if canal_radius<.000010:raise ValueError('The conus canal has not retained its explicit 10-micrometre minimum radius')
    angle=np.linspace(0,2*np.pi,16,endpoint=False);curves['canal']=center+canal_radius*np.c_[np.cos(angle),np.sin(angle)]
    target=np.empty((len(fixed),2))
    for i,control in enumerate(fixed):
        curve=curves[m['curves'][r['curve_kind'][control]]];edge=r['curve_edge'][control];w=r['curve_fraction'][control]
        target[i]=curve[edge]*(1-w)+curve[(edge+1)%len(curve)]*w
        original=r['curve:'+m['curves'][r['curve_kind'][control]]]
        original_point=original[edge]*(1-w)+original[(edge+1)%len(original)]*w
        target[i]+=(r['positions'][control]-original_point)*radii
    cache[y]=(target,dict(y=y,center=center.tolist(),radii=radii.tolist(),grey_scale=scale,
        maximum_grey_boundary_correction_m=grey_correction,grey_reference_posterior_length_fraction=grey_posterior_fraction,
        outer_inner_pial_profiles_literal=True,canal_radius_m=canal_radius));return cache[y]
face=np.sort(r['faces'],axis=1);reference=r['positions'][face]
a=reference[:,1]-reference[:,0];b=reference[:,2]-reference[:,0]
area=a[:,0]*b[:,1]-a[:,1]*b[:,0];expected=np.r_[-np.sign(area),np.sign(area),-np.sign(area)]
history=[]
binding=dict(template_sha256=digest_file(template/'template.npz'),source_sha256=digest_file(source/'cells.npz'),
    profile_code_sha256=hashlib.sha256(inspect.getsource(station).encode()).hexdigest(),
    regional_profile=args.regional_profile,measured_conus_profile=args.measured_conus_profile,segment=args.segment)
if args.station_cache:
    prior=Path(args.station_cache).resolve();meta=json.loads((prior/'stations.json').read_text())
    assert meta['binding']==binding,'A different profile/source cannot reuse this cache'
    prior_z=np.load(prior/'stations.npz');assert np.array_equal(prior_z['y'],ys)
    for y,target,record in zip(prior_z['y'],prior_z['targets'],meta['records']):cache[float(y)]=(target,record)
for refinement in range(8):
    if time.monotonic()-started>args.seconds:raise TimeoutError('Saved regional loft reached its native deadline')
    stage=out if refinement==0 else out/('profile-step-'+str(refinement))
    if refinement:stage.mkdir(exist_ok=False)
    data=[station(float(y))for y in ys]
    atomic_npz(stage/'stations.npz',y=ys,targets=np.array([p[0]for p in data]))
    atomic_json(stage/'stations.json',dict(binding=binding,records=[p[1]for p in data],
        complete_profile_stage=True,installable=False))
    try:
        reference_world=r['positions']*r['reference_radii']+r['reference_center']
        targets=np.array([p[0]for p in data])
        if radial_carrier:
            mapped=[radial_carrier_section(r,row[1]['inner_profile'],row[1]['outer_profile'],row[1]['center'],
                grey_faces,row[1]['desired_grey_fraction'],return_carrier=True)for row in data]
            fields=np.array([value[0]for value in mapped])
            section_proof=dict(all_triangles_positive=True,all_materials_share_one_map=True,
                topology_unchanged=True,common_piecewise_affine_carrier=True,
                minimum_double_area_m2=min(value[1]['minimum_double_area_m2']for value in mapped),
                regional_sections=[value[1]for value in mapped])
        elif args.successive_transport:
            fields,section_proof=transport_section_sequence(reference_world,r['faces'],fixed,targets,
                int(np.argmin(abs(ys-1.516826123))),interior_displacement_m=args.interior_section_repair_bound,
                wall_seconds=min(30.,args.seconds-(time.monotonic()-started)))
        else:fields,section_proof=harmonic_sections(reference_world,r['faces'],fixed,targets,weights='mean-value')
    except ValueError as error:
        if hasattr(error,'bad_section_faces'):
            atomic_npz(stage/'rejected-sections.npz',positions=error.section_fields,
                reference_positions=r['positions'],faces=r['faces'],bad_faces=error.bad_section_faces,
                fixed=fixed,y=ys,region=r['region'])
            atomic_json(stage/'rejected-sections.json',dict(passed=False,
                bad_faces=len(error.bad_section_faces),minimum_double_area_m2=error.minimum_section_double_area_m2,
                source_template_sha256=digest_file(template/'template.npz'),elapsed_seconds=time.monotonic()-started))
            if hasattr(error,'sequence_repairs'):atomic_json(stage/'rejected-kernel-repairs.json',error.sequence_repairs)
        if args.successive_transport and hasattr(error,'failed_sequence_station')and refinement<7:
            failed=error.failed_sequence_station;anchor=int(np.argmin(abs(ys-1.516826123)))
            previous=failed+1 if failed<anchor else failed-1
            gap=abs(ys[failed]-ys[previous])
            atomic_json(out/('rejected-step-'+str(refinement)+'.json'),dict(
                failed_y_m=float(ys[failed]),previous_y_m=float(ys[previous]),
                bad_faces=len(error.bad_section_faces),kernel_repairs=getattr(error,'sequence_repairs',[]),
                elapsed_seconds=time.monotonic()-started,installed=False))
            if gap<=.000020:raise ValueError('A shared profile step cannot be made positive within its declared 20-micrometre minimum interval')
            history.append(dict(round=refinement,stations=len(ys),failed_section_step=[float(ys[failed]),float(ys[previous])],
                maximum_interior_motion_m=args.interior_section_repair_bound))
            ys=np.sort(np.r_[ys,(ys[failed]+ys[previous])/2]);continue
        if not args.section_fit_limit or not hasattr(error,'section_fields'):raise
        affine=np.stack([r['positions']*row[1]['radii']+row[1]['center']for row in data])
        fields,fit=bound_section_deformations(affine,error.section_fields,r['faces'],maximum_displacement_m=args.section_fit_limit)
        atomic_json(stage/'section-fit.json',fit)
        if not fit['passed']:raise ValueError('Section fit exceeds its explicit source-target deviation bound: '+str(fit['maximum_target_deviation_m']))
        q=fields[:,r['faces']];aa=q[:,:,1]-q[:,:,0];bb=q[:,:,2]-q[:,:,0];area2=aa[:,:,0]*bb[:,:,1]-aa[:,:,1]*bb[:,:,0]
        section_proof=dict(all_triangles_positive=bool(np.all(area2>0)),minimum_double_area_m2=float(area2.min()),
            all_materials_share_one_map=True,topology_unchanged=True,bounded_fit=fit)
        assert section_proof['all_triangles_positive']
    positions=np.stack([np.c_[p[:,0],np.full(len(p),y),p[:,1]]for p,y in zip(fields,ys)])
    carrier_positions=(np.stack([np.c_[value[2][:,0],np.full(len(value[2]),y),value[2][:,1]]for value,y in zip(mapped,ys)])if radial_carrier else None)
    for end in end_sections.values():
        k=0 if end['y']==low else len(ys)-1
        positions[k,:,1]=end['plane']+positions[k,:,0]*end['slope'][0]+positions[k,:,2]*end['slope'][1]
        if radial_carrier:carrier_positions[k,:,1]=end['plane']+carrier_positions[k,:,0]*end['slope'][0]+carrier_positions[k,:,2]*end['slope'][1]
    if radial_carrier:
        gate_face=np.sort(r['carrier_faces'],axis=1);reference=r['carrier_positions'][gate_face]
        aa=reference[:,1]-reference[:,0];bb=reference[:,2]-reference[:,0]
        carrier_area=aa[:,0]*bb[:,1]-aa[:,1]*bb[:,0]
        gate_expected=np.r_[-np.sign(carrier_area),np.sign(carrier_area),-np.sign(carrier_area)]
        gate_positions=carrier_positions
    else:gate_face=face;gate_expected=expected;gate_positions=positions
    bad=[];minimum=np.inf;volume=0.
    for k in range(len(ys)-1):
        a,b,c=gate_positions[k,gate_face].transpose(1,0,2);dd,e,h=gate_positions[k+1,gate_face].transpose(1,0,2)
        determinant=np.r_[np.linalg.det(np.stack([b-a,c-a,h-a],axis=1)),
            np.linalg.det(np.stack([b-a,e-a,h-a],axis=1)),
            np.linalg.det(np.stack([dd-a,e-a,h-a],axis=1))]
        if np.any(np.sign(determinant)!=gate_expected):bad.append(k)
        minimum=min(minimum,float(abs(determinant).min()/6));volume+=float(abs(determinant).sum()/6)
    if not bad:
        print(json.dumps(dict(stage='positive-volume-carrier',stations=len(ys),seconds=time.monotonic()-started)),flush=True);break
    history.append(dict(round=refinement,stations=len(ys),bad_slabs=bad))
    if refinement==7:raise ValueError('Positive common prism interpolation was not established')
    ys=np.sort(np.r_[ys,[(ys[k]+ys[k+1])/2 for k in bad]])
else:raise AssertionError('No complete positive loft')
points=positions.reshape(-1,3);n=positions.shape[1];extra_points=[];extra_lookup={}
regions=m['regions'];features=m['features'];slab_labels=[];attributes={};label_cache={}
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
adjacency_edges=defaultdict(list)
for i,f in enumerate(r['faces']):
    for aa,bb in zip(f,np.roll(f,-1)):adjacency_edges[tuple(sorted((int(aa),int(bb))))].append(i)
adjacency=np.array([uses for uses in adjacency_edges.values()if len(uses)==2],int)
for a,b in zip(ys,ys[1:]):
    middle=(a+b)/2;segment=next((s['name']for s in segments if s['y_range_m'][0]<=middle<=s['y_range_m'][1]),None)
    assert segment, 'Every slab needs an explicit root-registered segment'
    if segment in label_cache:slab_labels.append(label_cache[segment]);continue
    index=levels.index(segment);names=[]
    for region in regions:
        active=next((features[i]for i in region['candidates']
            if any(levels.index(first)<=index<=levels.index(last)
                for first,last in features[i].get('level_ranges',[[features[i]['first'],features[i]['last']]]))),None)
        native=region['macro'];name=(active['name']if active else native)+'@'+segment
        row=dict(native_region=native,axial_region=segment,y_range_m=next(s['y_range_m']for s in segments if s['name']==segment),
            reference_estimate=bool(active),ta98_terms=([active['term'],*active['parents']]if active else []))
        if name in attributes:assert attributes[name]['native_region']==native, name
        attributes[name]=row;names.append(name)
    labels=np.asarray(names,dtype=object)[r['region']]
    matches=adjacency[labels[adjacency[:,0]]==labels[adjacency[:,1]]]
    graph=coo_matrix((np.ones(len(matches)*2),(np.r_[matches[:,0],matches[:,1]],np.r_[matches[:,1],matches[:,0]])),shape=(len(labels),len(labels)))
    _,component=connected_components(graph,directed=False)
    # Disconnected reference territories may meet at a single transverse
    # vertex. Keep their volumes in separately owned cells: extruding that
    # point as one material boundary would create a nonmanifold edge.
    for name in np.unique(labels):
        selected=np.flatnonzero(labels==name);pieces=np.unique(component[selected])
        if len(pieces)<2:continue
        info=attributes.pop(str(name))
        for number,piece in enumerate(pieces):
            key=str(name).split('@')[0]+'~island-'+str(number)+'@'+segment
            labels[selected[component[selected]==piece]]=key;attributes[key]=dict(info,reference_island=number)
    label_cache[segment]=labels;slab_labels.append(labels)
# One canonical owner pair per chart, with opposed sides on its two cells.
charts=defaultdict(list);cells=defaultdict(list)
def record(f,first,second):
    if first==second:return
    if first is None:first,second=second,None;f=f[:,[0,2,1]]
    if second is not None and first>second:first,second=second,first;f=f[:,[0,2,1]]
    charts[(first,second)].append(f)
edge_uses=defaultdict(list)
for i,f in enumerate(r['faces']):
    for a,b in zip(f,np.roll(f,-1)):edge_uses[tuple(sorted((int(a),int(b))))].append((i,int(a),int(b)))
assert all(len(u)<=2 for u in edge_uses.values())
if radial_carrier:
    edge_list=list(edge_uses);first_edge=np.array([edge_uses[e][0]for e in edge_list],int)
    edge_vertices=first_edge[:,1:];edge_first=first_edge[:,0]
    edge_second=np.array([edge_uses[e][1][0]if len(edge_uses[e])==2 else -1 for e in edge_list],int)
    edge_carrier=np.array([r['triangle_carrier'][i]for i in edge_first],int)
    carrier_sorted=np.sort(r['carrier_faces'],axis=1);cq=r['carrier_positions'][carrier_sorted[edge_carrier]]
    mat=np.stack([cq[:,1]-cq[:,0],cq[:,2]-cq[:,0]],axis=2)
    w=np.linalg.solve(mat[:,None],(r['positions'][edge_vertices]-cq[:,None,0])[...,None])[...,0]
    edge_bary=np.concatenate([(1-w.sum(2))[...,None],w],axis=2)
    extra_protected=[]
    physical=set(map(int,r['fixed_controls']))
    def mapped_edges(k,active):
        vertices=edge_vertices[active];owners=carrier_sorted[edge_carrier[active]]
        cuts,mapped=carrier_prism_edge_grid(edge_bary[active],carrier_positions[k,owners],carrier_positions[k+1,owners])
        ids=np.empty(cuts.shape,int);ids[:,:,0]=vertices+k*n;ids[:,:,3]=vertices+(k+1)*n
        middle=cuts[:,:,1:3];candidate=mapped[:,:,1:3];node=np.broadcast_to(vertices[:,:,None],middle.shape)
        flat_t=middle.ravel();flat_v=node.ravel();flat_p=candidate.reshape(-1,3)
        zero=flat_t<=1e-12;one=flat_t>=1-1e-12;inside=~(zero|one);mid_ids=np.empty(len(flat_t),int)
        mid_ids[zero]=flat_v[zero]+k*n;mid_ids[one]=flat_v[one]+(k+1)*n
        key=np.rec.fromarrays([flat_v[inside],np.round(flat_t[inside],12)],names='vertex,t')
        unique,index,inv=np.unique(key,return_index=True,return_inverse=True)
        offset=len(points)+len(extra_points)
        mid_ids[inside]=offset+inv
        for p in flat_p[inside][index]:extra_points.append(p)
        extra_protected.extend(-1 if int(vertex)in physical else 0 for vertex in unique.vertex)
        ids[:,:,1:3]=mid_ids.reshape(middle.shape)
        pieces=[]
        for j in range(3):pieces.extend([ids[:,[0,0,1],[j,j+1,j+1]],ids[:,[0,1,1],[j,j+1,j]]])
        triangles=np.stack(pieces,axis=1)
        valid=np.all(np.diff(np.sort(triangles,axis=2),axis=2)!=0,axis=2)
        return triangles,valid
for k,labels in enumerate(slab_labels):
    if time.monotonic()-started>args.seconds:raise TimeoutError('Regional shared-chart emission reached its native deadline')
    if radial_carrier:
        second=np.array([None if j<0 else str(labels[j])for j in edge_second],object)
        first=labels[edge_first];active=np.flatnonzero(first!=second)
        triangles,valid=mapped_edges(k,active);groups=defaultdict(list)
        for local,i in enumerate(active):groups[(str(first[i]),second[i])].append(local)
        for(owner,peer),ids0 in groups.items():
            record(triangles[ids0][valid[ids0]],owner,peer)
    else:
        for uses in edge_uses.values():
            i,a,b=uses[0];first=str(labels[i]);second=str(labels[uses[1][0]])if len(uses)==2 else None
            if first==second:continue
            aa,bb=a+k*n,b+k*n;dd,cc=aa+n,bb+n
            quad=(np.array([[aa,dd,cc],[aa,cc,bb]])if a<b else np.array([[aa,dd,bb],[dd,cc,bb]]))
            record(quad,first,second)
    lower=None if k==0 else slab_labels[k-1]
    groups=defaultdict(list)
    for i,name in enumerate(labels):
        below=None if lower is None else str(lower[i]);name=str(name)
        if below!=name:groups[(name,below)].append(i)
    for(first,second),ids0 in groups.items():record(r['faces'][ids0]+k*n,first,second)
    if k%32==0:print(json.dumps(dict(stage='shared-face-batch',slab=k,slabs=len(slab_labels),seconds=time.monotonic()-started)),flush=True)
last=slab_labels[-1];groups=defaultdict(list)
for i,name in enumerate(last):groups[str(name)].append(i)
for name,ids0 in groups.items():record((r['faces'][ids0]+(len(ys)-1)*n)[:,[0,2,1]],name,None)
if extra_points:points=np.vstack([points,np.array(extra_points)])
surfaces={}
for number,((first,second),pieces)in enumerate(sorted(charts.items(),key=lambda row:str(row[0]))):
    sid='regional-interface:'+str(number);surfaces[sid]=np.concatenate(pieces)
    cells[first].append(dict(surface=sid,side=1,neighbor=second))
    if second is not None:cells[second].append(dict(surface=sid,side=-1,neighbor=first))
checks={}
for name,uses in cells.items():
    if args.emit_stage:break
    if time.monotonic()-started>args.seconds:raise TimeoutError('Regional topology checks reached their native deadline')
    f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
    used,inv=np.unique(f,return_inverse=True);check=topology(Mesh(points[used],inv.reshape(-1,3),'female-shared-regional-cord'))
    checks[name]={key:check[key]for key in ['closed_oriented_manifold','open_edges','nonmanifold_edges','zero_area_faces','components','oriented_triangle_integral_m3']}
bad=[name for name,c in checks.items()if not c['closed_oriented_manifold']]
atomic_npz(out/'cells.npz',positions=points,**surfaces)
interior=np.zeros(n,int);interior[r['fixed_controls']]=-1
protected=np.r_[np.tile(interior,len(ys)),np.array(extra_protected,int)if radial_carrier else np.zeros(len(extra_points),int)]
atomic_npz(out/'domain.npz',positions=points,interior_material=protected)
atomic_json(out/'domain.json',dict(source_template_sha256=digest_file(template/'template.npz'),
    source_cells_sha256=digest_file(out/'cells.npz'),
    physical_boundary_controls_protected=True,calculation_volume_not_stored_twice=True))
atomic_json(out/'cells.json',dict(cells=dict(cells),surfaces={s:len(f)for s,f in surfaces.items()},
    attributes=attributes,topology=checks,passed=not bad and not args.emit_stage,
    completed_construction_stage=True,source_topology_pending=bool(args.emit_stage),segment=args.segment or 'ALL',
    bands=[[s['name'],max(low,s['y_range_m'][0]),min(high,s['y_range_m'][1])]for s in segments if s['y_range_m'][0]<high and s['y_range_m'][1]>low],
    stations=[row[1]for row in data],
    construction=dict(shared_cross_section_graph=True,source_template_sha256=digest_file(template/'template.npz'),
        source_fitted_tissues_sha256=digest_file(source/'cells.npz'),
        reciprocal_axial_and_radial_interfaces=True,central_canal_is_void=True,
        embedded_charts_follow_volume_carrier=radial_carrier,
        microscopic_borders_are_reference_estimates=True,all_positive_reference_jacobians=True,
        minimum_tetrahedron_volume_m3=minimum,source_volume_m3=volume,section_proof=section_proof,
        adaptive_axial_refinement=history,
        terminal_sections={key:{'y':e['y'],'maximum_source_cap_offset_m':e['maximum_source_cap_offset_m']}for key,e in end_sections.items()},
        complete_emitted_geometry_source_deviation_fixed_spine_and_anatomical_checks_pending=True),
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
assert not bad,bad[:10]
print(json.dumps(dict(passed=not args.emit_stage,construction_complete=True,topology_pending=args.emit_stage,cells=len(cells),stations=len(ys),positions=len(points),
    triangles=sum(len(f)for f in surfaces.values()),seconds=time.monotonic()-started)),flush=True)
