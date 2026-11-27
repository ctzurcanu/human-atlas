"""One bounded original partial floor-of-mouth mucosa/fascia construction.

The two adult primary studies establish a continuous submucosal fascial roof,
not a cord; they do not provide numerical layer thicknesses. Estimates here
are explicitly illustrative. This script never searches shape parameters,
modifies source anatomy or installs an incomplete cell.
"""
from pathlib import Path
import argparse,copy,json,sys,time,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
import triangle
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from shapely.geometry import Polygon,Point
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.tissue import pack,normal_words
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_continuity import boundary_loops,annular_rim,finite_wall_cells
from tongue_whole_input import WORK,fixed_bone_sources

SOURCE=WORK/'shared-coverings-exposed-apo-sl-v3/packed/tissue.json'
ATTACH=WORK/'floor-of-mouth-attachment-v2/attachment-witnesses.json'
CONTRACT=WORK/'floor-of-mouth-contract-v1/contract.json'
REFS=['https://onlinelibrary.wiley.com/doi/10.1002/ca.23343',
      'https://onlinelibrary.wiley.com/doi/full/10.1155/2020/1820978']

def bound(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=digest_file(p))

def require(x,s):
    if not x:raise ValueError(s)

def arc_parameter(p):
    d=np.r_[0.,np.linalg.norm(np.diff(p,axis=0),axis=1).cumsum()]
    require(np.all(np.diff(d)>0),'Repeated arc controls');return d/d[-1]

def plan(output):
    require(not output.exists()and output.is_relative_to(WORK),'Use a new local stage')
    a=json.loads(ATTACH.read_text());c=json.loads(CONTRACT.read_text());s=json.loads(SOURCE.read_text())
    require(s['passed']and digest_file(s['bank'])==s['bank_sha256'],'Source cell bank changed/failed')
    rim=np.load(c['tongue_reflection']['source']['path'],allow_pickle=False)
    # The actual left/right rim vertices just anterior to the gland-roof
    # witnesses bound this partial construction. They are computational review
    # cuts, not claimed anatomical posterior borders.
    ids=list(range(56,64))+list(range(0,9))
    output.mkdir(parents=True)
    result=dict(schema='tongue-partial-fom-plan/v1',source=bound(SOURCE),attachment=bound(ATTACH),contract=bound(CONTRACT),
        inventory=bound(ROOT/'reports/ta98-tongue-inventory.json'),primary_references=REFS,
        primary_constraints=dict(continuous_mucosa_and_submucosal_fascia=True,
            attachment='Broad inner mandibular periosteal arc to anterior/lateral ventral tongue; fascia lies above fixed glands/ducts/GG.',
            thickness='Adults vary. The anterior fascia is thicker and gradually becomes thinner posterolaterally; published studies do not give measured layer thicknesses.',
            posterior='Thin fascia continues posteriorly under the lateral tongue into the piriform region; therefore these finite lateral cuts are preview extents, not anatomical terminations.',
            frenulum='A central continuous fold of this mucosa/fascia field. Current neutral pose retained; no discrete AP cord or forced tongue movement.'),
        explicit_original_estimates=dict(mucosa_m=.00012,fascia_anterior_m=.00010,fascia_lateral_cut_m=.00006,
            local_duct_lip_radial_m=.00008,midline_frenulum_relief_m=.00025,
            status='Prototype values, not donor measurements or a reference-supported quantitative thickness claim'),
        tongue_seam=dict(source=c['tongue_reflection']['source'],loop_indices=ids,
            source_position_indices=rim['source_position_indices'][ids].tolist(),
            posterior_endpoints_m=rim['positions'][ids][[0,-1]].tolist(),
            rule='Reuse the literal source outer/inner mucosal rim facet pairs with reciprocal cell ownership; no duplicate cap or coordinate-only seam claim.'),
        field_extent='Anterior and bilateral lateral horseshoe field between the demonstrated inner gingival arc and this ventral tongue seam; posterior two side cuts remain explicit partial-review boundaries.',
        openings=[dict(part_id=r['part_id'],loop=r['source_loop'],
            representation='Open twelve-sided finite mucosal/fascial collar around the actual source oral endpoint; lumen remains empty',
            native_full_duct_wall_verified=False)for r in a['duct_endpoint_identity']if r['accepted_endpoint_identity']=='oral'],
        existing_gingiva_clashes=a['source_gland_gum_crossings'],
        fixed='All existing glands, nerves, muscles and skeleton literal. No whole-gingiva repair.',
        adjoining_interfaces=['Current tongue mucosa / FOM mucosa: literal reciprocal source rim',
            'FOM mucosa / FOM fascia: one common chart, opposite cell ownership',
            'Gingival arc: witnessed adjacency only; no claim of repaired native gingiva',
            'Duct mouth: original twelve-point loop retained, open local collar; full native duct wall remains unbuilt'],
        installed=False,anatomical_acceptance=False)
    atomic_json(output/'plan.json',result);print(json.dumps({'plan':str(output/'plan.json')}),flush=True)

def old_words(source,surfaces):
    bank=Path(source['bank']).read_bytes();words={}
    for cell,row in source['parts'].items():
        n=row['vertexCount'];f=np.ndarray((row['indexCount']//3,3),'<u4',buffer=bank,offset=row['indices'])
        nw=np.ndarray((n,3),'<i2',buffer=bank,offset=row['normals']);ns=np.ndarray((n,),'i1',buffer=bank,offset=row['normalSigns'])
        off=0
        for use in source['cells'][cell]:
            sid=use['surface'];amount=len(surfaces[sid]);w=nw[f[off:off+amount]].astype(np.int32)*ns[f[off:off+amount]][:,:,None]
            if use['side']==-1:w=-w[:,[0,2,1]]
            if sid in words:require(np.array_equal(words[sid],w),'Original reciprocal words differ')
            else:words[sid]=w.astype('<i2')
            off+=amount
    return words

def construct(plan_path):
    began=time.monotonic();output=plan_path.parent/'candidate';require(not output.exists(),'No retry or overwrite of this stage')
    p=json.loads(plan_path.read_text());require(digest_file(SOURCE)==p['source']['sha256']and digest_file(ATTACH)==p['attachment']['sha256'],'Plan binding changed')
    source=json.loads(SOURCE.read_text());a=json.loads(ATTACH.read_text());c=json.loads(CONTRACT.read_text())
    with np.load(source['arrangement'],allow_pickle=False)as raw:
        old=raw['positions'].copy();surfaces={k:raw[n].copy()for k,n in source['surfaces'].items()}
    words=old_words(source,surfaces);cells=copy.deepcopy(source['cells']);part_ids=dict(source['part_ids'])
    rr=np.load(c['tongue_reflection']['source']['path'],allow_pickle=False)
    inner_ids=rr['source_position_indices'][p['tongue_seam']['loop_indices']]
    inner=old[inner_ids];it=arc_parameter(inner);outer=np.array([r['point_m']for r in a['witnesses']]);ot=arc_parameter(outer)
    # Both source arcs start on the same positive-X side.
    require(inner[0,0]>0 and outer[0,0]>0,'Source arc ordering changed')
    oc=CubicSpline(ot,outer,axis=0,bc_type='natural');ic=CubicSpline(it,inner,axis=0,bc_type='natural')
    outer_set=set(np.unique(surfaces['mucosa:outer']));paired={}
    free=surfaces['mucosa:free'];join_edge_set={tuple(sorted((int(x),int(y))))for x,y in zip(inner_ids[:-1],inner_ids[1:])}
    seam_faces=[]
    for i,f in enumerate(free):
        flags=np.array([x in outer_set for x in f])
        if flags.sum()==2:
            j=int(np.flatnonzero(~flags)[0]);aa,bb,cc=int(f[(j+1)%3]),int(f[(j+2)%3]),int(f[j]);paired[bb]=cc
            if tuple(sorted((aa,bb)))in join_edge_set:
                seam_faces.append(i)
                # The source split diagonal is retained literally.
                match=[k for k,q in enumerate(free)if k!=i and set([aa,cc]).issubset(q)and sum(x in outer_set for x in q)==1]
                require(len(match)==1,'Source rim pair is ambiguous');seam_faces+=match
    seam_faces=np.unique(seam_faces);require(len(seam_faces)==2*(len(inner_ids)-1),'Source seam strip incomplete')
    ii=np.array([paired[int(q)]for q in inner_ids]);inner_depth=old[ii]-inner
    # Source outer normals averaged from the actual exterior faces.
    sn=np.zeros_like(old);tri=old[surfaces['mucosa:outer']];nn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    for k in range(3):np.add.at(sn,surfaces['mucosa:outer'][:,k],nn)
    sn/=np.maximum(np.linalg.norm(sn,axis=1),1e-30)[:,None];nc=CubicSpline(it,sn[inner_ids],axis=0,bc_type='natural')
    mouths=[]
    for row in p['openings']:
        z=np.load(row['loop']['path'],allow_pickle=False);q=z['positions'].copy();center=q.mean(0);normal=np.linalg.svd(q-center)[2][-1]
        if normal[1]<0:normal=-normal
        mouths.append(dict(row=row,points=q,center=center,normal=normal))
    vm=.45
    def midbase(u):return (oc(u)+ic(u))*.5
    for m in mouths:
        target=m['center'][0]
        m['u']=brentq(lambda u:float(midbase(u)[0]-target),.15,.85)
    midu=np.array([0.,mouths[0]['u'],.5,mouths[1]['u'],1.]);order=np.argsort(midu);midu=midu[order]
    middle=(mouths[0]['center']+mouths[1]['center'])/2;middle[0]=0
    midpoints=np.array([midbase(0),mouths[0]['center'],middle, mouths[1]['center'],midbase(1)])[order]
    require(np.all(np.diff(midu)>0),'Oral endpoint parameter ordering invalid')
    mc=CubicSpline(midu,midpoints,axis=0,bc_type='natural')
    nmid=CubicSpline(np.array([0.,mouths[0]['u'],mouths[1]['u'],1.])[np.argsort([0.,mouths[0]['u'],mouths[1]['u'],1.])],
        np.array([[0,1,0],mouths[0]['normal'],mouths[1]['normal'],[0,1,0]])[np.argsort([0.,mouths[0]['u'],mouths[1]['u'],1.])],axis=0,bc_type='natural')
    def h(x):return x*x*(3-2*x)
    def hermite(a,b,da,db,t,length):
        t=np.asarray(t);return ((2*t**3-3*t*t+1)[:,None]*a+(-2*t**3+3*t*t)[:,None]*b+
            (t**3-2*t*t+t)[:,None]*da*length+(t**3-t*t)[:,None]*db*length)
    def field(uv):
        u,v=uv.T;o=oc(u);i=ic(u);m=mc(u);du=ic(u,1);N=nc(u);N/=np.linalg.norm(N,axis=1)[:,None]
        di=np.cross(du,N);di/=np.linalg.norm(di,axis=1)[:,None];di*=.008
        nm=nmid(u);nm/=np.linalg.norm(nm,axis=1)[:,None]
        dm=np.cross(du,nm);dm/=np.linalg.norm(dm,axis=1)[:,None];dm*=.010
        do=(m-o)/vm;out=np.empty_like(o);take=v<=vm
        out[take]=hermite(o[take],m[take],do[take],dm[take],v[take]/vm,vm)
        out[~take]=hermite(m[~take],i[~take],dm[~take],di[~take],(v[~take]-vm)/(1-vm),1-vm)
        # Small continuous central reflection relief; no extra solid cord.
        ridge=p['explicit_original_estimates']['midline_frenulum_relief_m']*np.exp(-((u-.5)/.07)**2)*np.sin(np.pi*v)**4
        out[:,1]+=ridge
        return out
    # Literal source arc nodes define the boundary; constrained triangulation
    # retains them and both native 12-point oral loops.
    uv=[];segments=[];roles=[]
    def polygon(points,role):
        start=len(uv);uv.extend(points);roles.extend([role]*len(points))
        segments.extend([(start+j,start+(j+1)%len(points))for j in range(len(points))]);return np.arange(start,start+len(points))
    rect=np.r_[np.c_[ot,np.zeros(len(ot))],[[1.,v]for v in [.2,.4,.6,.8]],
        np.c_[it[::-1],np.ones(len(it))],[[0.,v]for v in [.8,.6,.4,.2]]]
    rectangle_ids=polygon(rect,'boundary');hole_polygons=[]
    for k,m in enumerate(mouths):
        at=np.array([[m['u'],vm]]);eps=1e-5
        J=np.c_[((field(at+[eps,0])-field(at-[eps,0]))/(2*eps))[0],((field(at+[0,eps])-field(at-[0,eps]))/(2*eps))[0]]
        local=(m['points']-m['center'])@np.linalg.pinv(J).T
        hp=local+at[0];require(np.all((hp>0)&(hp<1)),'Native mouth cannot fit this bounded field')
        m['uv']=hp;m['parameter_ids']=polygon(hp,'mouth:'+str(k));hole_polygons.append(Polygon(hp))
    for u in np.linspace(0,1,33)[1:-1]:
        for v in np.linspace(0,1,15)[1:-1]:
            q=np.array([u,v])
            if any(poly.distance(Point(q))<.025 for poly in hole_polygons):continue
            uv.append(q);roles.append('interior')
    uv=np.array(uv);holes=np.array([[m['u'],vm]for m in mouths]);domain=triangle.triangulate(dict(vertices=uv,segments=np.array(segments),holes=holes),'pQ')
    require(len(domain['vertices'])==len(uv)and np.array_equal(domain['vertices'],uv),'No source boundary Steiner controls permitted')
    f=domain['triangles'][:,[0,2,1]];top=field(uv)
    top[np.isclose(uv[:,1],0)]=oc(uv[np.isclose(uv[:,1],0),0])
    seam_local=np.array([int(np.flatnonzero((uv[:,0]==u)&(uv[:,1]==1))[0])for u in it])
    top[seam_local]=inner
    for m in mouths:top[m['parameter_ids']]=m['points']
    eps=1e-5;Bu=(field(uv+[eps,0])-field(uv-[eps,0]))/(2*eps);Bv=(field(uv+[0,eps])-field(uv-[0,eps]))/(2*eps)
    normals=-np.cross(Bu,Bv);normals/=np.maximum(np.linalg.norm(normals,axis=1),1e-30)[:,None]
    estimates=p['explicit_original_estimates'];mu=estimates['mucosa_m']
    fascia=estimates['fascia_lateral_cut_m']+(estimates['fascia_anterior_m']-estimates['fascia_lateral_cut_m'])*np.sin(np.pi*uv[:,0])**2
    first=top-mu*normals;second=first-fascia[:,None]*normals
    first[seam_local]=old[ii];second[seam_local]=first[seam_local]-fascia[seam_local,None]*nc(it)
    for m in mouths:
        ids=m['parameter_ids'];rel=m['points']-m['center'];rad=np.linalg.norm(rel,axis=1)
        inner_ring=m['center']+rel*(1-estimates['local_duct_lip_radial_m']/rad)[:,None]
        first[ids]=inner_ring-m['normal']*mu
        second[ids]=inner_ring-m['normal']*(mu+fascia[ids])[:,None]
    fields=np.array([top,first,second]);origin=np.array(json.loads((SOURCE.parent.parent/'construction.json').read_text())['storageOrigin'])
    fields=(fields-origin).astype('<f4').astype(float)+origin
    # Shared source coordinates must survive the exact emitted bank encoding.
    require(np.array_equal(fields[0,seam_local],inner)and np.array_equal(fields[1,seam_local],old[ii]),'Actual stored seam would move')
    n=len(uv);positions=np.r_[old,fields.reshape(-1,3)];mapping=np.arange(len(positions));base=len(old)
    mapping[base+seam_local]=inner_ids;mapping[base+n+seam_local]=ii
    fresh={};f0=mapping[base+f];f1=mapping[base+n+f];f2=mapping[base+2*n+f]
    fresh['FOM:oral']=f0;fresh['FOM:mucosa:fascia']=f1;fresh['FOM:deep']=f2
    loops=boundary_loops(f);rim_fascia=[];rim_mu=[]
    # Keep every true opening a hole. Annular rims span wall thickness only.
    for loop in loops:
        for aa,bb in zip(loop,np.roll(loop,-1)):
            if not(uv[aa,1]==1 and uv[bb,1]==1):
                rim_mu.extend([[bb,aa,aa+n],[bb,aa+n,bb+n]])
            rim_fascia.extend([[bb+n,aa+n,aa+2*n],[bb+n,aa+2*n,bb+2*n]])
    fresh['FOM:mucosa:rims']=mapping[base+np.array(rim_mu)];fresh['FOM:fascia:rims']=mapping[base+np.array(rim_fascia)]
    keep=np.setdiff1d(np.arange(len(free)),seam_faces);surfaces['mucosa:free']=free[keep]
    surfaces['tongue:FOM:literal-rim']=free[seam_faces];words['tongue:FOM:literal-rim']=words['mucosa:free'][seam_faces];words['mucosa:free']=words['mucosa:free'][keep]
    surfaces.update(fresh)
    cells['mucosa'].append(dict(surface='tongue:FOM:literal-rim',side=1,neighbor='FOM-mucosa'))
    cells['FOM-mucosa']=[dict(surface='FOM:oral',side=1,neighbor=None),dict(surface='FOM:mucosa:fascia',side=-1,neighbor='FOM-fascia'),dict(surface='FOM:mucosa:rims',side=1,neighbor=None),dict(surface='tongue:FOM:literal-rim',side=-1,neighbor='mucosa')]
    cells['FOM-fascia']=[dict(surface='FOM:mucosa:fascia',side=1,neighbor='FOM-mucosa'),dict(surface='FOM:deep',side=-1,neighbor=None),dict(surface='FOM:fascia:rims',side=1,neighbor=None)]
    part_ids.update({'FOM-mucosa':'TA98TONGUE:male:prototype:FOM-mucosa','FOM-fascia':'TA98TONGUE:male:prototype:FOM-fascia'})
    output.mkdir();atomic_npz(output/'prepack-arrangement.npz',positions=positions,**{'surface_'+str(i):x for i,x in enumerate(surfaces.values())},fields=fields,carrier_faces=f,parameters=uv)
    newmeshes={}
    for name in ['FOM-mucosa','FOM-fascia']:
        faces=np.r_[tuple(surfaces[x['surface']]if x['side']==1 else surfaces[x['surface']][:,[0,2,1]]for x in cells[name])]
        mesh=Mesh(positions,faces,source['frame']);mesh.save(output/(name+'.npz'));newmeshes[name]=mesh
    gate=finite_wall_cells(fields,f)
    atomic_json(output/'prepack.json',dict(plan=bound(plan_path),surfaces={k:'surface_'+str(i)for i,k in enumerate(surfaces)},cells=cells,part_ids=part_ids,frame=source['frame'],origin=origin.tolist(),finite_prism_checks=gate,
        seam_source_faces=seam_faces.tolist(),shared_source_rim_coordinates_literal=True,
        mouth_channels=[dict(part_id=m['row']['part_id'],source_loop=m['row']['loop'],source_outer_points_literal=np.array_equal(fields[0,m['parameter_ids']],m['points']),minimum_void_radius_m=float(np.linalg.norm(first[m['parameter_ids']]-np.mean(first[m['parameter_ids']],axis=0),axis=1).min()))for m in mouths]))
    # The exact original raw charts are retained; new charts get their actual
    # stored geometric normals. A failed prism/packing gate remains failed.
    words.update(normal_words((positions-origin).astype('<f4').astype(float),fresh))
    result=pack(positions,surfaces,cells,part_ids,source['frame'],output/'packed',origin=origin,corner_words=words,
        construction=dict(kind='bounded-source-bound-partial-floor-mouth',finite_wall_cells_pass=gate['finite_wall_cells_pass'],plan_sha256=digest_file(plan_path),anatomical_acceptance=False),
        max_pairs=100000,max_tests=1000000,checkpoint_directory=output/'checks')
    fixed=[]
    current=json.loads(CONTRACT.read_text())
    for row in current['source_exports']:
        if row['part_id']=='ZA:Tongue':continue
        fixed.append((row['part_id'],Path(row['mesh']['path'])))
    selected=WORK/'hyoglossus-origin-v7-nonplanar-rim-topology/selected-extrinsics.json'
    ex=json.loads(selected.read_text())
    fixed.extend((r['originalPartId'],Path(r['mesh']))for r in ex['parts'])
    fixed.append(('ZA:Gingiva',Path(json.loads((WORK/'floor-of-mouth-guides-v1/guides.json').read_text())['gingiva']['mesh']['path'])))
    fixed.extend((r['partId'],Path(r['meshSavePath']))for r in fixed_bone_sources()if r['partId']not in {q[0]for q in fixed})
    context=[]
    for name,m in newmeshes.items():
        for pid,path in fixed:
            _,s=pair_contacts(m,Mesh.load(path),0.,max_pairs=30000,max_tests=300000)
            context.append(dict(first=name,second=pid,source=bound(path),summary=s))
    clear=all(r['summary']['complete_within_pair_budget']and not(r['summary']['crossing_face_pairs']or r['summary']['coplanar_overlap_face_pairs'])for r in context)
    atomic_json(output/'construction.json',dict(schema='tongue-partial-fom-candidate/v1',plan=bound(plan_path),tissue_manifest=bound(output/'packed/tissue.json'),
        passed=bool(result['passed']and clear),new_cell_checks=result['checks'],finite_prism_checks=gate,fixed_context=context,
        existing_gingiva_clashes_separate=p['existing_gingiva_clashes'],shared_rim_source_faces=len(seam_faces),
        sources_modified=False,installed=False,anatomical_acceptance=False,complete_FOM=False,
        no_parameter_search=True,elapsed_seconds=time.monotonic()-began))
    atomic_json(output/'review-spec.json',dict(frame=source['frame'],resolution=768,wall_seconds=45,
        meshes=[dict(mesh=str(output/'FOM-mucosa.npz'),name='Partial finite FOM mucosa - uninstalled prototype',color=[.78,.42,.45,1]),dict(mesh=str(output/'FOM-fascia.npz'),name='Partial finite FOM fascia - estimated thickness',color=[.76,.70,.57,1])],
        views=[dict(name='oral-superior',direction=[0,1,0]),dict(name='oral-oblique',direction=[0,.7,1]),dict(name='sagittal-section',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.622,.066],clear_outer=True))],
        qualification='Bounded original uninstalled partial field. Primary sources give no numeric thickness. Existing fixed gingiva clashes reported separately; never install a failed candidate.'))
    print(json.dumps(dict(output=str(output),passed=result['passed']and clear,bank_passed=result['passed'],fixed_context_clear=clear,finite_prism_passed=gate['finite_wall_cells_pass'],elapsed_seconds=time.monotonic()-began)),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();sp=ap.add_subparsers(dest='action',required=True);q=sp.add_parser('plan');q.add_argument('--output',type=Path,required=True);q=sp.add_parser('construct');q.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    if args.action=='plan':plan(args.output.resolve())
    else:
        out=args.plan.resolve().parent/'candidate'
        try:construct(args.plan.resolve())
        except Exception as e:
            out.mkdir(exist_ok=True);atomic_json(out/'construction-failure.json',dict(passed=False,installed=False,no_parameter_search=True,error=type(e).__name__+': '+str(e),traceback=traceback.format_exc()));raise
