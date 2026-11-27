"""Canonical renal construction retaining the registered medullary topology.

This module creates candidates only. Geometry, encoded interfaces, sections,
recipient packing and source-bound local installation are separate gates.
"""
import time,json
import numpy as np
import manifold3d as md
from scipy.spatial import ConvexHull
from .renal_master import bound,OUT,FRAME
from .boolean_cells import BooleanCells
from .common import atomic_json,atomic_npz,digest_file,ROOT
from .model import Mesh
from .audit import audit

ORIGIN=np.array([.0625,1.1875,-.015625])


def hull_mesh(points):
    v=np.unique(points,axis=0);h=ConvexHull(v);f=h.simplices.copy();q=v[f];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);wrong=np.einsum('ij,ij->i',n,h.equations[:,:3])<0;f[wrong]=f[wrong][:,[0,2,1]]
    return v,f


def portal_mesh(points):
    cy,cz=(points[:,1:].min(0)+points[:,1:].max(0))/2;ry,rz=np.ptp(points[:,1:],axis=0)*.58*1.15
    a=np.linspace(0,2*np.pi,64,endpoint=False)
    v=np.array([[x,cy+ry*np.cos(t),cz+rz*np.sin(t)]for x in [-.02,points[:,0].max()+.0012]for t in a]);f=[]
    for i in range(64):j=(i+1)%64;f.extend([[i,j,j+64],[i,j+64,i+64]])
    for i in range(1,63):f.extend([[0,i+1,i],[64,64+i,64+i+1]])
    return v,np.array(f,int)


def build(version='441',*,fine_regions=True):
    out=OUT/('native-v'+version)
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);start=time.monotonic();b=BooleanCells();meta={};mr,mz=bound('medullary-collecting-cells-v300');fr,fz=bound('renal-ureter-fluid-v317')
    plan=dict(schema='renal-native-common-cells/v1',inventory_sha256=digest_file(ROOT/'TA98-KIDNEY-PARTS.md'),master='male-left',
              source_reports={n:digest_file(ROOT/'reports'/(n+'.json'))for n in ['medullary-collecting-cells-v300','renal-ureter-fluid-v317']},
              medullary_topology='Retain the actual registered eleven bulk pyramids and eleven papillae; no signed-field resampling.',
              cortex='Native inner capsule, excluding actual medullary/collecting occupancy and a declared central sinus compartment; one peripheral and inter-pyramid cortical arrangement.',
              sinus='Convex calyceal landmark envelope restricted to the 6.6 mm eroded core, continuing to the existing rounded hilum. Fat is a separate content, never a plug across urine.',
              region_limits='Gross zonal and cortical-depth boundaries are illustrative, based on histological organization rather than patient-specific measurements. Ray bundles are representative.',
              collecting_coats_m=[.00006,.00012],peripheral_cortical_depth_m=.0066,fixed_skeleton=True,installed=False,anatomical_acceptance=False)
    atomic_json(out/'plan.json',plan)
    def source(v,f,name):return b.source(v-ORIGIN,f,name)
    def derive(m,name,**kwargs):return b.derived_source(m,name,**kwargs)
    def stage(text):print(text,round(time.monotonic()-start,2),flush=True)
    K=source(mz['native_kidney_v'],mz['native_kidney_f'],'kidney-outer');I=source(mz['native_inner_v'],mz['native_inner_f'],'capsule-inner')
    z=np.load(OUT/'master-v428/source-fields.npz');pv,pf=portal_mesh(z['portal_points']);P=source(pv,pf,'rounded-hilum')
    O=source(fz['positions'],fz['outerClassifier'],'collecting-outer');F=source(fz['positions'],fz['innerClassifier'],'collecting-inner')
    W=source(fz['positions'],fz['wholeWallFaces'],'actual-complete-collecting-wall')
    sf={s['id']:mz[s['facesArray']]for s in mr['surfaces']};medulla={};exterior=[]
    for name,uses in mr['cells'].items():
        if name=='collecting-wall':continue
        f=np.concatenate([sf[u['surface']]if u['side']==1 else sf[u['surface']][:,[0,2,1]]for u in uses]);m=source(mz['positions'],f,name);medulla[name]=m
        for u in uses:
            if u.get('neighbor')in mr['cells']and u['neighbor']!='collecting-wall':continue
            exterior.append(sf[u['surface']]if u['side']==1 else sf[u['surface']][:,[0,2,1]])
    allM=source(mz['positions'],np.concatenate(exterior),'actual-medullary-exterior')
    Q=derive(allM+O,'occupied-medullary-and-collecting-envelope',retain_parents=True)
    stage('Original medullary cells bound')
    # Rebind derived classifier fields once, before any tissue partition, so
    # every owner uses the very same original planes and provenance IDs.
    ball=md.Manifold.sphere(1.,12)
    E=derive(I.minkowski_difference(ball.scale((.0066,)*3)),'cortical-base')
    internal_O=O^I;om=internal_O.to_mesh64();ov=np.asarray(om.vert_properties)[:,:3]
    bm=ball.to_mesh64();margin=np.asarray(bm.vert_properties)[:,:3]*.0015
    hv,hf=hull_mesh(np.concatenate([ov+d for d in margin]));S0=b.source(hv,hf,'calyceal-sinus-envelope');S=(S0^E)+P
    cortex=I-(Q+S);sinus=(S^I)-Q;capsule=(K-I)-P
    stage('Complete native cortex and declared sinus constructed')
    coats=[]
    for depth in [.00006,.00012]:coats.append(derive(O.minkowski_difference(ball.scale((depth,)*3)),'collecting-depth-'+str(depth)))
    A,B=coats
    solids={'capsule':capsule,'columns':(I^E)-(Q+S0+P),'cortex':(I-E)-(Q+P),'sinus-fat':sinus,
            'collecting-adventitia':W-A,'collecting-muscularis':(W^A)-B,'collecting-mucosa':W^B,**medulla}
    stage('Three nested collecting coats constructed')
    if fine_regions:
        # Planes divide only actual pyramid tissue. The canonical face overlay
        # also subdivides adjoining cortex/contact charts without duplicating
        # tissue or changing any source papilla.
        for letter in 'abcdefghijk':
            name='pyramid:'+letter;uses=mr['cells'][name];base=next(u for u in uses if u.get('neighbor')=='papilla:'+letter)
            tip=mz['positions'][np.unique(sf[base['surface']])].mean(0)-ORIGIN
            v=mz['positions'][np.unique(np.concatenate([sf[u['surface']]for u in uses]))]-ORIGIN
            direction=v.mean(0)-tip;direction/=np.linalg.norm(direction);length=float(((v-tip)@direction).max())
            # Build generous clipped halfspaces from actual local landmarks.
            def half(fraction):
                center=tip+direction*length*fraction;u=np.cross(direction,np.eye(3)[np.argmin(abs(direction))]);u/=np.linalg.norm(u);w=np.cross(direction,u)
                vertices=np.array([center+u*x+w*y+direction*z for x in [-.25,.25]for y in [-.25,.25]for z in [-.5,0.]])
                v0,f0=hull_mesh(vertices);return b.source(v0,f0,'zone-plane-'+letter+'-'+str(fraction))
            low,high=half(.5),half(.78);m=solids.pop(name)
            solids[name+':inner-zone']=m^low;solids[name+':inner-stripe']=(m-low)^high;solids[name+':outer-stripe']=m-high
        stage('Native medullary zones partitioned')
    # Save individual source-provenance cells even if the later interface gate
    # rejects their arrangement; failed experiments never replace the atlas.
    b.save(solids,out/'boolean-inputs.npz');raw=[]
    for name,m in solids.items():
        mm=m.to_mesh64();v=np.asarray(mm.vert_properties)[:,:3]+ORIGIN;f=np.asarray(mm.tri_verts,int)
        mesh=Mesh(v,f,FRAME,source=dict(plan_sha256=digest_file(out/'plan.json'),role=name));path=out/(name.replace(':','-')+'.npz');mesh.save(path)
        q=audit(mesh,require_closed=True,max_pairs=1000000);atomic_json(path.with_suffix('.audit.json'),q)
        row=dict(cell=name,faces=len(f),passed=q['passed'],genera=[c['closed_orientable_genus']for c in q['topology']['components']],failures=q['failures']);raw.append(row);print(row,flush=True)
    if not all(q['passed']for q in raw):
        atomic_json(out/'raw-gate.json',dict(passed=False,cells=raw));return dict(passed=False,cells=raw)
    r=b.assemble(solids,progress=lambda n,t,o:print('shared native charts',n,t,o,round(time.monotonic()-start,1),flush=True));r['positions']+=ORIGIN
    arrays={'positions':r['positions']};rows=[]
    for i,(sid,f)in enumerate(r['surfaces'].items()):arrays['surface_'+str(i)]=f;rows.append(dict(id=sid,array='surface_'+str(i)))
    atomic_npz(out/'arrangement.npz',**arrays)
    atomic_json(out/'result.json',dict(schema=plan['schema'],surfaces=rows,cells=r['cells'],raw_checks=raw,evidence=r['evidence'],seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
    return r


if __name__=='__main__':build()
