"""Source-bound canonical renal completion; never edits installed anatomy."""
from pathlib import Path
import json,time
import numpy as np
from scipy.spatial import ConvexHull
from .common import ROOT,CACHE,atomic_json,atomic_npz,digest_file
from .model import Mesh
from .audit import audit
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

FRAME='model:local-ta98:assembled-metres'
OUT=CACHE/'renal-completion'


def bound(name):
    p=ROOT/'reports'/(name+'.json');r=json.loads(p.read_text());archive=ROOT/r['candidate']
    if digest_file(archive)!=r['candidateSha256']:raise ValueError('Changed source candidate: '+name)
    return r,np.load(archive,allow_pickle=False)


def arrays(solid):
    m=solid.to_mesh64();return np.asarray(m.vert_properties)[:,:3].copy(),np.asarray(m.tri_verts,dtype='int64').copy()


def hull_solid(points):
    p=np.unique(np.asarray(points,float),axis=0);h=ConvexHull(p);f=h.simplices.copy();q=p[f];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);wrong=np.einsum('ij,ij->i',n,h.equations[:,:3])<0;f[wrong]=f[wrong][:,[0,2,1]]
    return make_solid(p,f)[0]


def trial(version='428'):
    out=OUT/('master-v'+version)
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);began=time.monotonic()
    mr,mz=bound('medullary-collecting-cells-v300');fr,fz=bound('renal-ureter-fluid-v317');cr,cz=bound('renal-cortex-native-v353');rr,rz=bound('renal-columns-native-v388')
    source_names=['medullary-collecting-cells-v300','renal-ureter-fluid-v317','renal-cortex-native-v353','renal-columns-native-v388']
    plan=dict(schema='renal-master-v'+version,inventory_sha256=digest_file(ROOT/'TA98-KIDNEY-PARTS.md'),sources={n:digest_file(ROOT/'reports'/(n+'.json'))for n in source_names},
              master='Male left; installed master and original sources retained',capsule='One convex hilar gateway joins existing four collecting/capsule crossings and extends medially outside the kidney. Classifier is not an emitted solid.',
              cortex='Union of registered source cortical and column guides, clipped by the same native inner capsule and excluding actual medullary/collecting occupancy.',
              construction_simplification_m=.00002,unchanged_native_outer=True,anatomical_acceptance=False,installed=False)
    atomic_json(out/'plan.json',plan)
    K=make_solid(mz['native_kidney_v'],mz['native_kidney_f'])[0];I=make_solid(mz['native_inner_v'],mz['native_inner_f'])[0];O=make_solid(fz['positions'],fz['outerClassifier'])[0]
    cap=K-I;crossing=cap^O;ports=[]
    for s in crossing.decompose():
        v,f=arrays(s);ports.append(v)
    if len(ports)!=4:raise ValueError('The declared four original collecting perforations no longer match the master')
    p=np.concatenate(ports);portal=hull_solid(np.r_[p,p+[-.07,0,0]])
    corrected_cap=cap-portal
    M_surfaces={s['id']:mz[s['facesArray']]for s in mr['surfaces']};mf=[]
    for cell,uses in mr['cells'].items():
        if cell=='collecting-wall':continue
        for use in uses:
            if use.get('neighbor')in mr['cells']and use['neighbor']!='collecting-wall':continue
            fs=M_surfaces[use['surface']];mf.append(fs if use['side']==1 else fs[:,[0,2,1]])
    M=make_solid(mz['positions'],np.concatenate(mf))[0]
    Cref=make_solid(cz['sourceCortex_v'],cz['sourceCortex_f'])[0].simplify(.00002)
    Rref=make_solid(rz['sourceCommonEncoded_v'],rz['sourceFaces'])[0].simplify(.00002)
    blocked=M+O
    full=((Cref+Rref)^I)-blocked
    columns,peripheral=full.split(Rref)
    results=[];solids={'capsule':corrected_cap,'cortex-whole':full,'columns':columns,'peripheral':peripheral,'sinus-unclassified':I-(M+O+full)}
    for name,solid in solids.items():
        v,f=arrays(solid);mesh=Mesh(v,f,FRAME,source=dict(plan_sha256=digest_file(out/'plan.json'),role=name));path=out/(name+'.npz');mesh.save(path)
        check=audit(mesh,max_pairs=750000);atomic_json(out/(name+'-audit.json'),check)
        row=dict(role=name,path=str(path),sha256=digest_file(path),geometry_sha256=mesh.geometry_sha,faces=len(f),occupied_components=len(solid.decompose()),volume_m3=float(solid.volume()),passed=check['passed'],boundary_genera=[g['closed_orientable_genus']for g in check['topology']['components']])
        results.append(row);print(row,flush=True)
    atomic_npz(out/'source-fields.npz',outer_v=mz['native_kidney_v'],outer_f=mz['native_kidney_f'],inner_v=mz['native_inner_v'],inner_f=mz['native_inner_f'],portal_points=p)
    report=dict(schema=plan['schema'],parts=results,plan_sha256=digest_file(out/'plan.json'),one_capsule_hilum=results[0]['boundary_genera']==[0],
                cortex_partition_volume_residual_m3=float(columns.volume()+peripheral.volume()-full.volume()),elapsed_seconds=time.monotonic()-began,installed=False,anatomical_acceptance=False)
    atomic_json(out/'result.json',report);return report


if __name__=='__main__':trial()


def common_field_master(version='433',spacing=.002):
    """One arrangement, including the real collecting continuation.

    The declared sinus envelope is a construction boundary around the calyceal
    tree and the actual hilar aperture, never a disk across a urine passage.
    Its extent and regional cortical depth require rendered anatomical review.
    """
    import point_cloud_utils as pcu
    from .implicit_cells import grid_tetrahedra,Arrangement
    out=OUT/('master-v'+version)
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);began=time.monotonic()
    mr,mz=bound('medullary-collecting-cells-v300');fr,fz=bound('renal-ureter-fluid-v317')
    z=np.load(OUT/'master-v428/source-fields.npz');ports=z['portal_points']
    cy,cz=(ports[:,1:].min(0)+ports[:,1:].max(0))/2
    ry,rz=np.ptp(ports[:,1:],axis=0)*.58*1.15;xend=ports[:,0].max()+.0012
    outer=fz['positions'][np.unique(fz['outerClassifier'])]
    lo=np.minimum(mz['native_kidney_v'].min(0),outer.min(0))-.0031
    hi=np.maximum(mz['native_kidney_v'].max(0),outer.max(0))+.0031
    def distance_at(points,v,f):
        used,inv=np.unique(f,return_inverse=True)
        return pcu.signed_distance_to_mesh(points,np.ascontiguousarray(v[used]),np.ascontiguousarray(inv.reshape(-1,3),dtype=np.int32))[0]
    def domain_distance(points):
        return np.minimum(distance_at(points,mz['native_kidney_v'],mz['native_kidney_f']),distance_at(points,fz['positions'],fz['outerClassifier']))
    points,tets,grid=grid_tetrahedra(lo,hi,spacing,domain_distance=domain_distance)
    def sdf(v,f):
        return distance_at(points,v,f)
    K=sdf(mz['native_kidney_v'],mz['native_kidney_f']);I=sdf(mz['native_inner_v'],mz['native_inner_f'])
    O=sdf(fz['positions'],fz['outerClassifier']);F=sdf(fz['positions'],fz['innerClassifier'])
    P=np.maximum(np.sqrt(((points[:,1]-cy)/ry)**2+((points[:,2]-cz)/rz)**2)-1,(points[:,0]-xend)/max(ry,rz))*.01
    sf={s['id']:mz[s['facesArray']]for s in mr['surfaces']};names=[];fields=[]
    for name,uses in mr['cells'].items():
        if name=='collecting-wall':continue
        faces=np.concatenate([sf[u['surface']]if u['side']==1 else sf[u['surface']][:,[0,2,1]]for u in uses])
        fields.append(sdf(mz['positions'],faces));names.append(name)
    # Common global field order is part of the construction identity.
    values=np.array([np.minimum(K,O),K,I,P,O,F,O*.67+F*.33,O*.33+F*.67,O-.0015,I+.0066,*fields]).T
    field_names=['domain','kidney','inner-capsule','hilum','collecting-outer','urine','coat-outer-third','coat-inner-third','sinus-envelope','peripheral-cortex-depth',*names]
    atomic_npz(out/'sampled-fields.npz',points=points,tets=tets,values=values)
    plan=dict(schema='renal-common-fields/v1',master='male-left',spacing_m=spacing,grid=grid,field_order=field_names,
              inventory_sha256=digest_file(ROOT/'TA98-KIDNEY-PARTS.md'),sources={n:digest_file(ROOT/'reports'/(n+'.json'))for n in ['medullary-collecting-cells-v300','renal-ureter-fluid-v317']},
              exact_source_guides=True,output_is_remeshed=True,cortical_peripheral_depth_m=.0066,sinus_envelope_offset_m=.0015,
              fixed_skeleton=True,installed=False,anatomical_acceptance=False)
    atomic_json(out/'plan.json',plan)
    # Picking the deepest overlapping medullary guide requires pairwise fields
    # at overlaps; source cells are a disjoint arrangement, so use a stable
    # explicit tissue priority instead. Their genuine common zero faces remain
    # reciprocal; no independently generated boolean union is used.
    def classify(v):
        if v[1]<0 and v[2]>=0 and v[3]>=0:return 'capsule'
        if v[1]<0:
            inside=np.flatnonzero(v[10:]<0)
            if len(inside):return names[int(inside[0])]
        if v[4]<0:
            if v[5]<0:return None
            return 'collecting-adventitia'if v[6]>=0 else 'collecting-muscularis'if v[7]>=0 else 'collecting-mucosa'
        if v[2]<0 and v[3]>=0 and v[8]>=0:return 'cortex'if v[9]>=0 else 'columns'
        return None
    print('field samples',len(points),'tetrahedra',len(tets),'seconds',round(time.monotonic()-began,2),flush=True)
    A=Arrangement(points,values,classify,vertex_snap_m=1e-12)
    p,s,c=A.build(tets,skip_fields=[0],progress=lambda n,t,v,l:print('arrangement',n,t,v,l,round(time.monotonic()-began,1),flush=True))
    arrays={'positions':p};rows=[]
    for i,(sid,f)in enumerate(s.items()):arrays['surface_'+str(i)]=f;rows.append(dict(id=sid,array='surface_'+str(i)))
    atomic_npz(out/'arrangement.npz',**arrays)
    checks=[]
    for name,uses in c.items():
        f=np.concatenate([s[u['surface']]if u['side']==1 else s[u['surface']][:,[0,2,1]]for u in uses]);mesh=Mesh(p,f,FRAME,source=dict(plan_sha256=digest_file(out/'plan.json'),role=name));mesh.save(out/(name.replace(':','-')+'.npz'))
        q=audit(mesh,require_closed=True,max_pairs=750000);atomic_json(out/(name.replace(':','-')+'-audit.json'),q)
        row=dict(cell=name,faces=len(f),passed=q['passed'],genera=[g['closed_orientable_genus']for g in q['topology']['components']],badlinks=len(q['topology']['invalid_vertex_links']),crossings=q['self_intersections']['defect_pairs']);checks.append(row);print(row,flush=True)
    result=dict(schema=plan['schema'],surfaces=rows,cells=c,checks=checks,seconds=time.monotonic()-began,passed=all(q['passed']for q in checks),installed=False,anatomical_acceptance=False)
    atomic_json(out/'result.json',result);return result
