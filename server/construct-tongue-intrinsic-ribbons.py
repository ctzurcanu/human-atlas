"""Original broad intrinsic fascicle ribbons, local candidate only.

The article's subdorsal longitudinal sheet, ventrolateral longitudinal bellies
and interleaved transverse/vertical lamellae guide the layout. Source-specific
internal borders and fibre thickness are not available; these remain explicitly
approximate representative anatomy, rather than copied donor segmentation.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time
import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.audit import topology

module_spec = importlib.util.spec_from_file_location('tongue_fascicles', ROOT/'server/construct-tongue-intrinsics.py')
base = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(base)


class UpperBodyEnvelope(base.Envelope):
    """Own only the dorsal interior interval; never bridge an exterior gap."""
    def __init__(self,mesh):
        super().__init__(mesh)
        self.multiple_interval_queries=0

    def interval(self,x,z):
        uv=np.einsum('nij,nj->ni',self.inverse,np.array([x,z])-self.base[self.valid])
        bary=np.c_[1-uv.sum(1),uv]
        hit=np.all(bary>=-1e-10,axis=1)
        yy=np.sort(np.einsum('ni,ni->n',bary[hit],self.tri[self.valid[hit],:,1]))
        yy=yy[np.r_[True,np.diff(yy)>1e-10]]
        if len(yy)<2 or len(yy)%2:
            raise ValueError(f'Unpaired native interior boundaries at {(x,z)}: {yy}')
        self.multiple_interval_queries+=len(yy)>2
        return yy[-2],yy[-1]


class FacetedUpperSurface:
    """Literal upper Y surface of an emitted neighbour, without smoothing it."""
    def __init__(self, mesh):
        tri=mesh.triangles
        xz=tri[:,:,[0,2]]
        matrix=np.stack([xz[:,1]-xz[:,0],xz[:,2]-xz[:,0]],axis=-1)
        valid=np.abs(np.linalg.det(matrix))>1e-15
        self.tri=tri[valid]
        self.origin=xz[valid,0]
        self.inverse=np.linalg.inv(matrix[valid])

    def upper(self,x,z):
        uv=np.einsum('nij,nj->ni',self.inverse,np.array([x,z])-self.origin)
        bary=np.c_[1-uv.sum(1),uv]
        hit=np.all(bary>=-1e-10,axis=1)
        if not np.any(hit):
            return None
        return float(np.max(np.einsum('ni,ni->n',bary[hit],self.tri[hit,:,1])))


def taper(u):
    return .3 + .7 * np.sin(np.pi*u)**.55


def core_limits(env, x, z):
    bottom, top = env.interval(x, z)
    height = top-bottom
    # A conservative corridor between the broad longitudinal fascicle regions.
    lower = bottom + min(.0055, .28*height) + min(.0026, .15*height) + .0008
    upper = top - min(.0031, .25*height) - min(.00125, .12*height) - .0008
    if upper <= lower:
        raise ValueError(f'No positive representative interlamellar corridor at {(x,z)}')
    return lower, upper


def vertical_limits(env,x,z):
    lower,upper=core_limits(env,x,z)
    owners=getattr(env,'inferior_owners',None)
    if owners is not None:
        boundary=owners['left' if x>=0 else 'right'].upper(x,z)
        if boundary is not None:
            # Representative interstitial spacing, not a measured histological gap.
            lower=max(lower,boundary+.00075)
    if upper<=lower:
        raise ValueError(f'Neighbour clipping leaves no vertical fascicle corridor at {(x,z)}')
    return lower,upper


def x_for(env, fraction, z):
    lo, hi = env.width(z)
    return (lo+hi)/2 + fraction*(hi-lo)/2


def finite_sweep(point, centre, steps, radial):
    vv=[]; ff=[]
    for u in np.linspace(0,1,steps):
        for theta in np.arange(radial)*2*np.pi/radial:
            vv.append(point(u, theta))
    vv.extend([centre(0.), centre(1.)])
    for j in range(steps-1):
        for k in range(radial):
            a=j*radial+k; b=j*radial+(k+1)%radial
            c=(j+1)*radial+k; d=(j+1)*radial+(k+1)%radial
            ff.extend([[a,b,c],[b,d,c]])
    c0,c1=len(vv)-2,len(vv)-1
    for k in range(radial):
        ff.extend([[c0,(k+1)%radial,k],[c1,(steps-1)*radial+k,(steps-1)*radial+(k+1)%radial]])
    v,f=np.array(vv),np.array(ff)
    t=v[f]-v.mean(0)
    if np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()<0:
        f=f[:,::-1]
    return v,f


def longitudinal(env, fraction, side=1., superior=True):
    half_fraction=.105 if superior else .20
    def point(u,theta):
        z=.030 + .0502*u
        f=fraction+side*half_fraction*np.cos(theta)*taper(u)
        x=x_for(env,f,z)
        bottom,top=env.interval(x,z); h=top-bottom
        y=(top-min(.0031,.25*h)+min(.00125,.12*h)*np.sin(theta)*taper(u)
           if superior else bottom+min(.0055,.28*h)+min(.0026,.15*h)*np.sin(theta)*taper(u))
        return [x,y,z]
    def centre(u):
        z=.030+.0502*u; x=x_for(env,fraction,z); bottom,top=env.interval(x,z); h=top-bottom
        return [x,top-min(.0031,.25*h) if superior else bottom+min(.0055,.28*h),z]
    return finite_sweep(point,centre,30,10)


def transverse(env, z0, side, level):
    def point(u,theta):
        z=z0+.0013*np.cos(theta)*taper(u)
        x=x_for(env,side*(.035+.725*u),z)
        floor,ceil=core_limits(env,x,z)
        y=floor+level*(ceil-floor)+min(.0014,.14*(ceil-floor))*np.sin(theta)*taper(u)
        return [x,y,z]
    def centre(u):
        x=x_for(env,side*(.035+.725*u),z0); floor,ceil=core_limits(env,x,z0)
        return [x,floor+level*(ceil-floor),z0]
    return finite_sweep(point,centre,12,8)


def vertical(env,z0,fraction):
    def point(u,theta):
        z=z0+.0014*np.sin(theta)*taper(u)
        x=x_for(env,fraction,z)+.00085*np.cos(theta)*taper(u)
        floor,ceil=vertical_limits(env,x,z)
        return [x,floor+u*(ceil-floor),z]
    def centre(u):
        x=x_for(env,fraction,z0); floor,ceil=vertical_limits(env,x,z0)
        return [x,floor+u*(ceil-floor),z0]
    return finite_sweep(point,centre,10,8)


def main(out,repair_from=None):
    began=time.monotonic()
    assert (ROOT/'reports/ta98-tongue-inventory.json').exists()
    a,native,reader=base.source()
    assert topology(native)['closed_oriented_manifold']
    env=UpperBodyEnvelope(native)
    previous=None; repair_plan=None
    if repair_from is not None:
        previous=json.loads((repair_from/'construction.json').read_text())
        prior_check=json.loads((repair_from/'checked-geometry.json').read_text())
        assert previous['source_catalogue_sha256']==digest_file(base.SOURCE)
        assert prior_check['construction_sha256']==digest_file(repair_from/'construction.json')
        assert prior_check['cross_muscle_defect_pairs']==14 and not prior_check['passed']
        for row in previous['meshes']:
            assert digest_file(row['mesh'])==row['sha256']
        env.inferior_owners={side:FacetedUpperSurface(Mesh.load(repair_from/('inferior-longitudinal-'+side+'.npz')))
                             for side in ['left','right']}
        repair_plan=dict(schema='tongue-intrinsic-local-repair/v1',
            previous_construction_sha256=digest_file(repair_from/'construction.json'),
            previous_check_sha256=digest_file(repair_from/'checked-geometry.json'),
            preserved_defect_witnesses=str(repair_from/'defect-witnesses.json'),
            preserved_defect_witnesses_sha256=digest_file(repair_from/'defect-witnesses.json'),
            original_defect_face_pairs=14,changed_groups=['vertical-left','vertical-right'],
            method='Clip only the inferior ends of vertical ribbons against actual emitted inferior-longitudinal upper triangles',
            illustrative_clearance_m=.00075,thresholds_relaxed=False,prior_rejected_evidence_preserved=True)
        out.mkdir(parents=True,exist_ok=True)
        atomic_json(out/'repair-plan.json',repair_plan)
    groups={
        'superior-longitudinal':[longitudinal(env,f) for f in [-.78,-.52,-.26,0.,.26,.52,.78]],
    }
    for side,sign in [('left',1.),('right',-1.)]:
        groups['inferior-longitudinal-'+side]=[longitudinal(env,sign*.62,sign,False)]
        groups['transverse-'+side]=[transverse(env,z,sign,level) for z in np.linspace(.034,.074,6) for level in [.30,.70]]
        groups['vertical-'+side]=[vertical(env,z,sign*f) for z in np.linspace(.038,.070,5) for f in [.22,.62]]
    terms={'superior-longitudinal':'A05.1.04.106','inferior-longitudinal':'A05.1.04.107','transverse':'A05.1.04.108','vertical':'A05.1.04.109'}
    colors={'superior-longitudinal':[.72,.30,.29,1.],'inferior-longitudinal':[.68,.34,.31,1.],
            'transverse':[.78,.43,.36,1.],'vertical':[.62,.28,.32,1.]}
    rows=[]; out.mkdir(parents=True,exist_ok=True)
    origin=np.array([0.,1.625,.055])
    limits=(
        'Original approximate fascicle ribbons guided by primary anatomical sections, not a measured donor segmentation. '
        'The broad subdorsal longitudinal group and curved ventrolateral bellies outline regional extents; the interdigitated '
        'transverse/vertical lamellae are representative, with illustrative spacing, count, thickness and terminal caps. '
        'The model is not a complete voxel/tissue partition or proof of muscular attachments. Native external-muscle '
        'envelopes may include interdigitating territories; their exclusive ownership is not established.'
    )
    for key,bundles in groups.items():
        vv=[]; ff=[]; components=[]; start=0
        for v,f in bundles:
            vv.append(v); ff.append(f+start); components.append(dict(vertex_start=start,vertex_count=len(v),triangle_count=len(f)))
            start+=len(v)
        v=np.concatenate(vv); f=np.concatenate(ff)
        emitted=(v-origin).astype(np.float32).astype(float)+origin
        n=igl.per_vertex_normals(emitted,f)
        kind=key.replace('-left','').replace('-right','')
        mesh=Mesh(emitted,f,base.FRAME,{'normals':n},dict(kind='original-intrinsic-fascicle-ribbons',
            source_catalogue_sha256=digest_file(base.SOURCE),source_tongue_geometry_sha256=native.geometry_sha,
            primary_reference=base.PAPER,source_mesh_vertices_copied=False,proprietary_geometry_copied=False,
            representation='approximate-fascicle-ribbons',complete_muscle_volume_segmentation=False,
            construction_origin_m=origin.tolist(),components=components,limitation=limits))
        target=out/(key+'.npz'); mesh.save(target)
        rows.append(dict(key=key,term=terms[kind],mesh=str(target),sha256=digest_file(target),geometry_sha256=mesh.geometry_sha,
            partId='TA98TONGUE:male:intrinsic:'+key,bundles=len(bundles),triangles=len(f),vertices=len(emitted),
            color=colors[kind],representation='approximate-fascicle-ribbons',coverageLimitation=limits,
            source_corners_copied=False,anatomical_acceptance=False))
    if previous is not None:
        previous_rows={r['key']:r for r in previous['meshes']}
        deltas=[]
        for row in rows:
            old=Mesh.load(previous_rows[row['key']]['mesh']); new=Mesh.load(row['mesh'])
            assert np.array_equal(old.faces,new.faces)
            displacement=np.linalg.norm(new.vertices-old.vertices,axis=1)
            changed=int(np.count_nonzero(displacement>0))
            if not row['key'].startswith('vertical-'):
                assert new.geometry_sha==old.geometry_sha and np.array_equal(new.attributes['normals'],old.attributes['normals'])
            deltas.append(dict(key=row['key'],changed_vertices=changed,max_displacement_m=float(displacement.max()),
                              preserved_geometry=new.geometry_sha==old.geometry_sha))
        atomic_json(out/'repair-deltas.json',dict(schema='tongue-intrinsic-local-repair-deltas/v1',meshes=deltas))
    report=dict(schema='tongue-original-intrinsic-ribbons/v1',source_catalogue_sha256=digest_file(base.SOURCE),
        source_tongue_geometry_sha256=native.geometry_sha,source_bank_sha256=native.source['bank_sha256'],frame=base.FRAME,
        primary_reference=base.PAPER,source_geometry_modified=False,installed=False,anatomical_acceptance=False,
        meshes=rows,total_meshes=len(rows),total_triangles=sum(r['triangles']for r in rows),
        total_raw_geometry_bytes=sum(r['vertices']*18+r['triangles']*12 for r in rows),
        construction_script_sha256=digest_file(__file__),checker_script_sha256=digest_file(base.__file__),
        native_ray_ownership='Upper interior interval beneath the dorsum; separate lower intervals and exterior gaps are never bridged',
        repair=repair_plan,
        multiple_native_interval_queries=env.multiple_interval_queries,
        all_canonical_muscle_terms=list(terms.values()),complete_muscle_volume_segmentation=False,
        state='constructed-awaiting-independent-audit-and-rendered-review',remaining_limits=[limits],elapsed_seconds=time.monotonic()-began)
    assert report['total_meshes']<=8 and report['total_triangles']<=15000
    atomic_json(out/'construction.json',report)
    source=json.loads((ROOT/'.local-models/geometry-workbench/tongue-completion/source-v1/source-audit.json').read_text())
    context=[dict(mesh=p['meshSavePath'],name='Fixed native context: '+p['name'],color=[.6,.7,.8,.12])for p in source['contextParts'][:5]]
    review=dict(frame=base.FRAME,meshes=[dict(mesh=r['mesh'],name=r['key'],color=r['color'])for r in rows],resolution=640,wall_seconds=60,
        views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[1,0,0]),dict(name='top',direction=[0,1,0]),
               dict(name='right-section',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.625,.055],clear_outer=True))])
    atomic_json(out/'review-spec.json',review)
    atomic_json(out/'review-with-context.json',{**review,'meshes':review['meshes']+context})
    print(json.dumps({k:report[k]for k in ['total_meshes','total_triangles','total_raw_geometry_bytes','elapsed_seconds']}),flush=True)
    base.OUT=out
    base.check()


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--repair-from',type=Path);args=ap.parse_args()
    out=args.output.resolve()
    if not out.is_relative_to(ROOT/'.local-models/geometry-workbench/tongue-completion'):
        raise ValueError('Keep original local candidates in the tongue workbench')
    main(out,args.repair_from.resolve() if args.repair_from else None)
