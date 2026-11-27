"""One bounded intralingual-only extrinsic-envelope intersection candidate.

The computational tongue envelope is not used to cut bony origins. A declared
floor-of-mouth classifier preserves muscle below the mucosal reflection, while
the exposed intralingual portion fits inside the inner mucosal carrier.
"""
from __future__ import annotations
import argparse,importlib.util,time,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,provenance
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts

spec=importlib.util.spec_from_file_location('exposed_mucosa',Path(__file__).with_name('propose-tongue-mucosal-patch.py'))
patch=importlib.util.module_from_spec(spec);spec.loader.exec_module(patch)
WORK=patch.WORK

def field_above_actual_facets(mesh,error=.00002):
    """Upper-bound H(X,Z)-Y over every actual finite-wall triangle.

    For each separable C1 smoothstep term, the interpolation error is bounded
    by max|H''| range²/8. Adaptive subdivision tightens this certificate;
    triangles and all acceptance predicates themselves are unchanged.
    """
    triangles=mesh.triangles;lower=float(np.max(patch.reflection_height(triangles.reshape(-1,3))-triangles[:,:,1].ravel()))
    pending=[(t,0)for t in triangles];upper=-np.inf;tested=0;maximum_depth=0
    while pending:
        t,depth=pending.pop();tested+=1;maximum_depth=max(maximum_depth,depth)
        if tested>100000:raise ValueError('Reflection-facet certificate work budget exceeded')
        values=patch.reflection_height(t)-t[:,1];sample=float(values.max());lower=max(lower,sample)
        remainder=.002*6/.016**2*np.ptp(t[:,0])**2/8+.0025*6/.017**2*np.ptp(t[:,2])**2/8
        bound=sample+remainder
        if bound<=lower+error:upper=max(upper,bound);continue
        if depth>=12:raise ValueError('Reflection-facet certificate depth budget exceeded')
        a,b,c=t;ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
        pending.extend((q,depth+1)for q in [np.array([a,ab,ca]),np.array([ab,b,bc]),np.array([ca,bc,c]),np.array([ab,bc,ca])])
    return dict(maximumActualFacetOvershootLowerM=lower,maximumActualFacetOvershootUpperM=float(upper),
                requestedErrorM=error,testedSubtriangles=tested,maximumDepth=maximum_depth,
                method='C1-smoothstep global Hessian bound plus adaptive subdivision; actual emitted facets')

def floor_classifier(actual_facet_bound):
    # The setback is derived from actual rim facets, plus a proven terrain
    # interpolation bound and an explicit 50 µm connective clearance.
    xs=np.linspace(-.06,.06,61);zs=np.linspace(-.03,.13,81)
    xx,zz=np.meshgrid(xs,zs,indexing='ij');p=np.c_[xx.ravel(),np.zeros(xx.size),zz.ravel()]
    interpolation_bound=(.002*6/.016**2*np.diff(xs).max()**2+.0025*6/.017**2*np.diff(zs).max()**2)/8
    setback=actual_facet_bound+interpolation_bound+.00005
    if setback>=.0003:raise ValueError('A certified classifier would alter protected source origins; stop rather than trim them')
    p[:,1]=patch.reflection_height(p)-setback
    top=p.copy();bottom=p.copy();bottom[:,1]=1.50;n=len(p);v=np.r_[top,bottom];f=[]
    for i in range(len(xs)-1):
        for j in range(len(zs)-1):
            a=i*len(zs)+j;b=a+len(zs);c=b+1;d=a+1
            f.extend([[a,d,c],[a,c,b],[a+n,c+n,d+n],[a+n,b+n,c+n]])
    ring=np.r_[np.arange(len(zs)),np.arange(1,len(xs))*len(zs)+len(zs)-1,
               (len(xs)-1)*len(zs)+np.arange(len(zs)-2,-1,-1),np.arange(len(xs)-2,0,-1)*len(zs)]
    for a,b in zip(ring,np.roll(ring,-1)):f.extend([[b,a,a+n],[b,a+n,b+n]])
    f=np.asarray(f,int);solid,_=make_solid(v,f)
    # Explicit bound from the maximum smoothstep curvature: interpolation
    # error ≤ max|h''|*step²/8 independently in the separable X and Z fields.
    return solid,float(interpolation_bound),float(setback)

def oriented_triangle(tri):
    q=[tuple(p)for p in tri];return min(tuple(q[i:]+q[:i])for i in range(3))

def decoded_candidate(result,original,original_id,frame):
    raw=result.to_mesh64();p=np.asarray(raw.vert_properties,float);f=np.asarray(raw.tri_verts,int)
    origin,parent=provenance(raw);tri=p[f,:3];normals=np.zeros_like(tri)
    properties=p[f,3:6]if p.shape[1]>=6 else np.zeros_like(tri)
    source=origin==original_id;length=np.linalg.norm(properties,axis=2)
    normals[source]=np.divide(properties[source],length[source,:,None],out=np.zeros_like(properties[source]),where=length[source,:,None]>0)
    fallback=(~source)|np.any(length==0,axis=1)
    geometric=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);geometric/=np.linalg.norm(geometric,axis=1)[:,None]
    normals[fallback]=geometric[fallback,None,:]
    words=np.rint(np.clip(normals,-1,1)*32767).astype('<i2')
    original_keys={oriented_triangle(t):i for i,t in enumerate(original.triangles)}
    literal=[]
    for i,t in enumerate(tri):
        old=original_keys.get(oriented_triangle(t))
        if old is not None:
            literal.append(old)
            for corner,point in enumerate(t):
                match=np.flatnonzero(np.all(original.triangles[old]==point,axis=1))
                words[i,corner]=original.attributes['normals'][original.faces[old,match[0]]]
    world=tri.reshape(-1,3);origin_m=np.array([0.,1.625,0.]);stored=(world-origin_m).astype('<f4');emitted=stored.astype(float)+origin_m
    out=Mesh(emitted,np.arange(len(world)).reshape(-1,3),frame,dict(normals=words.reshape(-1,3)),
        dict(derived_from_geometry=original.geometry_sha,method='intralingual-envelope-intersection; fixed below-reflection origin/courses'))
    return out,dict(literalSourceFacesBeforeEncoding=sorted(set(literal)),
        sourceFaceCount=len(original.faces),generatedIntersectionFaces=int((~source).sum()),
        maximumEncodingErrorM=float(np.linalg.norm(emitted-world,axis=1).max()),geometryOrigin=origin_m.tolist()),stored

def main(output):
    began=time.monotonic();output=output.resolve()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local tongue stage')
    inputs=patch.load_sources();whole=inputs['whole'][1];normal=whole.attributes['normals'].astype(float);normal/=np.linalg.norm(normal,axis=1)[:,None]
    thick=.00035+.00007*np.clip(normal[:,1],0,1)-.00005*np.clip(-normal[:,1],0,1)
    # 50 µm illustrative connective clearance, explicitly separate from the
    # mucosal thickness. This is not a microscopic insertion/enthesis model.
    inner=whole.vertices-(thick+.00005)[:,None]*normal
    mucosa_path=WORK/'exposed-mucosa-v1/packed/TA98TONGUE_male_mucosa.npz';mucosa=Mesh.load(mucosa_path)
    facet_certificate=field_above_actual_facets(mucosa)
    carrier,_=make_solid(inner,whole.faces);floor,interpolation_bound,setback=floor_classifier(facet_certificate['maximumActualFacetOvershootUpperM']);allowed=carrier+floor
    output.mkdir(parents=True);rows=[]
    bones={name:Mesh.load(WORK/f'source-v1/meshes/native-{i:02}.npz')for name,i in [('mandible',5),('hyoid-body',13),('lesser-left',14),('lesser-right',15),('greater-left',16),('greater-right',17)]}
    for name in ['GG-left','GG-right','HG-left','HG-right']:
        if time.monotonic()-began>50:raise TimeoutError('Stopped before 60-second construction deadline')
        path,original=inputs[name]
        host,oid=make_solid(original.vertices,original.faces,properties=original.attributes['normals'][original.faces])
        result=host^allowed
        candidate,retained,stored=decoded_candidate(result,original,oid,whole.frame)
        dest=output/(name+'.npz');candidate.save(dest)
        atomic_npz(output/(name+'-stored.npz'),positions=stored,faces=candidate.faces,geometryOrigin=np.array([0.,1.625,0.]),
            frame=np.array(whole.frame),attribute_normals=candidate.attributes['normals'])
        checks=audit(candidate,require_closed=True,max_pairs=100000,max_tests=1000000)
        arrays,cover=pair_contacts(candidate,mucosa,0.,max_pairs=100000,max_tests=1000000)
        atomic_npz(output/(name+'-mucosa-contact.npz'),**arrays)
        bone_rows={}
        for bone_name,bone in bones.items():
            _,s=pair_contacts(candidate,bone,0.,max_pairs=100000,max_tests=1000000);bone_rows[bone_name]=s
        protected=np.all(original.vertices[original.faces,1]<patch.reflection_height(original.vertices[original.faces].reshape(-1,3)).reshape(-1,3)-.0003,axis=1)
        retained_set=set(retained['literalSourceFacesBeforeEncoding']);protected_ids=set(np.flatnonzero(protected))
        if not protected_ids<=retained_set:raise ValueError('Protected below-reflection source face changed')
        retained.update(protectedBelowReflectionFaces=len(protected_ids),protectedBelowReflectionFacesLiteralBeforeEncoding=len(protected_ids&retained_set),
                        originalVolumeMl=float(host.volume()*1e6),candidateVolumeMl=float(result.volume()*1e6))
        row=dict(name=name,source=str(path),sourceSha256=digest_file(path),candidate=str(dest),candidateSha256=digest_file(dest),
                 preservation=retained,checks=checks,mucosaContact=cover,fixedBoneContacts=bone_rows)
        rows.append(row);atomic_json(output/(name+'-audit.json'),row)
    passed=all(row['checks']['passed']and row['mucosaContact']['crossing_face_pairs']==0 and row['mucosaContact']['coplanar_overlap_face_pairs']==0 for row in rows)
    report=dict(schema='tongue-intralingual-mucosal-fit/v1',candidatePassedNumericalSurfaceGate=passed,
        sourceWhole=str(inputs['whole'][0]),sourceWholeSha256=digest_file(inputs['whole'][0]),mucosa=str(mucosa_path),
        method='Intersection with inner tongue envelope OR below-reflection extralingual space; no individual muscle-shaped holes',
        fixedBoneGeometry=True,mucosalThicknessUnchanged=True,intralingualConnectiveClearanceM=.00005,
        floorClassifierSetbackM=setback,floorInterpolationErrorBoundM=interpolation_bound,actualMucosaFieldCertificate=facet_certificate,
        rows=rows,anatomicalAcceptance=False,installed=False,elapsedSeconds=time.monotonic()-began,
        remaining=['Actual surface/section review', 'Native HG bony-origin penetration predates this fit and remains a separate attachment defect',
                   'Coupled aponeurosis/intrinsic and future palatoglossal interfaces must be tested before installation'])
    atomic_json(output/'candidate.json',report)
    print(dict(passed=passed,elapsed=report['elapsedSeconds'],rows=[dict(name=r['name'],self=r['checks']['passed'],mucosa=r['mucosaContact']['crossing_face_pairs'],boneCrossings=sum(x['crossing_face_pairs']for x in r['fixedBoneContacts'].values()))for r in rows]),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
