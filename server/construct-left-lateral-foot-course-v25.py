"""One measured +4 mm lateral section-transport diagnostic, not installable.

The current fixed endpoint/fascia intersections are explicitly retained. No
radius reduction or obstacle movement is permitted to hide them.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.interpolate import CubicSpline
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_sections import redirect_local_sections,section_metrics,bend_diagnostics,path_tangents
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts,self_contact_summary
from geometry_tools.spatial import BVH
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
OUT=BASE/'left-lateral-proper-transport-v25'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def smooth(x):
    t=np.clip(x,0,1);return 6*t**5-15*t**4+10*t**3

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    contract=json.loads((BASE/'course-contract-v23/contract.json').read_bytes());replay=json.loads((BASE/'course-contract-v23/independent-mouth-array-replay.json').read_bytes());assert replay['passed']
    spec=next(r for r in contract['regions']if r['side']=='l'and r['kind']=='lateral');d=np.load(spec['stations'])
    source=Mesh.load(contract['sourceBindings'][spec['sourcePartId']]['mesh']);rings=d['rings'];v=d['physical_vertices'];f=d['physical_faces'];centers=d['centers'];inverse=d['source_to_physical']
    index=np.arange(len(rings));weight=smooth((index-2)/12)*(1-smooth((index-44)/25))
    delta=np.c_[.004*weight,np.zeros((len(weight),2))]
    target=centers+delta
    # The source/recipient circle profiles agree up to verified 12 nm
    # encoding error. Taper the 1.3456 degree rigid fit through source
    # stations 0..10, keeping every downstream section literal.
    mouthweight=1-smooth(index/10);Rmouth=Rotation.from_rotvec(mouthweight[:,None]*Rotation.from_matrix(d['arch_fit_rotation']).as_rotvec()).as_matrix()
    target+=mouthweight[:,None]*(d['first_ring_target'].mean(0)-centers[0])
    profiles=v[rings][None];redirected,Rcourse=redirect_local_sections(profiles,target,dtype=np.float64)
    transforms=Rcourse@Rmouth
    # Complete mouth rotation is prescribed by the actual arch, and cannot
    # be changed by a finite-difference endpoint-tangent approximation.
    transforms[0]=d['arch_fit_rotation'];target[0]=d['first_ring_target'].mean(0)
    moved=np.einsum('sij,svj->svi',transforms,profiles[0]-centers[:,None])+target[:,None]
    changed=np.any(delta!=0,axis=1)|(mouthweight>0)
    moved[~changed]=profiles[0,~changed]
    encoded=moved.astype('<f4').astype(float);encoded[~changed]=profiles[0,~changed]
    encoded[0]=d['first_ring_target']
    output=v.copy();output[rings]=encoded
    raw=output[inverse];normal=np.array(source.attributes['normals'],copy=True);st=np.full(len(v),-1,int)
    for j,ring in enumerate(rings):st[ring]=j
    for j in np.flatnonzero(changed):
        used=np.flatnonzero(st[inverse]==j);n=source.attributes['normals'][used].astype(float)/32767
        normal[used]=np.clip(np.rint((n@transforms[j].T)*32767),-32767,32767).astype('<i2')
    attrs={k:z.copy()for k,z in source.attributes.items()};attrs['normals']=normal
    mesh=Mesh(raw,source.faces,source.frame,attrs,dict(source.source,diagnosticCandidate=True,sourceGeometrySha256=source.geometry_sha,
        courseField='one 4 mm lateral C2 bump, proper source-section transport',anatomicalAcceptance=False))
    path=OUT/'candidate-ssv.npz';mesh.save(path)
    # Independent source-section and literal-protected-control evidence.
    orth=float(abs(transforms.transpose(0,2,1)@transforms-np.eye(3)).max());det=np.linalg.det(transforms)
    before=section_metrics(profiles);after=section_metrics(encoded[None]);areaerror=float(abs(after['area']-before['area']).max());relative=float((abs(after['area']-before['area'])/before['area']).max())
    diameter=float(abs(after['diameter']-before['diameter']).max());storage=float(np.linalg.norm(encoded-moved,axis=2).max())
    # Perimeter perturbation gives an explicit area roundoff allowance.
    per=np.linalg.norm(profiles[0]-np.roll(profiles[0],-1,axis=1),axis=2).sum(1);bound=2*storage*per+12*storage**2+256*np.finfo(float).eps*before['area'][0]
    assert np.all(abs(after['area'][0]-before['area'][0])<=bound)
    assert orth<1e-12 and det.min()>1-1e-12 and np.array_equal(raw[st[inverse]>=70],source.vertices[st[inverse]>=70])
    assert np.array_equal(mesh.faces,source.faces)and np.array_equal(mesh.attributes['uvs'],source.attributes['uvs'])
    assert np.array_equal(normal[st[inverse]>=70],source.attributes['normals'][st[inverse]>=70])
    # Continuous C2 cubic guide through the actual emitted section centers.
    # Bernstein convex bounds for its derivative/acceleration supply a
    # conservative whole-interval radius lower bound; mesh gates stay separate.
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(centers,axis=0),axis=1))];guide=CubicSpline(s,encoded.mean(1));bounds=[]
    for j,h in enumerate(np.diff(s)):
        aa,bb,cc,_=guide.c[:,j];db=np.array([cc,cc+bb*h,cc+2*bb*h+3*aa*h*h]);forward=target[j+1]-target[j];forward/=np.linalg.norm(forward)
        speed=float((db@forward).min());acc=float(max(np.linalg.norm(2*bb),np.linalg.norm(2*bb+6*aa*h)));radius=speed**2/acc if speed>0 and acc else(1e10 if acc==0 else 0.)
        bounds.append(dict(interval=j,speedLowerBound=speed,accelerationUpperBound=acc,bendRadiusLowerBound_m=radius))
    selfgate,_=self_contact_summary(mesh,max_pairs=250000,max_tests=2000000);assert selfgate['complete_within_pair_budget']
    topology_gate=topology(mesh)
    atomic_npz(OUT/'field.npz',rotation=transforms,source_rings=rings,source_physical=v,source_to_physical=inverse,target_centers=target,encoded_sections=encoded,
        displacement_weight=weight,mouth_weight=mouthweight,source_profile_area=before['area'],candidate_profile_area=after['area'],area_roundoff_bound=bound,
        unchanged_station_mask=~changed,source_faces=source.faces)
    runtime=json.loads(Path(contract['fixedCurrentContext']).read_bytes());base=json.loads((BASE/'left-whole-context-baseline-v24/evidence.json').read_bytes());old={r['partId']:r for r in base['rows']}
    ids=np.array(spec['exactSourceFaceOwnership'],int);region=Mesh(mesh.vertices,mesh.faces[ids],mesh.frame);tree=BVH(region);bb=region.bounds;contacts=[]
    for p in runtime['rows'][0]['rows']:
        if not p['retained']or p['id']==spec['sourcePartId']:continue
        obstacle=Mesh.load(contract['sourceBindings'][p['id']]['mesh']);box=obstacle.bounds
        if np.any(box[1]<bb[0])or np.any(box[0]>bb[1]):gate=None;pairs=np.empty((0,2),int);rels=[]
        else:
            detail,gate=pair_contacts(region,obstacle,0.,first_bvh=tree,max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget'];pairs=detail['pairs'];rels=detail['relation'].tolist()
        prior=old[p['id']]['gate'];beforecross=prior['crossing_face_pairs'] if prior else 0;nowcross=gate['crossing_face_pairs'] if gate else 0
        contacts.append(dict(partId=p['id'],system=p['system'],sourceCrossings=beforecross,candidateCrossings=nowcross,gate=gate,pairs=pairs.tolist(),relations=rels,
            strictCountNonRegression=nowcross<=beforecross))
    arch=Mesh.load(contract['sourceBindings'][spec['recipientPartId']]['mesh']);joinedv=np.vstack([mesh.vertices,arch.vertices]);joinedf=np.vstack([mesh.faces,arch.faces+len(mesh.vertices)])
    joined=Mesh(joinedv,joinedf,mesh.frame);joinedtop=topology(joined);joinedself,_=self_contact_summary(joined,max_pairs=250000,max_tests=2000000);assert joinedself['complete_within_pair_budget']
    # Actual diagnostic sections through unchanged bone and all moved source
    # geometry; no cap is added and original Blender scene stays untouched.
    viewmesh=Mesh(mesh.vertices,mesh.faces[ids],mesh.frame,mesh.attributes,dict(mesh.source,diagnosticFootOnly=True));viewpath=OUT/'candidate-foot.npz';viewmesh.save(viewpath)
    orig=Mesh(source.vertices,source.faces[ids],mesh.frame,source.attributes);origpath=OUT/'source-foot.npz';orig.save(origpath)
    context=['ZA:Calcaneus.l:Bone-5','ZA:Fibularis brevis muscle.l:Tendon','ZA:Fibularis longus muscle.l:Tendon','ZA:Superior fibular retinaculum.l','ZA:Crural fascia.l']
    entries=[dict(mesh=str(viewpath),name='One proper-transport candidate; held for fascia endpoints',color=[.05,.8,.95,1.])]
    for pid in context:entries.append(dict(mesh=contract['sourceBindings'][pid]['mesh'],name='Fixed '+pid,color=[.72,.72,.68,.6]if'Calcaneus'in pid else[.7,.4,.3,.5]))
    center=encoded[:73].mean((0,1));review=dict(frame=mesh.frame,meshes=entries,views=[dict(name='ONE candidate left foot, fixed whole context; lateral',direction=[1,0,0]),
        dict(name='ONE candidate left foot; actual uncapped axial cut',direction=[0,1,0],section=dict(normal=[0,1,0],point=center.tolist(),clear_outer=True))],resolution=1400,wall_seconds=45)
    atomic_json(OUT/'review-spec.json',review);request=prepare(None,review)
    result=dict(schema='one-left-lateral-foot-source-transport-diagnostic/v25',sourceCatalogueSha256=contract['sourceCatalogueSha256'],contractSha256=digest_file(BASE/'course-contract-v23/contract.json'),
        candidate=str(path),candidateSha256=digest_file(path),fieldSha256=digest_file(OUT/'field.npz'),properRotationOrthogonality=orth,minimumProperDeterminant=float(det.min()),
        sourceProfilesRetainedWithinFloat32Rounding=True,maximumRelativeSectionAreaError=relative,maximumDiameterChange_m=diameter,
        maximumSectionPointRoundoff_m=storage,maximumControlDisplacement_m=float(np.linalg.norm(raw-source.vertices,axis=1).max()),
        allFacesAndUVLiteral=True,allStations70OnwardPositionAndNormalsLiteral=True,allOtherMeshesFixed=True,sourceMouthCoordinatesNowExact=True,
        sampledBend=bend_diagnostics(encoded.mean(1),encoded[None]),continuousCubicGuideIntervalBounds=bounds,
        candidateTopology=topology_gate,candidateSelfContact=selfgate,joinedArchOuterTopology=joinedtop,joinedArchOuterSelfContact=joinedself,
        all104CurrentNeighbors=contacts,allNeighborCountsNonRegressing=all(q['strictCountNonRegression']for q in contacts),
        fixedEndpointFasciaIntersectionsRemain=True,acceptedForInstallation=False,finiteWallOrBloodLumenEstablished=False,
        review=request,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],review=request['request_id'],selfCrossings=selfgate['crossing_face_pairs'],
        joinedCrossings=joinedself['crossing_face_pairs'],areaRelative=relative,contacts=[dict(id=q['partId'],before=q['sourceCrossings'],after=q['candidateCrossings'])for q in contacts if q['sourceCrossings']or q['candidateCrossings']])))

if __name__=='__main__':main()
