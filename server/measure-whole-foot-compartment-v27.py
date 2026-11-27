"""Read-only fascia/skin ray witnesses and nerve relations of foot sources.

Ray intervals use actual open source surfaces and source section support.
They are sampled geometric opportunities, not occupied-volume certificates
or a selected rigid transformation.
"""
from pathlib import Path
import json,sys,time,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME='model:local-ta98:assembled-metres'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def ray_triangles(triangles,origin,direction,minimum=-.030,maximum=.040):
    # General Moller--Trumbore line/triangle intersection, keeping both
    # positive and negative signed parameters and actual barycentric data.
    t=triangles;e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0];h=np.cross(np.broadcast_to(direction,e2.shape),e2);det=np.einsum('ij,ij->i',e1,h)
    valid=abs(det)>1e-16;s=origin-t[:,0];u=np.divide(np.einsum('ij,ij->i',s,h),det,out=np.full(len(t),np.nan),where=valid);q=np.cross(s,e1)
    v=np.divide(q@direction,det,out=np.full(len(t),np.nan),where=valid);mu=np.divide(np.einsum('ij,ij->i',e2,q),det,out=np.full(len(t),np.nan),where=valid)
    keep=valid&(u>=0)&(v>=0)&(u+v<=1)&(mu>=minimum)&(mu<=maximum);ids=np.flatnonzero(keep)
    return [dict(face=int(i),distance_m=float(mu[i]),point_m=(origin+mu[i]*direction).tolist(),barycentric=[float(1-u[i]-v[i]),float(u[i]),float(v[i])])for i in ids]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--side',choices=['l','r'],required=True);ap.add_argument('--group',choices=['dorsal','marginal','deep'],required=True);ap.add_argument('--mode',choices=['rays','contacts'],required=True);args=ap.parse_args()
    OUT=BASE/('whole-foot-compartment-v29-'+args.side+'-'+args.group+'-'+args.mode)
    start=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    a=json.loads(CAT.read_bytes());ps={p['id']:p for p in a['parts']};reader=AtlasReader();cache={};bindings={}
    contract=json.loads((BASE/'course-contract-v23/contract.json').read_bytes());base=json.loads((BASE/'whole-foot-native-interfaces-v26/evidence.json').read_bytes())
    def get(pid):
        if pid not in cache:
            m=reader.read_part(a,ps[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED);cache[pid]=m
            bindings[pid]=dict(recordSha256=identity(ps[pid]),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'])
        return cache[pid]
    sample_rows=[];contacts=[]
    # Execute existing body-loader retirement rules by using its previously
    # verified default-hidden IDs. No fascia/skin/nerve is discarded by name.
    runtimePath=BASE/'whole-foot-runtime-qualification-v28.json';runtime=json.loads(runtimePath.read_bytes());assert runtime['sourceCatalogueSha256']==EXPECTED;known={r['id']:r['retained']for r in runtime['rows']}
    hidden=set(a.get('defaultHidden',[]))
    def retained(p):
        if p['id']in known:return known[p['id']]
        return not p.get('suppressed')and p['id']not in hidden
    for side,sign in [(args.side,1 if args.side=='l'else-1)]:
        skins=[p for p in a['parts']if p['system']=='integumentary'and retained(p)and p.get('bounds')and p['bounds'][0][1]<.14 and p['bounds'][1][1]>.003 and
            p['bounds'][0][0]<( .18 if sign==1 else-.05)and p['bounds'][1][0]>(.05 if sign==1 else-.18)]
        nerves=[p for p in a['parts']if p['system']=='nervous'and retained(p)and p.get('bounds')and p['bounds'][0][1]<.14 and p['bounds'][1][1]>.003 and
            p['bounds'][0][0]<( .18 if sign==1 else-.05)and p['bounds'][1][0]>(.05 if sign==1 else-.18)]
        fascia=get('ZA:Crural fascia.'+side)
        skintris=np.concatenate([get(p['id']).triangles for p in skins]);skinids=[p['id']for p in skins];skinoff=np.r_[0,np.cumsum([len(get(pid).faces)for pid in skinids])]
        nervetris=np.concatenate([get(p['id']).triangles for p in nerves]);nerveids=[p['id']for p in nerves];nerveoff=np.r_[0,np.cumsum([len(get(pid).faces)for pid in nerveids])]
        for prefix,kind,direction in [('ZA:Dorsal venous arch of foot.','dorsal',[0,1,0]),('ZA:Dorsal metatarsal veins.','dorsal',[0,1,0]),
                ('ZA:Dorsal digital veins of foot.','dorsal',[0,1,0]),('ZA:Small saphenous vein.','lateral',[sign,0,0]),('ZA:Great saphenous vein.','medial',[-sign,0,0]),
                ('ZA:Intercapitular veins of foot.','cross-compartment',[0,1,0]),('ZA:Plantar venous arch.','plantar',[0,-1,0]),
                ('ZA:Plantar metatarsal veins.','plantar',[0,-1,0]),('ZA:Plantar digital veins.','plantar',[0,-1,0]),
                ('ZA:Lateral plantar veins.','plantar',[0,-1,0]),('ZA:Medial plantar veins.','plantar',[0,-1,0])]:
            group='dorsal'if kind=='dorsal'else('marginal'if kind in ['lateral','medial']else'deep')
            if group!=args.group:continue
            atomic_json(OUT/('progress-'+identity(prefix+side)[:12]+'.json'),dict(sourcePartId=prefix+side,mode=args.mode,startedElapsedSeconds=time.monotonic()-start))
            pid=prefix+side;m=get(pid);v,iv=np.unique(m.vertices,axis=0,return_inverse=True);f=iv[m.faces];components=recover_station_rings(v,f)
            # Every source-connected cylinder contributes actual source
            # station samples; the original material/group name is not a
            # fitted centreline or a container that closes the network.
            selected_faces=[]
            for componentOrdinal,component in enumerate(components):
                rings=component['rings'];cc=np.array([v[q].mean(0)for q in rings]);foot=np.flatnonzero(cc[:,1]<.13)
                selected_faces.extend(component['face_ids'][np.all(cc[component['station_of_vertex'][f[component['face_ids']]],1]<.13,axis=1)].tolist())
                if not len(foot)or args.mode!='rays':continue
                ordinal=np.unique(np.r_[foot[0],foot[::12],foot[-1]])
                for j in ordinal:
                    origin=cc[j];ring=v[rings[j]];radius=float(np.linalg.norm(ring-origin,axis=1).max());n=np.array(direction,float)
                    fh=ray_triangles(fascia.triangles,origin,n);sh=ray_triangles(skintris,origin,n);nh=ray_triangles(nervetris,origin,n)
                    for hits,ids,offsets in [(sh,skinids,skinoff),(nh,nerveids,nerveoff)]:
                        for hit in hits:
                            k=int(np.searchsorted(offsets,hit['face'],side='right')-1);hit['partId']=ids[k];hit['face']-=int(offsets[k])
                    sh.sort(key=lambda z:z['distance_m']);fh.sort(key=lambda z:z['distance_m']);nh.sort(key=lambda z:z['distance_m'])
                    positive=[h for h in sh if h['distance_m']>0];outer=positive[0]if positive else None
                    below=[h for h in fh if outer and h['distance_m']<outer['distance_m']];deep=below[-1]if below else None
                    interval=[deep['distance_m']+radius,outer['distance_m']-radius]if deep and outer else None
                    sample_rows.append(dict(side=side,sourcePartId=pid,kind=kind,sourceComponent=componentOrdinal,sourceStation=int(j),
                        point_m=origin,direction=n,actualSourceRingCircumsphereRadius_m=radius,fasciaHits=fh,skinHits=sh,nerveHits=nh,
                        selectedFasciaWitness=deep,selectedSkinWitness=outer,coarseCentreTranslationIntervalAlongRay_m=interval,
                        positiveWidthIfBothHits=bool(interval and interval[1]>interval[0]),
                        currentSectionFitsBothRayBoundaries=bool(interval and interval[0]<=0<=interval[1]),
                        currentSkinIsOpenSurfaceNotVolumeCertificate=True,coarseRayIntervalIsNotAClearanceProof=True))
            # Actual complete source foot subsets versus every current
            # neighboring nerve, including unchanged data on plantar sides.
            if args.mode!='contacts':continue
            selected_faces=np.unique(selected_faces);region=Mesh(m.vertices,m.faces[selected_faces],FRAME);bb=region.bounds;tree=BVH(region)
            for p in nerves:
                nm=get(p['id']);nb=nm.bounds
                if np.any(nb[1]<bb[0])or np.any(nb[0]>bb[1]):continue
                detail,gate=pair_contacts(region,nm,0.,first_bvh=tree,max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget']
                contacts.append(dict(side=side,veinSourcePartId=pid,nerveSourcePartId=p['id'],gate=gate,
                    veinOriginalSourceFaceIds=selected_faces[detail['pairs'][:,0]],nerveOriginalSourceFaceIds=detail['pairs'][:,1],
                    closestSurfaceNotRecomputed='Actual surface-crossing gate and source face witnesses only; ray stage retains directional nerve distances.'))
    assert digest_file(CAT)==EXPECTED
    result=dict(schema='whole-foot-native-compartment-coarse-witnesses/v29',sourceCatalogueSha256=EXPECTED,stage=dict(side=args.side,group=args.group,mode=args.mode),
        actualRuntimeQualificationSha256=digest_file(runtimePath),sourceBindings=bindings,sourceInterfaceInventorySha256=digest_file(BASE/'whole-foot-native-interfaces-v26/evidence.json'),
        samples=sample_rows,actualNeighborNerveContacts=contacts,
        noRigidTranslationRotationOrFieldSelected=True,noGeometryConstructedOrModified=True,
        scope='Actual native outer skin and one-sided fascia ray intersections at literal source stations. Plantar and intercapitular source vessels are kept anatomically distinct from the superficial dorsal system.',
        limitations=['Ray samples and an open source skin patch do not establish an enclosing body volume or a continuous vascular corridor.',
            'A free ray interval can be obstructed by actual nerves, tendons, arteries, ligaments or neighboring veins. All must be checked on a proposed whole course.',
            'Native current source boundary contacts are not literal annuli or finite blood-flow interfaces. No topology is silently closed.',
            'For nearby foot records outside the previously executed loader context, literal source suppressed/defaultHidden predicates are recorded; no additional source preference is inferred.'],
        anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),samples=len(sample_rows),nervePairs=len(contacts),sourceBindings=len(bindings),elapsedSeconds=result['elapsedSeconds'],
        bilateralNerveCrossings=[dict(side=q['side'],source=q['veinSourcePartId'],nerve=q['nerveSourcePartId'],crossings=q['gate']['crossing_face_pairs'])for q in contacts if q['gate']['crossing_face_pairs']])))

if __name__=='__main__':main()
