"""Read-only literal component/rib/recipient ownership for two named casts."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/posterior-intercostal-native-ownership-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
ORDINAL=['First','Second','Third','Fourth','Fifth','Sixth','Seventh','Eighth','Ninth','Tenth','Eleventh','Twelfth']

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def bounds_hit(a,b):return bool(np.all(a[1]>=b[0])and np.all(b[1]>=a[0]))

def pair(a,b,tree=None):
    if not bounds_hit(a.bounds,b.bounds):return dict(crossing_face_pairs=0,coplanar_overlap_face_pairs=0,complete_within_pair_budget=True,broadPhaseDisjoint=True),[]
    details,gate=pair_contacts(a,b,0.,first_bvh=tree,max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget'];bad=np.isin(details['relation'],['crossing','coplanar-overlap']);return gate,details['pairs'][bad]

def main():
    args=argparse.ArgumentParser();args.add_argument('--side',choices=['l','r'],required=True);side=args.parse_args().side
    start=time.monotonic();OUT=BASE/side;OUT.mkdir(parents=True,exist_ok=False);assert digest_file(CAT)==EXPECTED
    atlas=json.loads(CAT.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();cache={};bindings={}
    def get(pid):
        if pid not in cache:
            p=parts[pid];m=reader.read_part(atlas,p,frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED);cache[pid]=m;bindings[pid]=dict(partId=pid,sourceRecord=p,recordSha256=identity(p),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'],bankPath=m.source['bank_path'])
        return cache[pid]
    pid='ZA:'+('Left'if side=='l'else'Right')+' superior intercostal vein';source=get(pid);v,inv=np.unique(source.vertices,axis=0,return_inverse=True);f=inv[source.faces];courses=recover_station_rings(v,f);rows=[];centrepaths={};componentMeshes={}
    recipientIds=['ZA:Azygos vein','ZA:Hemi-azygos vein','ZA:Accessory hemi-azygos vein','ZA:Left brachiocephalic vein','ZA:Right brachiocephalic vein','ZA:Superior vena cava']
    recipients={q:(get(q),BVH(get(q)))for q in recipientIds}
    ribIds=[q for q in parts if q.startswith('ZA:')and q.endswith('.'+side+':Bone-'+q.rsplit('Bone-',1)[-1])and any(q.startswith('ZA:'+n+' rib.')for n in ORDINAL)]
    ribs={q:(get(q),BVH(get(q)))for q in ribIds}
    spineIds=['ZA:Vertebra T'+str(n)+':Bone-5'for n in range(1,13)];spine={q:get(q)for q in spineIds}
    for index,course in enumerate(courses):
        faces=np.array(course['face_ids']);rings=np.array(course['rings']);centers=v[rings].mean(1);centrepaths[index]=centers;component=Mesh(source.vertices,source.faces[faces],source.frame);componentMeshes[index]=component;tree=BVH(component);sampleIds=np.unique(np.rint(np.linspace(0,len(centers)-1,21)).astype(int));sample=centers[sampleIds];ribRows=[]
        for rid,(bone,btree)in ribs.items():
            closest,meta=btree.closest_points(sample,backend='bvh');gate,pairs=pair(component,bone,tree)
            ribRows.append(dict(partId=rid,medianCentreToSourceSurfaceDistance_m=float(np.median(closest['distance_m'])),minimumCentreToSourceSurfaceDistance_m=float(closest['distance_m'].min()),maximumCentreToSourceSurfaceDistance_m=float(closest['distance_m'].max()),sampleSourceStationIds=sampleIds,closestSourceFaceIds=closest['face_ids'],closestSourcePoints_m=closest['points'],verticalCentreMinusClosestBone_m=sample[:,1]-closest['points'][:,1],actualContactGate=gate,actualCrossingSourceFacePairs=[[int(faces[a]),int(b)]for a,b in pairs]))
        ribRows.sort(key=lambda q:q['medianCentreToSourceSurfaceDistance_m']);spineRows=[]
        for rid,bone in spine.items():
            if not bounds_hit(component.bounds,bone.bounds):continue
            gate,pairs=pair(component,bone,tree);spineRows.append(dict(partId=rid,gate=gate,actualCrossingSourceFacePairs=[[int(faces[a]),int(b)]for a,b in pairs]))
        mouths=[]
        for end,j in [('first',0),('last',len(rings)-1)]:
            ring=rings[j];rim=v[ring];distance=[];endpointCollarVertices=set(rings[:3].ravel()if j==0 else rings[-3:].ravel());ff=np.flatnonzero(np.isin(f[faces],list(endpointCollarVertices)).any(1));collar=Mesh(source.vertices,source.faces[faces[ff]],source.frame);ct=BVH(collar)
            for rid,(recipient,rtree)in recipients.items():
                closest,meta=rtree.closest_points(rim,backend='bvh');center=rtree.closest_point(centers[j]);gate,pairs=pair(collar,recipient,ct);literal=[int(np.any(np.all(recipient.vertices==q,axis=1)))for q in rim]
                distance.append(dict(recipientId=rid,centreToActualSourceSurface=center,minimumRimToSurfaceDistance_m=float(closest['distance_m'].min()),maximumRimToSurfaceDistance_m=float(closest['distance_m'].max()),literalSharedRimCoordinates=sum(literal),actualRimSourceFaceIds=closest['face_ids'],actualRimSourceBarycentric=closest['barycentric'],actualRimClosestPoints_m=closest['points'],endCollarFaceIds=faces[ff],endCollarContactGate=gate,endCollarCrossingSourceFacePairs=[[int(faces[ff[a]]),int(b)]for a,b in pairs],noMatchedRecipientAnnulusEstablished=True))
            distance.sort(key=lambda q:q['centreToActualSourceSurface']['distance_m']);mouths.append(dict(end=end,sourcePhysicalVertexIds=ring,center_m=centers[j],sourceRimCoordinates_m=rim,nearestActualRecipients=distance))
        nearestSpine=[]
        for mouth in mouths:
            p=np.array(mouth['center_m']);nearest=[(BVH(bone).closest_point(p),rid)for rid,bone in spine.items()];nearest.sort(key=lambda q:q[0]['distance_m']);nearestSpine.append(dict(end=mouth['end'],nearestActualSpinalSources=[dict(partId=rid,witness=q)for q,rid in nearest[:2]]))
        rows.append(dict(component=index,sourceComponentId=course['source_component'],sourceFaceIds=faces,faces=len(faces),stations=len(rings),bounds_m=component.bounds,sourceCentres_m=centers,rankedRibCompatibility=ribRows,actualSpinalContacts=spineRows,mouths=mouths,mouthSpinalLevelWitnesses=nearestSpine))
        atomic_json(OUT/('component-'+str(index)+'.json'),json.loads(json.dumps(rows[-1],default=clean)))
    mutual=[]
    for i,aa in centrepaths.items():
        for j,bb in centrepaths.items():
            if j<=i:continue
            d1=cKDTree(aa).query(bb)[0];d2=cKDTree(bb).query(aa)[0];haus=float(max(d1.max(),d2.max()))
            if haus<.003:
                gate,pairs=pair(componentMeshes[i],componentMeshes[j]);mutual.append(dict(first=i,second=j,symmetricCentrelineHausdorff_m=haus,actualContactGate=gate,notAutomaticallyDuplicateByCentrelineProximity=True))
    refs=[dict(pmid='37681447',url='https://pubmed.ncbi.nlm.nih.gov/37681447/',basis='66-cadaver superior-collector composition and destinations; collector is not all thoracic intercostal courses'),dict(pmid='23813753',url='https://pubmed.ncbi.nlm.nih.gov/23813753/',basis='21-cadaver actual posterior intercostal/azygos/intervertebral connections'),dict(pmid='23303031',url='https://pubmed.ncbi.nlm.nih.gov/23303031/',basis='48-cadaver azygos-system levels and variant connections')]
    result=dict(schema='literal-native-intercostal-component-and-recipient-ownership/v1',side=side,sourceCatalogueSha256=EXPECTED,term='A12.3.07.014',sourcePartId=pid,sourceBindings=bindings,components=len(rows),componentEvidenceFiles=[str(OUT/('component-'+str(i)+'.json'))for i in range(len(rows))],nearDuplicateCourseDiagnostics=mutual,primaryReferences=refs,
        wholeSourceMustNotBeAliasedToNewTerm=True,noPatencyFromCentresOrOverlapCasts=True,noSourceGeometryLabelsMetadataOrInstallationChanged=True,ownershipAndAnatomyNotAcceptedYet=True,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],components=len(rows),summary=[dict(component=q['component'],faces=q['faces'],rib=q['rankedRibCompatibility'][0]['partId'],medianDistance_mm=q['rankedRibCompatibility'][0]['medianCentreToSourceSurfaceDistance_m']*1000,allRibCrossingPairs=sum(r['actualContactGate']['crossing_face_pairs']for r in q['rankedRibCompatibility']),spineCrossingPairs=sum(r['gate']['crossing_face_pairs']for r in q['actualSpinalContacts']),recipient=[(m['nearestActualRecipients'][0]['recipientId'],m['nearestActualRecipients'][0]['centreToActualSourceSurface']['distance_m']*1000)for m in q['mouths']])for q in rows],nearDuplicate=mutual)))

if __name__=='__main__':main()
