"""Complete bounded major-route/bone triangle contacts before/after one fit."""
from pathlib import Path
from collections import Counter,defaultdict
import hashlib,json,sys,time
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
start=time.monotonic();B=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1';O=B/'cranial-route-fit-v1/actual-bone-check-v1'
if O.exists():raise SystemExit('Refusing to overwrite actual-bone-check-v1')
O.mkdir();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();fixed=B/'regional-interface-v1/fixed-context';bone=Mesh.load(fixed/'full.npz');a=np.load(fixed/'full.npz');part_index=a['face_owner_part_index'];source_face=a['face_owner_source_face_id'];owners=json.loads((fixed/'context.json').read_bytes())['partOwnership'];take=np.array([i for i,row in enumerate(owners)if row['region']=='headneck'or'Vertebra C7'in row['id']],int);face_ids=np.flatnonzero(np.isin(part_index,take));bone=Mesh(bone.vertices,bone.faces[face_ids],'active-male-metres');tree=BVH(bone);bone_keys=[owners[int(part_index[i])]['id']for i in face_ids]
names=['Left common carotid artery','Right common carotid artery','Left internal carotid artery','Right internal carotid artery','Left vertebral artery','Right vertebral artery','Basilar artery'];results={};overall=[]
for stage,directory in [('before_v3',B/'regional-fit-v3/candidate-v1'),('after_cranial_v1',B/'cranial-route-fit-v1/candidate-v1')]:
 mesh=Mesh.load(directory/'arterial.npz');rows=json.loads((directory/'instance-ownership.json').read_bytes())['instances'];values=[];faces=[];face_owners=[];source_candidate_faces=[];offset=0
 for name in names:
  row=next(r for r in rows if r['name']==name);vs=row['candidateVertexStart'];vc=row['candidateVertexCount'];fs=row['candidateFaceStart'];fc=row['candidateFaceCount'];values.append(mesh.vertices[vs:vs+vc]);faces.append(mesh.faces[fs:fs+fc]-vs+offset);face_owners.extend([name]*fc);source_candidate_faces.extend(range(fs,fs+fc));offset+=vc
 major=Mesh(np.concatenate(values),np.concatenate(faces),'active-male-metres');contacts,summary=pair_contacts(major,bone,0.,second_bvh=tree,max_pairs=100000,max_tests=600000);strict=np.flatnonzero(contacts['relation']=='crossing');pairs=contacts['pairs'][strict];by_vessel=Counter(face_owners[int(p[0])]for p in pairs);by_pair=Counter((face_owners[int(p[0])],bone_keys[int(p[1])])for p in pairs);counts=[dict(vessel=k[0],bonePart=k[1],strictTrianglePairWitnesses=n)for k,n in by_pair.most_common()]
 contacts.update(first_candidate_arterial_face_ids=np.array(source_candidate_faces,int),second_full_context_face_ids=face_ids,second_original_bone_part_face_ids=source_face[face_ids]);np.savez_compressed(O/(stage+'-contacts.npz'),**contacts)
 results[stage]=dict(summary=summary,byVessel={name:int(by_vessel[name])for name in names},byVesselAndBone=counts,firstMeshVertices=len(major.vertices),firstMeshFaces=len(major.faces),fixedBoneFaces=len(bone.faces),surfaceContactNotTissueOccupancy=True)
 overall.append(dict(stage=stage,strictTrianglePairs=summary['crossing_face_pairs'],complete=summary['complete_within_pair_budget'],byVessel=results[stage]['byVessel']))
 (O/'checkpoint.json').write_text(json.dumps(dict(stage=stage,complete=summary['complete_within_pair_budget'],elapsedSeconds=time.monotonic()-start),indent=2)+'\n')
report=dict(schema='actual-major-cranial-route-fixed-bone-check/v1',results=results,comparableCompleteTrianglePairQuery=all(r['summary']['complete_within_pair_budget']for r in results.values()),fixedSkeletonSha256=sha(fixed/'full.npz'),sourceBeforeOwnershipSha256=sha(B/'regional-fit-v3/candidate-v1/instance-ownership.json'),candidateAfterOwnershipSha256=sha(B/'cranial-route-fit-v1/candidate-v1/instance-ownership.json'),strictWitnessesAreActualTrianglePairs=True,notComparableTo32FaceSampleCounts=True,containmentAndPatencyUntested=True,allSevenNamedRoutesIncluded=True,nutrientEmissaryAndOtherBranchesNotAutomaticallyRejected=True,noGeometryAlteredByThisCheck=True,anatomicalAcceptance=False,installed=False,elapsedSeconds=time.monotonic()-start)
(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(results=overall,elapsedSeconds=report['elapsedSeconds'])))
