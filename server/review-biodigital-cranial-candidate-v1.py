"""Saved actual triangle-area and projection review, without spatial queries."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
from PIL import Image,ImageDraw
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from geometry_tools.model import Mesh
start=time.monotonic();B=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1';O=B/'cranial-route-fit-v1';out=O/'review-v1';out.mkdir(exist_ok=False)
def area(t):return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
def stats(x):return dict(count=len(x),minimum=float(np.min(x)),maximum=float(np.max(x)),p01=float(np.quantile(x,.01)),median=float(np.median(x)),p99=float(np.quantile(x,.99)),belowHalf=int(np.sum(x<.5)),outsideTenPercent=int(np.sum((x<.9)|(x>1.1))))
ownership=json.loads((O/'candidate-v1/instance-ownership.json').read_bytes())['instances'];combined=[];additional=[];selected=[];critical=['Left common carotid artery','Right common carotid artery','Left internal carotid artery','Right internal carotid artery','Left vertebral artery','Right vertebral artery','Basilar artery','Brachiocephalic trunk','Left brachiocephalic vein','Right brachiocephalic vein','Left subclavian artery','Right subclavian artery','Left subclavian vein','Right subclavian vein']
for system in ['arterial','venous','arterial-remnant-companion']:
 old=Mesh.load(B/'regional-fit-v3/candidate-v1'/(system+'.npz'));new=Mesh.load(O/'candidate-v1'/(system+'.npz'));assert np.array_equal(old.faces,new.faces);oa=area(old.triangles);na=area(new.triangles);ratio=na/oa;prior=np.load(B/'regional-fit-v3/skeleton-samples-v1'/(system+'-triangle-area-ratios.npz'))['candidate_face_area_ratio'];total=ratio*prior;assert np.isfinite(total).all();combined.append(total);additional.append(ratio)
 for row in ownership:
  if row['candidateMesh']!=system+'.npz'or row['name']not in critical:continue
  sl=slice(row['candidateFaceStart'],row['candidateFaceStart']+row['candidateFaceCount']);selected.append(dict(name=row['name'],sourceInstanceId=row['sourceInstanceId'],additionalAreaRatioVsV3=stats(ratio[sl]),totalAreaRatioVsG=stats(total[sl])))
report=dict(schema='cranial-candidate-actual-triangle-area-review/v1',allRetainedFaces=sum(len(x)for x in combined),additionalAreaRatioVsV3=stats(np.concatenate(additional)),totalAreaRatioVsG=stats(np.concatenate(combined)),criticalInstances=selected,qualification='Actual stored triangle surface areas and identical source faces. These are not lumen calibre, tissue volume or clinical tolerances.',installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start);(out/'triangle-area-review.json').write_text(json.dumps(report,indent=2)+'\n')
# Project the exact basilar/paired bone triangles that establish the new defect.
contacts=np.load(O/'actual-bone-check-v1/after_cranial_v1-contacts.npz');pairs=contacts['pairs'];strict=contacts['relation']=='crossing';source_faces=contacts['first_candidate_arterial_face_ids'];bone_faces=contacts['second_full_context_face_ids'];row=next(r for r in ownership if r['name']=='Basilar artery');source_global=source_faces[pairs[:,0]];keep=strict&(source_global>=row['candidateFaceStart'])&(source_global<row['candidateFaceStart']+row['candidateFaceCount']);cross_bone_ids=np.unique(bone_faces[pairs[keep,1]]);fixed=Mesh.load(B/'regional-interface-v1/fixed-context/full.npz');bone_tri=fixed.triangles[cross_bone_ids];before=Mesh.load(B/'regional-fit-v3/candidate-v1/arterial.npz').triangles[row['candidateFaceStart']:row['candidateFaceStart']+row['candidateFaceCount']];after=Mesh.load(O/'candidate-v1/arterial.npz').triangles[row['candidateFaceStart']:row['candidateFaceStart']+row['candidateFaceCount']];center=np.concatenate([before.reshape(-1,3),after.reshape(-1,3)]).mean(0);image=Image.new('RGB',(1200,850),(14,20,28));draw=ImageDraw.Draw(image);scale=11000
for col,basis in enumerate([np.array([[1.,0,0],[0,1.,0],[0,0,1.]]),np.array([[0.,0,1],[0,1.,0],[-1,0,0]])]):
 for r,(label,triangles)in enumerate([('Before V3',before),('After cranial candidate',after)]):
  offset=np.array([col*600+300,r*425+240]);objects=[]
  for tri in bone_tri:
   projected=(tri-center)@basis.T;objects.append((float(projected[:,2].mean()),projected,(98,108,124),(130,140,155)))
  for tri in triangles:
   projected=(tri-center)@basis.T;objects.append((float(projected[:,2].mean()),projected,(190,66,48),(255,145,110)))
  for _,p,color,edge in sorted(objects,key=lambda x:x[0]):
   q=p[:,:2]*scale;q[:,1]*=-1;q+=offset;draw.polygon([tuple(x)for x in q],fill=color,outline=edge)
  draw.text((col*600+12,r*425+12),label+(' / frontal'if col==0 else' / sagittal'),fill='white');draw.text((col*600+12,r*425+32),'Exact basilar triangles and bone faces in 56 strict after-pairs',fill='white')
image.save(out/'basilar-regression-projections.png')
history=dict(preparationFailure=dict(runtimeReceipt='runtime/a671e9a61eac4541843217e32b36beea.json',elapsedSeconds=1.02,reason='Dilation/twist identity outside support recomputed c + (p - c), changing last bits of fixed controls.'),implementationFix='Literal out=p and J=I outside support for every compact primitive.',focusedTests=dict(runtimeReceipt='runtime/1ab5e0abd27f4d689cb5223fcf94f1ed.json',tests=6,passed=True,elapsedSeconds=.22),preparationSuccess=dict(runtimeReceipt='runtime/b19119045b04445e981ecdd3b2ba2af2.json',elapsedSeconds=.71),fullCandidate=dict(runtimeReceipt='runtime/93c37fcb56644721801a98afa41f32c1.json',elapsedSeconds=15.61),actualMajorBoneCheck=dict(runtimeReceipt='runtime/5e6afe30b46349cda08aac8a0df6dba3.json',elapsedSeconds=2.20),currentNativeGuideCheck=dict(runtimeReceipt='runtime/44296a828abb45e69cc9a786a995d548.json',elapsedSeconds=2.04),oneAnatomicalCandidateOnly=True,sourceBanksAndTargetUnchanged=True,newBasilarConflictMeansCandidateUnaccepted=True);(out/'stage-history.json').write_text(json.dumps(history,indent=2)+'\n');print(json.dumps(dict(area=report['totalAreaRatioVsG'],additional=report['additionalAreaRatioVsV3'],elapsedSeconds=time.monotonic()-start,image=str(out/'basilar-regression-projections.png'))))
