"""Bounded sampled crossing and actual triangle-area diagnostics; no acceptance.

Compressed arrays used in loops are cached once. Saved reference witnesses are
replayed without another spatial query when --reference-samples is supplied.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'scripts')]
from biodigital_blood_source import BloodSource
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.distance import triangle_pairs

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--candidate-stage',default='regional-fit-v1')
parser.add_argument('--output',default='skeleton-samples-v1',help='New directory relative to candidate stage, or absolute path')
parser.add_argument('--reference-samples',help='Existing witness directory relative to workbench, or absolute path')
parser.add_argument('--additional-reference-samples',action='append',default=[],help='Other saved candidate witnesses to compare without another query')
parser.add_argument('--critical-strain-reference-stage',action='append',default=[],help='Also measure brachiocephalic/subclavian actual area strain from this saved candidate stage')
parser.add_argument('--area-strain',action='store_true')
args=parser.parse_args();started=time.monotonic();base=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';stage=base/args.candidate_stage;candidate=stage/'candidate-v1';output=stage/args.output
output.mkdir(exist_ok=False);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
reference=(base/args.reference_samples)if args.reference_samples else None
context=base/'regional-interface-v1/fixed-context';source=BloodSource(base/'source-v2');common=np.array(json.loads((base/'registration-v2/transform.json').read_bytes())['matrix'])
bone=Mesh.load(context/'full.npz');fixed=np.load(context/'full.npz');bone_owner=fixed['face_owner_part_index'];bone_source_face=fixed['face_owner_source_face_id'];owners=json.loads((context/'context.json').read_bytes())['partOwnership'];bone_triangles=bone.triangles
part_rows=json.loads((candidate/'instance-ownership.json').read_bytes())['instances'];reference_rows={r['sourceInstanceId']:r for r in json.loads((reference.parent/'candidate-v1/instance-ownership.json').read_bytes())['instances']}if reference else{}
used=np.unique(bone.faces);center=(bone.vertices[used].min(0)+bone.vertices[used].max(0))/2;extent=np.ptp(bone.vertices[used],axis=0).max();native=((bone.vertices-center)/extent).astype(np.float32).astype(float);nt=native[bone.faces];mask=np.linalg.norm(np.cross(nt[:,1]-nt[:,0],nt[:,2]-nt[:,0]),axis=1)>0;face_map=np.flatnonzero(mask);index=SurfaceIndex(Mesh(bone.vertices,bone.faces[mask],bone.frame));del native,nt
np.savez_compressed(output/'proposal-face-ownership.npz',original_bone_face_ids=face_map,excluded_from_native_proposal_only=np.flatnonzero(~mask))
checkpoint=output/'checkpoint.json'
def save_checkpoint(status,**extra):checkpoint.write_text(json.dumps(dict(status=status,candidateStage=args.candidate_stage,originalBoneFaces=len(bone.faces),nativeProposalFaces=len(face_map),sourceGeometryUnchanged=True,elapsedSeconds=time.monotonic()-started,**extra),indent=2)+'\n')
save_checkpoint('native-index-ready')

def area2(tri):return np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
def distribution(values):
 values=np.asarray(values);good=values[np.isfinite(values)];return dict(count=int(len(values)),finiteCount=int(len(good)),minimum=float(good.min())if len(good)else None,maximum=float(good.max())if len(good)else None,quantiles={str(p):float(v)for p,v in zip([.001,.01,.05,.5,.95,.99,.999],np.quantile(good,[.001,.01,.05,.5,.95,.99,.999]))}if len(good)else{},belowHalf=int(np.sum(good<.5)),aboveDouble=int(np.sum(good>2)),outsideTenPercent=int(np.sum((good<.9)|(good>1.1))))
sample_meta=[];before=[];after=[];area_reports=[];all_ratios=[];all_instance_area=[];critical_reference_area=[]
critical_reference_rows={s:{r['sourceInstanceId']:r for r in json.loads((base/s/'candidate-v1/instance-ownership.json').read_bytes())['instances']}for s in args.critical_strain_reference_stage}
for system in ['arterial','venous','arterial-remnant-companion']:
 mesh=Mesh.load(candidate/(system+'.npz'));triangles=mesh.triangles;az=np.load(candidate/(system+'-ancestry.npz'));source_face=az['candidate_face_source_face'];source_instance=az['candidate_face_source_instance'];ratios=np.full(len(mesh.faces),np.nan);ref_source_face=None
 critical_meshes={s:Mesh.load(base/s/'candidate-v1'/(system+'.npz'))for s in args.critical_strain_reference_stage};critical_triangles={s:m.triangles for s,m in critical_meshes.items()};critical_source_faces={s:np.load(base/s/'candidate-v1'/(system+'-ancestry.npz'))['candidate_face_source_face']for s in args.critical_strain_reference_stage}
 if reference:ref_source_face=np.load(reference.parent/'candidate-v1'/(system+'-ancestry.npz'))['candidate_face_source_face']
 for row in part_rows:
  if row['candidateMesh']!=system+'.npz':continue
  start=row['candidateFaceStart'];count=row['candidateFaceCount'];ids=start+np.unique(np.linspace(0,count-1,min(32,count)).astype(int));sf=source_face[ids];instance=source.by_id[row['sourceInstanceId']]
  if reference:
   rr=reference_rows[row['sourceInstanceId']];ri=rr['candidateFaceStart']+np.unique(np.linspace(0,rr['candidateFaceCount']-1,min(32,rr['candidateFaceCount'])).astype(int));assert np.array_equal(sf,ref_source_face[ri]);assert count==rr['candidateFaceCount']
  if not reference:
   original=source.world(instance);p=original['positions']@common[:3,:3].T+common[:3,3];before.extend(p[original['indices'][sf]])
  after.extend(triangles[ids]);sample_meta.extend(dict(instance=row['sourceInstanceId'],name=row['name'],candidate_mesh=system,candidate_face=int(i),source_face=int(s))for i,s in zip(ids,sf))
  if args.area_strain:
   geometry=source.geometry(instance['geometryId']);world=np.array(instance['worldMatrix']).reshape(4,4,order='F');affine=common[:3,:3]@world[:3,:3];translation=common[:3,:3]@world[:3,3]+common[:3,3];registered=geometry['positions'].astype(float)@affine.T+translation;ff=geometry['indices'][source_face[start:start+count]]
   if instance['reverseWinding']:ff=ff[:,::-1]
   oldarea=area2(registered[ff]);newarea=area2(triangles[start:start+count]);valid=(oldarea>0)&np.isfinite(oldarea)&np.isfinite(newarea);local=np.divide(newarea,oldarea,out=np.full(count,np.nan),where=valid);ratios[start:start+count]=local;all_instance_area.append(dict(instance=row['sourceInstanceId'],name=row['name'],candidateMesh=system,candidateFaceStart=start,faceCount=count,baselineZeroOrInvalidFaces=int(np.sum(~valid)),areaRatio=distribution(local)))
   if any(term in row['name'].lower()for term in ['brachiocephalic','subclavian']):
    for s in args.critical_strain_reference_stage:
     cr=critical_reference_rows[s][row['sourceInstanceId']];cs=cr['candidateFaceStart'];assert cr['candidateFaceCount']==count;assert np.array_equal(critical_source_faces[s][cs:cs+count],source_face[start:start+count]);ca=area2(critical_triangles[s][cs:cs+count]);ratio=np.divide(ca,oldarea,out=np.full(count,np.nan),where=valid);critical_reference_area.append(dict(candidateStage=s,instance=row['sourceInstanceId'],name=row['name'],faceCount=count,areaRatio=distribution(ratio)))
 if args.area_strain:
  np.savez_compressed(output/(system+'-triangle-area-ratios.npz'),candidate_face_area_ratio=ratios,candidate_face_source_face=source_face,candidate_face_source_instance=source_instance)
  area_reports.append(dict(candidateMesh=system,areaRatio=distribution(ratios)));all_ratios.append(ratios)
 del mesh,triangles,az,critical_meshes,critical_triangles,critical_source_faces
after=np.array(after);before=np.array(before);assert len(sample_meta)==55108;assert len({m['instance']for m in sample_meta})==1724
np.savez_compressed(output/'sample-ownership.npz',source_instance_ids=np.array([m['instance']for m in sample_meta]),candidate_mesh=np.array([m['candidate_mesh']for m in sample_meta]),candidate_face=np.array([m['candidate_face']for m in sample_meta]),source_face=np.array([m['source_face']for m in sample_meta]))
save_checkpoint('samples-and-area-strain-ready',sampledBloodFaces=len(sample_meta),areaStrainComplete=args.area_strain)
results={};stage_arrays={}
if reference:
 ref_report=json.loads((reference/'report.json').read_bytes())
 for label in ['initial_similarity','regional_candidate']:
  rz=np.load(reference/(label+'-witnesses.npz'));cached={k:rz[k]for k in rz.files};assert len(cached['relation'])==len(sample_meta);stage_arrays['initial_similarity'if label=='initial_similarity'else'previous_regional_candidate']=cached
 results['initial_similarity']={**ref_report['results']['initial_similarity'],'reusedSavedWitnesses':True,'sourceWitnessDirectory':str(reference)}
 np.savez_compressed(output/'initial_similarity-witnesses.npz',**stage_arrays['initial_similarity'])
for additional in args.additional_reference_samples:
 directory=base/additional;rz=np.load(directory/'regional_candidate-witnesses.npz');cached={k:rz[k]for k in rz.files};assert len(cached['relation'])==len(sample_meta);stage_arrays['saved_candidate:'+directory.parent.name]=cached
pending=[('regional_candidate',after)]if reference else[('initial_similarity',before),('regional_candidate',after)]
for label,tri in pending:
 closest,evidence=index.closest_points(tri.mean(axis=1));bf=face_map[closest['face_ids']];actual=triangle_pairs(tri,bone_triangles[bf]);relation=actual['relation'];crossing=np.flatnonzero(relation=='crossing');examples=[]
 for i in crossing[:80]:
  face=int(bf[i]);pi=int(bone_owner[face]);examples.append(dict(**sample_meta[i],bone_part=owners[pi]['id'],compound_bone_face=face,original_bone_part_face=int(bone_source_face[face]),strict_relation=str(relation[i]),vessel_triangle_m=tri[i].tolist(),bone_triangle_m=bone_triangles[face].tolist()))
 arrays=dict(sampled_triangles=tri,bone_face_ids=bf,relation=relation,crossing_sample_indices=crossing);np.savez_compressed(output/(label+'-witnesses.npz'),**arrays);stage_arrays[label]=arrays
 results[label]=dict(sampledBloodFaces=len(tri),strictCrossingWitnessCount=len(crossing),distinctVesselInstancesWithCrossingWitnesses=len({sample_meta[i]['instance']for i in crossing}),examples=examples,nativeProposalEvidence=evidence)
save_checkpoint('strict-witness-query-complete',sampledBloodFaces=len(sample_meta))

def counted(rows,key):return[{key:k,'sampledStrictWitnesses':n}for k,n in Counter(r[key]for r in rows).most_common()]
comparisons={};witness_rows={}
for label,arrays in stage_arrays.items():
 rows=[];bf=arrays['bone_face_ids'];centroids=arrays['sampled_triangles'].mean(1)
 for i in arrays['crossing_sample_indices']:
  face=int(bf[i]);owner=owners[int(bone_owner[face])];rows.append(dict(**sample_meta[i],sampleIndex=int(i),anatomicalRegion=owner['region'],bonePart=owner['id'],originalBonePartFace=int(bone_source_face[face]),centroid_m=centroids[i].tolist()))
 witness_rows[label]=rows;comparisons[label]=dict(sampledStrictWitnesses=len(rows),byRegion=counted(rows,'anatomicalRegion'),byBonePart=counted(rows,'bonePart'),byInstance=[dict(instance=k,name=next(r['name']for r in rows if r['instance']==k),sampledStrictWitnesses=n)for k,n in Counter(r['instance']for r in rows).most_common()])
pair_counts=defaultdict(Counter)
for label,rows in witness_rows.items():
 for r in rows:pair_counts[(r['instance'],r['name'],r['anatomicalRegion'],r['bonePart'])][label]+=1
pair_rows=[dict(instance=k[0],name=k[1],anatomicalRegion=k[2],bonePart=k[3],**{stage:int(v.get(stage,0))for stage in stage_arrays},deltaVsPrevious=int(v.get('regional_candidate',0)-v.get('previous_regional_candidate',0)),deltaVsG=int(v.get('regional_candidate',0)-v.get('initial_similarity',0)))for k,v in pair_counts.items()];pair_rows.sort(key=lambda r:(-r['regional_candidate'],r['name']))
comparison=dict(stages=comparisons,pairs=pair_rows,sameDeterministicSampleSourceFaceIds=True,sampledBloodFaces=len(sample_meta),completeCrossingInventory=False,noAnatomicalAcceptance=True)
(output/'comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
if args.area_strain:
 strain=dict(schema='actual-triangle-area-ratios/v1',baseline='Immutable G applied to source-v2 controls with the same source faces and placement arithmetic as candidate builder',candidateStage=args.candidate_stage,allRetainedTriangleCount=sum(len(x)for x in all_ratios),areaRatio=distribution(np.concatenate(all_ratios)),systems=area_reports,instances=all_instance_area,qualification='Actual triangle area after/before, dimensionless; surface triangle strain, not lumen calibre, wall thickness, vessel diameter, volume or anatomical acceptance.')
 (output/'triangle-area-strain.json').write_text(json.dumps(strain,indent=2)+'\n')
 if critical_reference_area:(output/'critical-reference-area-strain.json').write_text(json.dumps(dict(referenceInstances=critical_reference_area,qualification=strain['qualification']),indent=2)+'\n')
report=dict(schema='biodigital-fixed-skeleton-sampled-crossings/v2',candidateStage=args.candidate_stage,results=results,scope='32 deterministic retained source faces per vessel instance; one nearest native-proposed bone face per sample, replayed with original Float64 triangles',allVesselInstancesRepresented=1724,originalFixedBoneFaces=len(bone.faces),omittedNativeProposalFaces=int((~mask).sum()),positiveCrossingWitnessesAreReal=True,zeroCountIsNotClearanceProof=True,completeCrossingInventory=False,containmentTested=False,anatomicalAcceptance=False,installed=False,sourceModified=False,referenceSamples=str(reference)if reference else None,sameSampleSourceFaceIdsVerified=bool(reference),builderSha256=sha(__file__),candidateOwnershipSha256=sha(candidate/'instance-ownership.json'),fixedContextSha256=sha(context/'full.npz'),elapsedSeconds=time.monotonic()-started)
(output/'report.json').write_text(json.dumps(report,indent=2)+'\n');save_checkpoint('complete-sampled-crossings-and-area-diagnostics')
print(json.dumps(dict(path=str(output/'report.json'),stages={k:{'strictWitnesses':v['sampledStrictWitnesses'],'regions':v['byRegion']}for k,v in comparisons.items()},areaRatio=strain['areaRatio']if args.area_strain else None,elapsedSeconds=report['elapsedSeconds'])))
