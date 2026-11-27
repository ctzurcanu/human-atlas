"""Create one bounded, common positive thoracic correction contract.

Whole-bone surface centroids are coarse region guides, never foramina or
vessel-contact targets. No raw donor, global registration or skeleton changes.
"""
from pathlib import Path
import hashlib,importlib.util,json,sys,time
import numpy as np
from scipy.optimize import least_squares
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R),str(R/'scripts')]
from geometry_tools.model import AtlasReader
from geometry_tools.algorithms.compact_positive_fields import compose_compact_radial_maps
from biodigital_blood_source import BloodSource
start=time.monotonic();B=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1';O=B/'regional-fit-v2'
if O.exists():raise SystemExit('Refusing to overwrite regional-fit-v2')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('bone_landmarks',R/'server/biodigital-registration-landmarks.py');lm=importlib.util.module_from_spec(spec);spec.loader.exec_module(lm)
lp=B/'landmarks/landmarks.json';landmarks=json.loads(lp.read_bytes());Gpath=B/'registration-v2/transform.json';G=np.array(json.loads(Gpath.read_bytes())['matrix']);old=json.loads((B/'regional-fit-v1/contract.json').read_bytes());cp=Path(landmarks['target_catalog']);assert sha(cp)==old['target_catalogue_sha256'];target=json.loads(cp.read_bytes())
source_catalog=Path(landmarks['source_catalog']);src=json.loads(source_catalog.read_bytes());assert sha(source_catalog)==landmarks['source_catalog_sha256']
entries={n['objectId']:(n,m,anc)for n,m,anc in lm.flatten(src['sceneGraph'],src['records'])};reader=AtlasReader(root=R)
vispath=B.parent/'circulatory-completion/native-mca-controls-v3/proposal-v2/resolved-installed.json';vis={x['id']:x for x in json.loads(vispath.read_bytes())['effectivePartVisibility']};rows=[];arrays={}
for label in ['T1','T7','T12']:
 row=next(x for x in landmarks['landmarks']if x['label']==label);p=np.array(row['source']['area_weighted_surface_centroid']);rows.append(dict(label=label,source_registered_m=(p@G[:3,:3].T+G[:3,3]).tolist(),target_m=row['target']['area_weighted_surface_centroid'],qualification='Area-weighted complete named bone exterior surface region centroid. Not a vertebral body centre, foraminal point, vessel corridor or anatomical contact.',existingEvidenceLabel=label,existingEvidenceSha256=sha(lp)))
for side,suffix in [('left','l'),('right','r')]:
 sid='human_20_male_skeletal_system-'+side+'_first_rib_ID';node,matrix,anc=entries[sid];key=lm.resource_id(node['geometryId'],'geometry');ref=src['geometry'][key]
 with (source_catalog.parent/ref['file']).open('rb')as stream:stream.seek(ref['offset']);raw=stream.read(ref['length'])
 p,f,version=lm.decode_bda(raw);world=p.astype(float)@matrix[:3,:3].T+matrix[:3,3];measure=lm.measure(world,f);registered=world@G[:3,:3].T+G[:3,3];sp=np.array(measure['area_weighted_surface_centroid'])@G[:3,:3].T+G[:3,3]
 parts=[p for p in target['parts']if p.get('conceptId')=='ZA:First rib.'+suffix and not vis[p['id']]['suppressed']and not vis[p['id']]['defaultHidden']];assert len(parts)==2
 tv=[];tf=[];controls=[];offset=0
 for part in parts:
  mesh=reader.read_part(target,part,frame='active-male-metres');tv.append(mesh.vertices);tf.append(mesh.faces+offset);offset+=len(mesh.vertices);controls.append(dict(id=part['id'],bankPath=mesh.source['bank_path'],bankSha256=mesh.source['bank_sha256'],record=part))
 tv=np.concatenate(tv);tf=np.concatenate(tf);tm=lm.measure(tv,tf);label=side+' first rib';arrays[side+'_source_registered_vertices_m']=registered;arrays[side+'_source_faces']=f;arrays[side+'_target_vertices_m']=tv;arrays[side+'_target_faces']=tf
 rows.append(dict(label=label,source_registered_m=sp.tolist(),target_m=tm['area_weighted_surface_centroid'],qualification='Complete named first-rib area-weighted exterior surface region centroids. Source whole named bone versus active target bone/cartilage patches. Not a subclavian groove, foraminal point, contact patch, or vessel anchor.',sourceControl=dict(id=sid,geometryKey=key,bdaSlice=ref,bdaSha256=hashlib.sha256(raw).hexdigest(),bdaVersion=version,worldMatrix=matrix.tolist(),faceCount=len(f),sourceMeasure=measure),targetControls=controls,targetMeasure=tm))
labels=[r['label']for r in rows];p=np.array([r['source_registered_m']for r in rows]);q=np.array([r['target_m']for r in rows]);basefield=old['fields'][0];base,_,_=compose_compact_radial_maps(p,[basefield])
root_source=np.array(old['interface_constraint']['source']['center_m']);root_target=np.array(old['interface_constraint']['target']['center_m']);root_evidence=np.load(B/'regional-interface-v1/aortic-root/evidence.npz');root_points=root_evidence['source-proximal-annular-cap_boundary_0_points_m'];root_core=np.vstack([root_source,root_points]);root_mapped,_,_=compose_compact_radial_maps(root_core,[basefield]);assert np.linalg.norm(root_mapped[0]-root_target)<1e-12
fieldspecs=[];maxnorm=[];metadata=[]
for i,label in enumerate(labels):
 center=base[i];minimum_root_distance=float(np.linalg.norm(root_mapped-center,axis=1).min());support=min({'T1':.130,'T7':.105,'T12':.160}.get(label,.115),minimum_root_distance-.006);core=min(.020,support*.2);assert support>core+.025
 fieldspecs.append(dict(center=center.tolist(),core_m=core,support_m=support,kind='translation'))
 maximum=(support-core)*.74/2.109375;maxnorm.append(maximum);metadata.append(dict(fieldIndex=i+1,label=label,coordinates='after the original aortic field; every correction is composed in this order',minimumDistanceToPreservedAorticOuterRim_m=minimum_root_distance,rootSupportMargin_m=minimum_root_distance-support,translationNormMaximum_m=maximum,guideQualification=rows[i]['qualification']))
maxnorm=np.array(maxnorm)
def makefields(x):
 z=np.asarray(x).reshape(-1,3);d=maxnorm[:,None]*z/np.sqrt(1+np.sum(z*z,axis=1))[:,None];return [basefield]+[{**f,'translation_m':di.tolist()}for f,di in zip(fieldspecs,d)]
initial=q-base;norm=np.linalg.norm(initial,axis=1);initial*=np.minimum(1,.9*maxnorm/np.maximum(norm,1e-30))[:,None];z0=initial/np.sqrt(np.maximum(maxnorm**2-np.sum(initial*initial,axis=1),1e-30))[:,None]
def residual(x):return (compose_compact_radial_maps(p,makefields(x))[0]-q).ravel()
solve=least_squares(residual,z0.ravel(),max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12);fields=makefields(solve.x);mapped,jac,proof=compose_compact_radial_maps(p,fields);errors=np.linalg.norm(mapped-q,axis=1)
for i,row in enumerate(rows):row.update(aorticOnly_m=base[i].tolist(),candidate_m=mapped[i].tolist(),aorticOnlyResidual_m=float(np.linalg.norm(base[i]-q[i])),candidateResidual_m=float(errors[i]),candidateMinusAorticOnly_m=(mapped[i]-base[i]).tolist())
root_after,rj,rproof=compose_compact_radial_maps(root_core,fields);assert np.array_equal(root_after,root_mapped);assert errors.max()<.0015;assert solve.success
assert proof['analytic_determinant_lower_bound']>.01
source=BloodSource(B/'source-v2');assert len(source.instances)==1724;minimum_det=np.inf;maxdisp=0.;affected=0;sample_count=0;finite_normals=True;sample_rows=[];witness_preservation=[]
for i,ins in enumerate(source.instances):
 mesh=source.geometry(ins['geometryId']);v=mesh['positions'];ids=np.unique(np.linspace(0,len(v)-1,min(16,len(v))).astype(int));matrix=np.array(ins['worldMatrix']).reshape(4,4,order='F');affine=G[:3,:3]@matrix[:3,:3];t=G[:3,:3]@matrix[:3,3]+G[:3,3];registered=v[ids].astype(float)@affine.T+t;out,j,pr=compose_compact_radial_maps(registered,fields);det=np.linalg.det(j);assert np.isfinite(out).all()and np.isfinite(j).all()and det.min()>=pr['analytic_determinant_lower_bound']-1e-11
 normals=mesh['viewer_normals'][ids].astype(float)@np.linalg.inv(affine);transported=np.einsum('ni,nij->nj',normals,np.linalg.inv(j));assert np.isfinite(transported).all();displacement=float(np.linalg.norm(out-registered,axis=1).max());minimum_det=min(minimum_det,float(det.min()));maxdisp=max(maxdisp,displacement);sample_count+=len(ids);affected+=int(displacement>0)
 sample_rows.append(dict(instanceIndex=i,sourceInstanceId=ins['id'],sampleVertices=len(ids),minimumJacobianDeterminant=float(det.min()),maximumDisplacement_m=displacement))
arrays.update(guide_source_registered_m=p,guide_target_m=q,guide_aortic_only_m=base,guide_candidate_m=mapped,root_before_m=root_core,root_aortic_only_m=root_mapped,root_candidate_m=root_after)
contract={**old,'schema':'biodigital-common-regional-correction/v2','status':'numerically valid common field contract; full candidate geometry and anatomical acceptance pending','fields':fields,'analytic_proof':proof,'guide_metadata':rows,'field_metadata':metadata,'root_preservation':dict(aorticOuterRimCenterError_m=float(np.linalg.norm(root_after[0]-root_target)),entireOuterRimExactlyEqualToAorticOnly=True,allSubsequentSupportSpheresExcludeOuterRim=True,sourceRootNormsAndCalibreAfterAorticOnlyPreserved=True),'finite_tolerances':dict(maximumCoarseBoneRegionCentroidResidual_m=.0015,maximumAorticOuterRimCenterResidual_m=1e-10,minimumContinuousJacobianDeterminant=.01,minimumSampledJacobianDeterminant=.01,maximumSubsequentRootCoordinateChange_m=0.),'optimization':dict(method='bounded-norm parameterization of common compact C2 translations; least squares on five qualified bone region guides',success=bool(solve.success),nfev=solve.nfev,maximumGuideResidual_m=float(errors.max()),noWitnessToBoneCentroidCorrespondenceInvented=True),'all_instance_sample_validation':dict(instances=1724,sampledVertices=sample_count,instancesMovedInSample=affected,minimumJacobianDeterminant=minimum_det,maximumTotalDisplacement_m=maxdisp,normalTransportFinite=True,discreteFaceQualityPending=True),'binding_evidence':dict(sourceManifestSha256=sha(B/'source-v2/manifest.json'),globalRegistrationSha256=sha(Gpath),boneLandmarksSha256=sha(lp),sourceCatalogueSha256=sha(source_catalog),targetCatalogueSha256=sha(cp),visibilitySha256=sha(vispath),savedCrossingAggregationSha256=sha(B/'regional-fit-v1/skeleton-samples-v1/aggregation-v1/summary.json'),builderSha256=sha(__file__)),'rationale':'The aortic-only field worsened sampled torso crossings 663 to 707. This candidate preserves its entire aligned aortic outer rim and adds a common thoracic/cervicothoracic correction constrained by actual bone region guides. Crossing reduction remains to be measured.','required_followups':['Apply all six ordered fields to every selected source instance from immutable source-v2 after immutable G; do not apply v2 onto v1 baked meshes.','Transport normals with the inverse transpose of the full chain-rule Jacobian.','Run all-source nonzero-face, orientation, ancestry and exact aortic-rim preservation checks.','Replay the same saved skeleton sample ownership with the new candidate and compare thoracic/trunk witnesses.','Review actual surfaces and sections; bone centroids are coarse guides, not certified vascular corridors or foramina.','No default installation, anatomy acceptance or skeleton movement.'],'installed':False,'anatomical_acceptance':False,'elapsed_seconds':time.monotonic()-start}
contract['exact_core_translation']=False
contract['outside_300mm_support_unchanged']=False
contract['source_root_orientation_and_calibre_preserved_in_core']=False
contract['root_preservation']['qualification']='The complete measured outer-rim support is exactly preserved after the original translation. No claim is made that its entire surrounding 40 mm ball is rigid under the composition.'
contract['application_order']='G @ immutable source-world coordinates, then all contract.fields in list order. Each field uses coordinates produced by the preceding field. Apply identically to all 1724 instances; no per-vessel masks. Identity outside the union of the ordered support spheres.'
O.mkdir();np.savez_compressed(O/'guide-evidence.npz',**arrays);(O/'guide-validation.json').write_text(json.dumps(dict(rows=rows,fieldMetadata=metadata,optimization=contract['optimization'],allInstanceSamples=contract['all_instance_sample_validation'],sampleInstanceRows=sample_rows,elapsed_seconds=contract['elapsed_seconds']),indent=2)+'\n');(O/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
print(json.dumps(dict(contract=str(O/'contract.json'),fields=len(fields),guideResidualMm={r['label']:r['candidateResidual_m']*1000 for r in rows},aorticOnlyGuideResidualMm={r['label']:r['aorticOnlyResidual_m']*1000 for r in rows},analyticDeterminantLowerBound=proof['analytic_determinant_lower_bound'],minimumSampledDeterminant=minimum_det,instanceSamples=sample_count,elapsedSeconds=contract['elapsed_seconds'])))
