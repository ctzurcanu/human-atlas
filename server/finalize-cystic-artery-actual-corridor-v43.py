"""Final 60s stage: bounded proper-translation correction and emitted replay."""
from pathlib import Path
import sys,time,json,types,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.open3d_backend import SurfaceIndex,library
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.contact import pair_contacts
from geometry_tools.review import prepare

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1';OUT=B/'artery-final-source-frame-v43';CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345';FRAME='model:local-ta98:assembled-metres'

def load(name):
    sp=importlib.util.spec_from_file_location(name.replace('-','_'),ROOT/'server'/name);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def refine(began):
    v10=load('correct-cystic-profiles-and-host-mouths-v10.py');s=np.load(B/'admissible-departure-v14/sections.npz');old=s['outer_sections'];source_centers=s['centers'];base=np.load(B/'artery-valid-envelope-corridor-v41/sections.npz')['centers'];centers=base.copy();N=len(base);variable=np.arange(8,22);origin=source_centers[0];oracle=Mesh.load(B/'biliary-actual-envelope-v39/actual-biliary-outer-computational-envelope.npz');ds=make_solid((oracle.vertices-origin)*1000,oracle.faces)[0];index=SurfaceIndex(oracle);oracle_normals=np.cross(oracle.triangles[:,1]-oracle.triangles[:,0],oracle.triangles[:,2]-oracle.triangles[:,0]);oracle_normals/=np.linalg.norm(oracle_normals,axis=1)[:,None]
    bary=np.array([[1,0,0],[0,1,0],[0,0,1],[.5,.5,0],[0,.5,.5],[.5,0,.5],[1/3,1/3,1/3],[.6,.2,.2],[.2,.6,.2],[.2,.2,.6]])
    D=np.diff(np.eye(N),n=2,axis=0);history=[];passed=False
    atomic_json(OUT/'assembly-contract.json',dict(schema='final-local-artery-translation-contract/v43',sourceCatalogueSha256=EXPECTED,movableRings=variable.tolist(),affectedSideIntervals=list(range(7,22)),protectedRootRings=list(range(8)),protectedTerminalRings=list(range(22,N)),sourceProfilesReceiveProperTranslationsOnly=True,sourceCalibreAndShapeLiteral=True,maximumAdditionalSectionTranslation_m=.00075,actualGlobalDuctOracleProofSha256=digest_file(B/'biliary-actual-envelope-v39/proof.json'),allOtherSourcesAndVeinLiteral=True,stage=3,maximumSavedStages=3,installed=False))
    for iteration in range(16):
        if time.monotonic()-began>10:break
        outer=old+(centers-source_centers)[:,None];v,f,_=v10.swept_profiles(outer);bo=v10.previous.solid(v,f,origin);volume=float((bo^ds).volume());side=f[:(N-1)*36];sel=side[7*36:22*36];tri=v[sel];q=np.einsum('bc,tcj->tbj',bary,tri).reshape(-1,3);hit,_=index.closest_points(q,certify=False);inside=index.scene.compute_occupancy(library().core.Tensor(((q-index.origin)/index.scale).astype(np.float32)),nthreads=2).numpy()>.5;signed=hit['distance_m']*np.where(inside,-1,1);bad=signed<.0002;history.append(dict(iteration=iteration,completeOuterEnvelopeOccupation_mm3=volume,minimumProposedSampleClearance_m=float(signed.min()),additionalSectionDisplacementMax_m=float(np.linalg.norm(centers-base,axis=1).max())))
        if volume==0 and signed.min()>.00009:passed=True;break
        if not bad.any():break
        toward=q-hit['points'];length=np.linalg.norm(toward,axis=1);direction=np.zeros_like(toward);valid=length>1e-12;direction[valid]=toward[valid]/length[valid,None];direction[inside]*=-1;direction[~valid]=oracle_normals[hit['face_ids'][~valid]]
        displacement=direction[bad]*(.00023-signed[bad])[:,None]
        rings=sel//18;A=np.zeros((len(sel),len(bary),N))
        for k in range(3):
            for j in range(len(sel)):A[j,:,rings[j,k]]+=bary[:,k]
        A=A.reshape(-1,N)[bad][:,variable];lhs=np.r_[A,.7*D[:,variable],.1*np.eye(len(variable))];rhs=np.r_[displacement,np.zeros((len(D)+len(variable),3))];step=np.linalg.lstsq(lhs,rhs,rcond=None)[0];scale=min(1.,.00024/max(np.linalg.norm(step,axis=1).max(),1e-30));centers[variable]+=step*scale
        if np.linalg.norm(centers-base,axis=1).max()>.00075:break
    atomic_json(OUT/'bounded-outward-translation-history.json',dict(history=history,passedProposalVolume=passed,elapsedSeconds=time.monotonic()-began,installed=False));atomic_npz(OUT/'refined-center-field.npz',centers=centers)
    return passed

def construct_and_pack():
    text=(ROOT/'server/finish-cystic-artery-profile-transport-v42.py').read_text().replace('artery-source-profile-transport-v42','artery-final-source-frame-v43/candidate')
    text=text.replace("candidate=np.load(B/'artery-valid-envelope-corridor-v41/sections.npz')","candidate=np.load(B/'artery-final-source-frame-v43/refined-center-field.npz')")
    text=text.replace('for alpha in [0.,.25,.5,.75]:','for alpha in [0.]:').replace("weights.insert(0,('source-root-frame-smoothly-transported-after-ring11',taper))","# Source profile frames stay literal; no tangent reorientation.")
    text=text.replace('stage=2,maximumSavedStages=3','stage=3,maximumSavedStages=3')
    mod=types.ModuleType('final_source_frame_construct');mod.__file__=str(ROOT/'server/finish-cystic-artery-profile-transport-v42.py');exec(compile(text,mod.__file__,'exec'),mod.__dict__);mod.main()

def triangle_word_lookup(mesh):
    words=mesh.attributes['normals'][mesh.faces].astype(np.int64)*mesh.attributes['normalSigns'][mesh.faces,None];rows={}
    for tri,word in zip(mesh.triangles,words):
        keys=[tuple(map(tuple,np.roll(tri,k,axis=0)))for k in range(3)];k=min(range(3),key=keys.__getitem__);assert keys[k]not in rows,'Duplicate oriented source triangle';rows[keys[k]]=np.roll(word,k,axis=0)
    return rows

def literal_surface_words(first,second):
    a=triangle_word_lookup(first);b=triangle_word_lookup(second)
    return a.keys()==b.keys()and all(np.array_equal(a[k],b[k])for k in a)

def emitted(began):
    folder=OUT/'candidate/bank/artery';manifest=json.loads((folder/'manifest.json').read_text());records=json.loads((folder/'records.json').read_text());a=json.loads(CAT.read_text());current={p['id']:p for p in a['parts']};reader=AtlasReader();bank=folder/'cells-and-regions.bin';catalogue={'chunks':[dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')))]};r={pid:dict(p,chunk=0)for pid,p in records.items()};parent=reader.read_part(catalogue,r['TA98HEPATIC:male:artery'],frame=FRAME);branch=reader.read_part(catalogue,r['TA98CYSTICBLOOD:male:artery'],frame=FRAME);oldparent=reader.read_part(a,current['TA98HEPATIC:male:artery'],frame=FRAME,catalogue_path=CAT)
    words=lambda m:m.attributes['normals'][m.faces].astype(np.int64)*m.attributes['normalSigns'][m.faces,None]
    parentliteral=literal_surface_words(parent,oldparent);assert parentliteral,'Protected complete parent oriented triangles/corner words changed'
    z=np.load(OUT/'candidate/owned/artery/partition.npz');owners=z['parent_owner_ids'];pids=z['parent_source_face_ids'];shared=set(map(int,z['shared_parent_faces']));take=np.array([i not in shared for i in range(len(pids))]);order=[]
    for owner in sorted(set(owners[take])):
        for sheet in sorted(set((pids//1000000)[take&(owners==owner)])):order.extend(np.flatnonzero(take&(owners==owner)&(pids//1000000==sheet)))
    order.extend(z['shared_parent_faces']);order=np.array(order);native=np.isin(pids[order]//1000000,[0,2]);assert np.array_equal(z['parent_native_corner_words'][order[native]],words(parent)[native])
    regionchecks=[];count=0
    for pid,p in r.items():
        if pid in ['TA98HEPATIC:male:artery','TA98CYSTICBLOOD:male:artery']:continue
        m=reader.read_part(catalogue,p,frame=FRAME);old=reader.read_part(a,current[pid],frame=FRAME,catalogue_path=CAT);literal=literal_surface_words(m,old);count+=len(m.faces);regionchecks.append(dict(id=pid,triangles=len(m.faces),literalOldGeometryAndCornerWords=literal,oldSurfaceSha256=old.surface_sha,newSurfaceSha256=m.surface_sha));assert literal
    assert count==len(parent.faces)and len(regionchecks)==10 and len(records)==12
    replay=load('replay-cystic-complete-duct-occupancy-v40.py');replay.OUT=OUT/'emitted-validation';replay.OUT.mkdir(exist_ok=True);envelopes,recovery=replay.recover(branch,OUT/'candidate/owned/artery/partition.npz','candidate-artery');origin=np.array([0.,1.25,0.]);os=make_solid((envelopes[0].vertices-origin)*1000,envelopes[0].faces)[0];fs=make_solid((envelopes[1].vertices-origin)*1000,envelopes[1].faces)[0];occupation={}
    for side in ['outer','fluid']:
        duct=Mesh.load(B/('biliary-actual-envelope-v39/actual-biliary-'+side+'-computational-envelope.npz'));ds=make_solid((duct.vertices-origin)*1000,duct.faces)[0];occupation[side]=float((os^ds).volume())
    assert all(v==0 for v in occupation.values())
    checks=[]
    for pid in ['TA98CYSTICBLOOD:male:vein','TA98HEPATIC:male:portal','TA98HEPATIC:male:hepatic_vein','TA98HEPATIC:male:biliary','FJ3080','TA98PART:FJ3079:A05.8.01.061']:
        m=reader.read_part(a,current[pid],frame=FRAME,catalogue_path=CAT);data,c=pair_contacts(branch,m,0.,max_pairs=500000,max_tests=5000000);assert c['complete_within_pair_budget']and c['crossing_face_pairs']==c['coplanar_overlap_face_pairs']==0,(pid,c);atomic_npz(replay.OUT/(pid.replace(':','_')+'-contacts.npz'),**data);checks.append(dict(id=pid,geometrySha256=m.geometry_sha,bankSha256=m.source['bank_sha256'],**c))
    inner=Mesh.load(B/'first-similarity/fixed-mucosa-free-surface.npz');_,c=pair_contacts(branch,inner,0.,max_pairs=500000,max_tests=5000000);assert c['complete_within_pair_budget']and c['crossing_face_pairs']==c['coplanar_overlap_face_pairs']==0;checks.append(dict(id='GB-innermost',**c))
    currentbranch=reader.read_part(a,current['TA98CYSTICBLOOD:male:artery'],frame=FRAME,catalogue_path=CAT);part=np.load(B/'owned-parent-cells-v18/artery/partition.npz');oldids=part['branch_source_face_ids'];bs=part['shared_branch_faces'];tk=np.ones(len(oldids),bool);tk[bs]=False;oldorder=[]
    for side in sorted(set((oldids//1000000)[tk])):oldorder.extend(np.flatnonzero(tk&(oldids//1000000==side)))
    oldorder.extend(bs);oldlookup=triangle_word_lookup(currentbranch);ids=z['branch_source_face_ids'];bs=z['shared_branch_faces'];tk=np.ones(len(ids),bool);tk[bs]=False;neworder=[]
    for side in sorted(set((ids//1000000)[tk])):neworder.extend(np.flatnonzero(tk&(ids//1000000==side)))
    neworder.extend(bs);literaloutside=0;changed=[]
    for fid,tri,word in zip(ids[neworder],branch.triangles,words(branch)):
        role=int(fid//1000000);local=int(fid%1000000);affected=role in [1,3]and local<(71-1)*36 and 7<=local//36<=21
        if not affected:
            keys=[tuple(map(tuple,np.roll(tri,k,axis=0)))for k in range(3)];k=min(range(3),key=keys.__getitem__);assert keys[k]in oldlookup and np.array_equal(np.roll(word,k,axis=0),oldlookup[keys[k]]),'Protected distal/root branch triangle or normals changed';literaloutside+=1
        else:changed.append(int(fid))
    proof=dict(schema='corrected-cystic-artery-complete-emitted/v43',sourceCatalogueSha256=EXPECTED,bank=str(bank),bankSha256=digest_file(bank),records=list(records),regionalControls=10,physicalOwners=2,allTenControlsRebuiltOrProvenLiteral=True,protectedNativeParentGeometryAndNormalsLiteral=parentliteral,parentArrayOrderMayDifferButExactOrientedTrianglesAndCornerWordsLiteral=True,regionalChecks=regionchecks,protectedBranchFacesAndCornerWordsLiteral=literaloutside,affectedBranchSourceFaces=changed,completeActualOuterEnvelopeDuctOccupation_mm3=occupation,actualOuterFluidWallReconstruction=recovery,globalActualWallContextChecks=checks,wholeEmittedGatesPassed=True,sourceRootAndTipLiteral=True,noRadiusShrink=True,allOtherRecordsLiteral=True,installed=False,nativeBlenderReviewPending=True,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'emitted-final-proof.json',proof)
    return parent,branch,proof

def queue_review(parent,branch,proof):
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();entries=[]
    for name,m,color in [('Corrected emitted cystic artery',branch,[.86,.12,.10,1]),('Literal native right-hepatic arterial parent',parent,[.72,.20,.13,.32])]:
        path=OUT/(name.replace(' ','-')+'.npz');m.save(path);entries.append(dict(mesh=str(path),name=name,color=color))
    for pid,name,color in [('TA98CYSTICBLOOD:male:vein','Fixed installed cystic vein',[.12,.24,.86,1]),('TA98HEPATIC:male:biliary','Actual fixed biliary paired wall',[.65,.43,.27,.4])]:
        m=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT);path=OUT/(name.replace(' ','-')+'.npz');m.save(path);entries.append(dict(mesh=str(path),name=name,color=color))
    sec=np.load(OUT/'candidate/sections.npz');centers=sec['centers'];t=centers[8]-centers[6];t/=np.linalg.norm(t);vein=Mesh.load(OUT/'Fixed-installed-cystic-vein.npz');normal=np.cross(t,np.array([0.,1.,0.]));normal/=np.linalg.norm(normal)
    views=[dict(name='whole-corrected-course-surface-front',direction=[0,0,1]),dict(name='whole-corrected-course-surface-oblique',direction=[1,.3,1]),dict(name='whole-course-longitudinal-section',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=centers[15].tolist())),dict(name='native-shared-root-longitudinal-section',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=centers[7].tolist()))]
    review=prepare(None,dict(frame=FRAME,meshes=entries,views=views,resolution=1000,wall_seconds=40));atomic_json(OUT/'emitted-native-review.json',dict(**review,bankSha256=proof['bankSha256'],sourceCatalogueSha256=EXPECTED,requiredNativeReviewViews=[v['name']for v in views],installed=False));return review

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED
    if '--finish-emitted-only' in sys.argv:
        import copy
        folder=OUT/'candidate/bank/artery';manifest=json.loads((folder/'manifest.json').read_text());assert manifest['emittedGeometryChecksPassed']and digest_file(folder/'cells-and-regions.bin')==manifest['bankSha256']
        a=json.loads(CAT.read_text());current={p['id']:p for p in a['parts']};generated=json.loads((folder/'records.json').read_text());geometry={'positions','normals','normalSigns','indices','indexCount','vertexCount','geometryOrigin','triangles','bounds','tissueDrawGroups','pairedTissue','sectionSolid'};records={}
        for pid,p in generated.items():
            q=copy.deepcopy(current[pid])
            for k in ['geometryScale','geometryInstance','sourceOffset','uvs','colors']:q.pop(k,None)
            for k in geometry:
                if k in p:q[k]=p[k]
            records[pid]=q
        atomic_json(folder/'records.json',records,overwrite=True);atomic_json(OUT/'candidate/metadata-preservation.json',dict(sourceCatalogueSha256=EXPECTED,allSourceRecordIds=list(records),currentSourceMetadataRetained=True,geometryFieldsOnlyReplaced=sorted(geometry),termCountUnchanged=True,bankUnmodifiedDuringBookkeepingResume=True),overwrite=True)
        parent,branch,proof=emitted(began);review=queue_review(parent,branch,proof);atomic_json(OUT/'bounded-stage3-resume.json',dict(firstCommandElapsedSeconds=22.98,resumedMetadataAndReplaySeconds=time.monotonic()-began,totalStageRuntimeSeconds=22.98+time.monotonic()-began,noGeometryRetry=True,installed=False));print(json.dumps(dict(passed=True,bank=proof['bank'],review=review,elapsedSeconds=time.monotonic()-began)));return
    assert not OUT.exists();OUT.mkdir()
    if not refine(began):print(json.dumps(dict(stopped='No complete admissible corrected course within final bounded search',elapsedSeconds=time.monotonic()-began)));return
    construct_and_pack();folder=OUT/'candidate/bank/artery'
    if not (folder/'manifest.json').exists():print(json.dumps(dict(stopped='Final profile/context/pack gates rejected; no further stage authorized',elapsedSeconds=time.monotonic()-began)));return
    parent,branch,proof=emitted(began)
    if time.monotonic()-began<52:review=queue_review(parent,branch,proof);print(json.dumps(dict(passed=True,bank=proof['bank'],review=review,elapsedSeconds=time.monotonic()-began)))
    else:print(json.dumps(dict(passed=True,reviewHeldForDeadline=True,elapsedSeconds=time.monotonic()-began)))
    assert digest_file(CAT)==EXPECTED
if __name__=='__main__':main()
