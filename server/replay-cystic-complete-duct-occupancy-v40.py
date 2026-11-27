"""Full continuous occupied-volume audit against the recovered actual duct."""
from pathlib import Path
import sys,time,json,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import topology
from geometry_tools.contact import self_contact_summary,pair_contacts
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.finite_tissue_boolean import make_solid,same_oriented_triangles

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
PRE=B/'biliary-actual-envelope-v39';OUT=B/'cystic-complete-duct-occupancy-v40'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'

def recover(actual,partition,label):
    z=np.load(partition);ids=z['branch_source_face_ids'];shared=z['shared_branch_faces'];take=np.ones(len(ids),bool);take[shared]=False;order=[]
    for sheet in sorted(set((ids//1000000)[take])):order.extend(np.flatnonzero(take&(ids//1000000==sheet)))
    order.extend(shared);order=np.array(order,int);assert len(order)==len(actual.faces)
    roles=ids[order]//1000000;assert set(roles)=={1,3,6}
    pos,inv=np.unique(actual.vertices,axis=0,return_inverse=True);outer=inv[actual.faces[np.isin(roles,[1,6])]];inner=inv[actual.faces[roles==3]]
    ol=boundary_loops(outer);il=boundary_loops(inner);assert len(ol)==len(il)==2
    assert {frozenset(l)for l in ol}=={frozenset(l)for l in il}
    points=pos.tolist();caps=[]
    for ring in ol:
        center=len(points);points.append(pos[ring].mean(0));caps.extend([[center,int(b),int(a)]for a,b in zip(ring,np.roll(ring,-1))])
    points=np.array(points);caps=np.array(caps,int);of=np.r_[outer,caps];iff=np.r_[inner[:,::-1],caps]
    assert same_oriented_triangles(np.r_[of,iff[:,::-1]],np.r_[inv[actual.faces],caps,caps[:,::-1]])
    checks=[];meshes=[]
    for side,faces in [('outer',of),('fluid',iff)]:
        m=Mesh(points,faces,actual.frame);t=topology(m);assert t['closed_oriented_manifold']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==0,(label,side,t)
        sc,_=self_contact_summary(m,max_pairs=500000,max_tests=5000000,stop_on_defect=True);assert sc['complete_within_pair_budget']and sc['crossing_face_pairs']==sc['coplanar_overlap_face_pairs']==0,(label,side,sc)
        m.save(OUT/(label+'-actual-'+side+'-computational-envelope.npz'));meshes.append(m);checks.append(dict(side=side,geometrySha256=m.geometry_sha,topology=t,selfIntersections=sc))
    origin=np.array([0.,1.25,0.]);os,_=make_solid((points-origin)*1000,of);fs,_=make_solid((points-origin)*1000,iff);ws,_=make_solid((actual.vertices-origin)*1000,actual.faces);difference=os-fs;error=float(((difference-ws)+(ws-difference)).volume());escape=float((fs-os).volume());assert error<1e-9 and escape==0,(label,error,escape)
    return meshes,dict(actualBankSha256=actual.source['bank_sha256'],actualGeometrySha256=actual.geometry_sha,partitionSha256=digest_file(partition),sourceBoundaryRoleTags=[1,3,6],actualOwnerStartsAtSharedRing7=True,numericalCapTriangles=len(caps),numericalCapsCancelExactlyFromWall=True,wallReconstructionSymmetricVolume_mm3=error,fluidOutsideOuter_mm3=escape,envelopes=checks)

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not OUT.exists();OUT.mkdir()
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';origin=np.array([0.,1.25,0.])
    duct=reader.read_part(a,parts['TA98HEPATIC:male:biliary'],frame=frame,catalogue_path=CAT);outer=Mesh.load(PRE/'actual-biliary-outer-computational-envelope.npz');fluid=Mesh.load(PRE/'actual-biliary-fluid-computational-envelope.npz');proof=json.loads((PRE/'proof.json').read_text());assert proof['passed']and duct.geometry_sha==proof['actualWallGeometrySha256']and duct.source['bank_sha256']==proof['actualBankSha256']
    solids={}
    for name,m in [('outer',outer),('fluid',fluid),('wall',duct)]:solids[name]=make_solid((m.vertices-origin)*1000,m.faces)[0]
    spec=importlib.util.spec_from_file_location('v10',ROOT/'server/correct-cystic-profiles-and-host-mouths-v10.py');v10=importlib.util.module_from_spec(spec);spec.loader.exec_module(v10)
    results=[]
    for label,pid,partition,sections in [('installed-artery','TA98CYSTICBLOOD:male:artery',B/'owned-parent-cells-v18/artery/partition.npz',B/'admissible-departure-v14/sections.npz'),('fixed-installed-vein','TA98CYSTICBLOOD:male:vein',B/'vein-recipient-partition-v28/vein/partition.npz',B/'vein-recipient-course-v27/sections.npz')]:
        actual=reader.read_part(a,parts[pid],frame=frame,catalogue_path=CAT);envelopes,recovery=recover(actual,partition,label);os,fs=[make_solid((m.vertices-origin)*1000,m.faces)[0]for m in envelopes]
        occupancy={name:float((os^s).volume())for name,s in solids.items()};occupancy['bloodInsideBileFluid_mm3']=float((fs^solids['fluid']).volume())
        contacts,check=pair_contacts(actual,duct,0.,max_pairs=1000000,max_tests=8000000);assert check['complete_within_pair_budget'];atomic_npz(OUT/(label+'-global-actual-wall-crossings.npz'),**contacts)
        src=np.load(sections);profiles=src['outer_sections'];v,f,_=v10.swept_profiles(profiles);whole=Mesh(v,f,frame);whole.save(OUT/(label+'-full-source-course-outer-computational-envelope.npz'));t=topology(whole);assert t['closed_oriented_manifold']and len(t['components'])==1 and t['components'][0]['closed_orientable_genus']==0
        sc,_=self_contact_summary(whole,max_pairs=500000,max_tests=5000000,stop_on_defect=True);assert sc['complete_within_pair_budget']and sc['crossing_face_pairs']==sc['coplanar_overlap_face_pairs']==0
        whole_solid=make_solid((v-origin)*1000,f)[0];full_volume=float((whole_solid^solids['outer']).volume())
        # Replay every finite interval's full volume, not centers or a sampled
        # profile disk. Independent interval unions must agree with the full
        # continuous source-course Boolean (within arithmetic residue).
        intervals=[];interval_union=None
        for j in range(len(profiles)-1):
            iv,iff,_=v10.swept_profiles(profiles[j:j+2]);s=make_solid((iv-origin)*1000,iff)[0];hit=s^solids['outer'];volume=float(hit.volume())
            if volume>1e-12:intervals.append(dict(interval=j,outerEnvelopeOccupation_mm3=volume))
            interval_union=hit if interval_union is None else interval_union+hit
        interval_error=float(((interval_union-(whole_solid^solids['outer']))+((whole_solid^solids['outer'])-interval_union)).volume());assert interval_error<1e-8
        results.append(dict(id=pid,label=label,recovery=recovery,actualBranchEnvelopeOccupation_mm3=occupancy,completeActualGlobalDuctWallCrossings=check,fullOriginalSourceCourse=dict(sectionsSha256=digest_file(sections),stationCount=len(profiles),fixedRootRings0To7Included=True,fullContinuousOuterVolumeOccupation_mm3=full_volume,occupiedFiniteIntervals=intervals,intervalUnionVsWholeCourseSymmetricVolume_mm3=interval_error,selfIntersections=sc,topology=t),passedNoUnapprovedDuctOccupation=all(v==0 for v in occupancy.values())and full_volume==0 and check['crossing_face_pairs']==check['coplanar_overlap_face_pairs']==0))
        atomic_json(OUT/(label+'-occupancy.json'),results[-1]);print(json.dumps(dict(label=label,actualOccupation=occupancy,fullCourseOccupation_mm3=full_volume,crossings=check['crossing_face_pairs'],occupiedIntervals=intervals)),flush=True)
    result=dict(schema='complete-cystic-actual-duct-occupancy/v40',sourceCatalogueSha256=EXPECTED,actualDuctOracleProofSha256=digest_file(PRE/'proof.json'),wholeGlobalActualDuctBoundaryUsed=True,bothInnerAndOuterActualWallCrossingsComplete=True,allCourseIntervalsAndActualBranchTrianglesChecked=True,oracleCapsAreComputationalOnly=True,noGeometryFitOrReplacement=True,geometryModified=False,installed=False,results=results,oracleRecoveredAndValidated=True,allProtectedCoursesClear=all(r['passedNoUnapprovedDuctOccupation']for r in results),elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',result);assert digest_file(CAT)==EXPECTED;print(json.dumps(dict(oracleValid=True,allProtectedCoursesClear=result['allProtectedCoursesClear'],elapsedSeconds=result['elapsedSeconds'])))
if __name__=='__main__':main()
