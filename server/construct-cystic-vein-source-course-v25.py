"""Preserve the passing source course and transport a smooth lateral displacement."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
from scipy.optimize import lsq_linear
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.review import prepare
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
sp=importlib.util.spec_from_file_location('v23',ROOT/'server/plan-cystic-vein-local-entry-v23.py');v23=importlib.util.module_from_spec(sp);sp.loader.exec_module(v23);v14=v23.v14;v10=v14.v10
WORK=v14.WORK;OUT=WORK/'vein-source-course-v25';CURRENT=v23.CURRENT;EXPECTED=v23.EXPECTED;FRAME=v14.FRAME


def main():
    began=time.monotonic();assert digest_file(CURRENT)==EXPECTED;assert not OUT.exists();OUT.mkdir();a=json.loads(CURRENT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();artery=reader.read_part(a,parts['TA98CYSTICBLOOD:male:artery'],frame=FRAME,catalogue_path=CURRENT);ab=BVH(artery);cover=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz');cb=BVH(cover);inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');ib=BVH(inner);feas=json.loads((WORK/'vein-local-entry-v23/feasibility.json').read_text());assert feas['feasible'];endpoint=np.load(WORK/'vein-local-entry-v23/selected-endpoint.npz')
    old=np.load(WORK/'cystic-vein-v16/sections.npz');old_outer=old['outer_sections'];old_centers=old['centers'];N=len(old_centers);fixed=7
    # Never reproject the portal departure onto the distant recipient. The
    # previous surface-projection trial moved it 14 mm despite a 4 mm corridor.
    # A single proper displacement field acts on the passing v16 course.
    tangent=np.gradient(old_centers,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None];delta=endpoint['shift'];reference=delta/np.linalg.norm(delta)
    lateral=reference[None]-tangent*(tangent@reference)[:,None];lateral/=np.linalg.norm(lateral,axis=1)[:,None]
    end_amplitude=float(lateral[-1]@delta);end_remainder=delta-end_amplitude*lateral[-1];blend=np.clip((np.arange(N)-50)/20,0,1);blend=blend**3*(10-15*blend+6*blend*blend);extra=blend[:,None]*end_remainder
    az=np.load(WORK/'admissible-departure-v14/sections.npz');av,af,_=v10.swept_profiles(az['outer_sections'][7:]);origin=np.array(parts['TA98CYSTICBLOOD:male:artery']['geometryOrigin']);av=(av-origin).astype('<f4').astype(float)+origin;outer_artery=Mesh(av,af,FRAME)
    disk_faces=np.array([[18,i,(i+1)%18]for i in range(18)]);offsets=np.arange(21)*.0002;lower=np.zeros(N);selection=[];failed=None
    # Finite profile-aware station feasibility supplies lower bounds to one
    # smooth convex course solve, rather than repeated whole-mesh fitting.
    for j in range(fixed+1,N-1):
        found=None
        if time.monotonic()-began>38:failed=dict(station=j,reason='Construction feasibility time budget reached before completing every station');break
        for k,amount in enumerate(offsets):
            center=old_centers[j]+extra[j]+amount*lateral[j];profile=old_outer[j]+center-old_centers[j]
            if ab.closest_point(center)['distance_m']<.00025:continue
            inside,amb=ray_contains(center[None],outer_artery.triangles)
            if inside[0]or amb[0]:continue
            disk=Mesh(np.r_[profile,center[None]],disk_faces,FRAME);_,c=pair_contacts(disk,artery,.00025,second_bvh=ab,max_pairs=250000,max_tests=2000000);assert c['complete_within_pair_budget']
            if c['returned_triangle_pairs']:continue
            _,lumen=pair_contacts(disk,inner,0.,second_bvh=ib,max_pairs=250000,max_tests=2000000);assert lumen['complete_within_pair_budget']
            if lumen['crossing_face_pairs']or lumen['coplanar_overlap_face_pairs']:continue
            found=float(amount);break
        if found is None:failed=dict(station=j,reason='No source-profile-safe station within the single 4 mm lateral corridor');break
        lower[j]=found;selection.append(dict(station=j,minimumTestedLateralOffset_m=found))
    if failed:
        atomic_json(OUT/'gates.json',dict(schema='cystic-source-course/v25',sourceCatalogueSha256=EXPECTED,stage=1,maximumStages=3,stoppedBeforeConstruction=True,failure=failed,stationFeasibility=selection,installed=False,elapsedSeconds=time.monotonic()-began));print(json.dumps(dict(stoppedBeforeConstruction=True,failure=failed,elapsedSeconds=time.monotonic()-began)));return
    # One least-squares solution: stay above the explicit profile-clearance
    # bounds while keeping the displacement small and its turns smooth.
    amplitude=np.zeros(N);amplitude[-1]=end_amplitude;variable=np.arange(fixed+1,N-1);D=np.diff(np.eye(N),n=2,axis=0);A=np.r_[np.eye(N),np.sqrt(96)*D]
    fixed_mm=np.zeros(N);fixed_mm[-1]=end_amplitude*1000;rhs=-A@fixed_mm;fit=lsq_linear(A[:,variable],rhs,bounds=(lower[variable]*1000,4.01),max_iter=60,tol=1e-9);amplitude[variable]=fit.x/1000
    centers=old_centers+amplitude[:,None]*lateral+extra;centers[:fixed+1]=old_centers[:fixed+1];centers[-1]=endpoint['center'];outer=old_outer.copy()
    newtangent=np.gradient(centers,axis=0);newtangent/=np.linalg.norm(newtangent,axis=1)[:,None]
    for j in range(fixed+1,N-1):
        R=v23.rotation(tangent[j],newtangent[j]);outer[j]=(old_outer[j]-old_centers[j])@R.T+centers[j]
    target_end=v14.frame(endpoint['outer_profile'],old_centers[-1]-old_centers[-2])
    outer[-1]=endpoint['outer_profile'];fluid=centers[:,None]+.82*(outer-centers[:,None]);fluid[:fixed+1]=old['fluid_sections'][:fixed+1];fluid[-1]+=target_end[:,2]*.001
    ov,of,_=v10.swept_profiles(outer);iv,iff,_=v10.swept_profiles(fluid);side=of[:-36];fn=np.cross(ov[side[:,1]]-ov[side[:,0]],ov[side[:,2]]-ov[side[:,0]]);radial=ov[side].mean(1)-np.repeat(centers,36,axis=0)[:len(side)]
    if np.median(np.sum(fn*radial,axis=1))<0:side=side[:,::-1]
    nominal=centers[:,None]+.82*(outer-centers[:,None]);prism=StationPrismConstraints(np.array([outer.reshape(-1,3),nominal.reshape(-1,3)]),side,np.repeat(np.arange(N),18),np.ones(N,bool),allow_source_faults=True);bad=np.flatnonzero(~prism.sourceHealthy[0]);invalid=np.unique(np.repeat(np.arange(N),18)[side[bad]].min(1)).tolist();area=np.array([projected_area(p)for p in outer]);oldarea=np.array([projected_area(p)for p in old_outer]);diam=lambda p:np.linalg.norm(p[:,:,None]-p[:,None,:],axis=3).max((1,2));area_error=float(abs(area/oldarea-1).max());diam_error=float(abs(diam(outer)/diam(old_outer)-1).max())
    atomic_npz(OUT/'sections.npz',outer_sections=outer,fluid_sections=fluid,centers=centers,source_outer_sections=old_outer,station_offsets=amplitude,feasibility_lower_bounds=lower,old_tangents=tangent,new_tangents=newtangent)
    proof=dict(schema='cystic-source-course/v25',sourceCatalogueSha256=EXPECTED,stage=1,maximumStages=3,continuousPrismInvalidIntervals=invalid,maximumContinuousPrismDeterminant=float(prism.sourceMaximum.max()),maximumExternalSectionAreaRelativeChange=area_error,maximumExternalSectionDiameterRelativeChange=diam_error,sourceSectionCalibreErrorUpperBound=.04495321282425779+diam_error,sourceCalibreTolerance=.15,fixedPortalRootAndFirst8SectionsLiteral=bool(np.array_equal(outer[:8],old_outer[:8])),distalEntryShift_m=float(np.linalg.norm(centers[-1]-old_centers[-1])),maximumRouteCenterShift_m=float(np.linalg.norm(centers-old_centers,axis=1).max()),stationFeasibility=selection,smoothOffsetSolveSuccess=bool(fit.success),arteryAndAllOrgansFixed=True,installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began)
    if invalid or area_error>.15 or diam_error>.15:
        proof['stoppedBeforeCSG']=True;atomic_json(OUT/'gates.json',proof);print(json.dumps(dict(stoppedBeforeCSG=True,invalidPrismIntervals=invalid,areaError=area_error,diameterError=diam_error,elapsedSeconds=time.monotonic()-began)));return
    root=centers[0];host=Mesh.load(WORK/'fixed-target-21.npz');current_host=reader.read_part(a,parts['TA98HEPATIC:male:portal'],frame=FRAME,catalogue_path=CURRENT);assert current_host.geometry_sha==host.geometry_sha;blood=Mesh.load(WORK/'native-parent-envelopes-v9/vein-blood.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/vein-inner-lineage.json').read_text());_,ho,hb,reconstruction=v10.native_outer(host,blood,lineage,root);bo=v10.previous.solid(ov,of,root);bb=v10.previous.solid(iv,iff,root);joined=(ho+bo)-(hb+bb);jv,jf=v10.previous.arrays(joined,root);network=Mesh(jv,jf,FRAME);network.save(OUT/'vein-joined-wall.npz');bv,bf=v10.previous.arrays((bo-hb)-bb,root);branch=Mesh(bv,bf,FRAME);branch.save(OUT/'vein-branch-wall.npz');top=topology(network);blood_top=topology(Mesh(*v10.previous.arrays(hb+bb,root),FRAME));outside=[]
    for i,component in enumerate((bb-bo).decompose()):
        v,f=v10.previous.arrays(component,root);nearest=int(np.linalg.norm(centers-v.mean(0),axis=1).argmin());outside.append(dict(index=i,volume_m3=float(component.volume())/1e9,nearestStation=nearest))
    lateral=[r for r in outside if r['nearestStation']<N-2 and r['volume_m3']>1e-23];checks=[]
    contexts=[('fixed-cystic-artery',artery),('GB-innermost',inner)]+[(pid,reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CURRENT))for pid in ['TA98HEPATIC:male:artery','FJ3080','TA98PART:FJ3079:A05.8.01.061','TA98HEPATIC:male:biliary','TA98HEPATIC:male:hepatic_vein']if pid in parts]
    for name,context in contexts:
        contacts,check=pair_contacts(branch,context,0.,max_pairs=500000,max_tests=5000000);assert check['complete_within_pair_budget'];atomic_npz(OUT/(name.replace(':','_')+'-contacts.npz'),**contacts);checks.append(dict(name=name,geometrySha256=context.geometry_sha,**check))
    overlap=(hb^bb).decompose();passed=not lateral and all(c['crossing_face_pairs']==0 and c['coplanar_overlap_face_pairs']==0 for c in checks)and top['closed_oriented_manifold']and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==1 and len(overlap)==1 and len(blood_top['components'])==1 and blood_top['components'][0]['closed_orientable_genus']==0
    proof.update(joinedTopology=top,composedBloodTopology=blood_top,actualParentBloodOverlapComponents=len(overlap),lateralBloodEscapeComponents=lateral,allBloodOutsideComponents=outside,contextChecks=checks,allRequiredGeometryGatesPassed=passed,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'gates.json',proof)
    if passed:
        point=centers[-1];entries=[]
        for name,mesh,color in [('fixed-artery',artery,[.86,.13,.10,1]),('corrected-vein',branch,[.15,.27,.83,1]),('GB-outer',cover,[.35,.58,.16,1]),('GB-inner',inner,[.93,.65,.25,1])]:
            tri=mesh.triangles;take=np.all(tri.max(1)>=point-.012,axis=1)&np.all(tri.min(1)<=point+.012,axis=1)
            if not take.any():continue
            file=OUT/(name+'-tip-crop.npz');Mesh(mesh.vertices,mesh.faces[take],FRAME,source=dict(diagnosticUncappedCrop=True)).save(file);entries.append(dict(mesh=str(file),name=name+' actual source crop',color=color))
        spec=dict(frame=FRAME,meshes=entries,views=[dict(name='separated-vein-and-fixed-artery',direction=[0,0,1]),dict(name='wall-entry-cut',direction=[0,0,1],section=dict(normal=[0,0,1],point=point.tolist(),clear_outer=True))],resolution=900,wall_seconds=40);review=prepare(None,spec);atomic_json(OUT/'review.json',review)
    print(json.dumps(dict(passed=passed,invalidPrismIntervals=invalid,newArteryCrossings=next(c['crossing_face_pairs']for c in checks if c['name']=='fixed-cystic-artery'),innermostCrossings=next(c['crossing_face_pairs']for c in checks if c['name']=='GB-innermost'),lateralEscapes=len(lateral),distalShift_mm=proof['distalEntryShift_m']*1000,elapsedSeconds=time.monotonic()-began)))


if __name__=='__main__':main()
