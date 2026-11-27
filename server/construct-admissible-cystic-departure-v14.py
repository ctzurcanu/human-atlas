"""One deterministic full-connector transport from the feasible donor port."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
from geometry_tools.review import prepare
sp=importlib.util.spec_from_file_location('v12',Path(__file__).with_name('round-cystic-station13-v12.py'));v12=importlib.util.module_from_spec(sp);sp.loader.exec_module(v12);v10=v12.v10
WORK=v10.WORK;OUT=WORK/'admissible-departure-v14';FRAME=v10.FRAME

def frame(profile,forward):
    q=profile-profile.mean(0);normal=np.cross(q,np.roll(q,-1,axis=0)).sum(0);normal/=np.linalg.norm(normal)
    if normal@forward<0:normal=-normal
    x=q[0]-normal*(q[0]@normal);x/=np.linalg.norm(x);return np.column_stack([x,np.cross(normal,x),normal])

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();feasibility=json.loads((WORK/'departure-feasibility-v13/feasibility.json').read_text());assert feasibility['feasibleLocalPrerequisites']
    z=np.load(WORK/'profile-mouth-correction-v10/artery-sections.npz');old=z['outer_sections'];old_centers=z['centers'];root=old_centers[0];outer=old.copy();centers=old_centers.copy();fluid=z['fluid_sections'].copy();port,last=7,20
    direction=np.array(feasibility['admissibleProposedWorldDirection']);rotation=np.array(feasibility['wholeProfileProperRotation']);distance=feasibility['externalProfilePlacementTrials'][-1]['distanceFromBloodRoot_m'];port_center=root+distance*direction;literal_profile=(old[port]-old_centers[port])@rotation.T+port_center
    target_end=frame(old[last],old_centers[last+1]-old_centers[last-1]);start_frame=frame(literal_profile,direction);length=np.linalg.norm(old_centers[last]-port_center)
    path=CubicHermiteSpline([0,1],np.array([port_center,old_centers[last]]),np.array([direction,target_end[:,2]])*length)
    source_arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(old_centers[port:last+1],axis=0),axis=1))];t=source_arc/source_arc[-1];new_centers=path(t);tangent=path(t,1);tangent/=np.linalg.norm(tangent,axis=1)[:,None];frames=rotation_minimizing_frames(new_centers,tangent,initial_reference=start_frame[:,0])
    angle=np.arctan2(target_end[:,0]@frames[-1,:,1],target_end[:,0]@frames[-1,:,0])
    for j in range(port,last+1):
        i=j-port;basis=frames[i];theta=t[i]*angle;c,s=np.cos(theta),np.sin(theta);roll=np.array([[c,-s,0],[s,c,0],[0,0,1]]);target=basis@roll;source=frame(old[j],old_centers[min(j+1,len(old)-1)]-old_centers[max(j-1,0)])
        centers[j]=new_centers[i];outer[j]=(old[j]-old_centers[j])@source@target.T+centers[j]
    outer[port]=literal_profile;centers[port]=port_center;outer[last]=old[last];centers[last]=old_centers[last]
    initial_radius=np.linalg.norm(old[0]-root,axis=1).max();full_radius=np.linalg.norm(literal_profile-port_center,axis=1).max();initial_scale=initial_radius/full_radius
    for j in range(port):
        fraction=j/port;blend=fraction*fraction*(3-2*fraction);factor=initial_scale+(1-initial_scale)*blend;centers[j]=root+fraction*distance*direction;outer[j]=centers[j]+factor*(literal_profile-port_center)
    fluid[:last]=centers[:last,None]+.82*(outer[:last]-centers[:last,None]);assert np.array_equal(outer[last:],old[last:]);assert np.array_equal(fluid[last:],z['fluid_sections'][last:]);assert np.array_equal(centers[0],root)
    area=lambda a:np.array([projected_area(p)for p in a]);diameter=lambda a:np.linalg.norm(a[:,:,None]-a[:,None,:],axis=3).max((1,2));area_error=float(abs(area(outer)/area(old)-1).max());diameter_error=float(abs(diameter(outer)/diameter(old)-1).max())
    ov,of,oc=v10.swept_profiles(outer);iv,iff,ic=v10.swept_profiles(fluid);side=of[:-36];radial=ov[side].mean(1)-np.repeat(centers,36,axis=0)[:len(side)];fn=np.cross(ov[side[:,1]]-ov[side[:,0]],ov[side[:,2]]-ov[side[:,0]])
    if np.median(np.sum(radial*fn,axis=1))<0:side=side[:,::-1]
    nominal_inner=centers[:,None]+.82*(outer-centers[:,None]);prisms=StationPrismConstraints(np.array([outer.reshape(-1,3),nominal_inner.reshape(-1,3)]),side,np.repeat(np.arange(len(outer)),18),np.ones(len(outer),bool),allow_source_faults=True);invalid_faces=np.flatnonzero(~prisms.sourceHealthy[0]);invalid=np.unique(np.repeat(np.arange(len(outer)),18)[side[invalid_faces]].min(1)).tolist()
    atomic_npz(OUT/'sections.npz',outer_sections=outer,fluid_sections=fluid,centers=centers,source_profile=z['source_profile'],fixed_distal_start=np.array(last),moved_intermediate_port_station=np.array(port),original_profile_new_pose=literal_profile)
    # Preserve failed fields and stop immediately instead of performing CSG
    # on a continuously inverted paired wall.
    base=dict(schema='cystic-admissible-departure/v14',departureFeasibilitySha256=digest_file(WORK/'departure-feasibility-v13/feasibility.json'),continuousPrismInvalidIntervals=invalid,continuousPrismPassed=not invalid,maximumContinuousPrismDeterminant=float(prisms.sourceMaximum.max()),original18ControlProfileShapePreserved=True,originalProfileWorldPoseMovedAsConstructedIntermediate=True,fixedBloodRootLiteral=True,fixedDistalSectionsLiteral=True,maximumExternalAreaRelativeChange=area_error,maximumExternalDiameterRelativeChange=diameter_error,sourceDiameterErrorUpperBound=.03756634960524918+diameter_error,sourceCalibreTolerance=.15,installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began)
    if invalid or area_error>.15 or diameter_error>.15:
        base['stoppedBeforeCSG']=True;atomic_json(OUT/'gates.json',base);print(json.dumps(dict(output=str(OUT),stoppedBeforeCSG=True,invalidIntervals=invalid,areaRelativeChange=area_error,diameterRelativeChange=diameter_error,elapsedSeconds=base['elapsedSeconds'])));return
    host=Mesh.load(WORK/'fixed-target-20.npz');blood=Mesh.load(WORK/'native-parent-envelopes-v9/artery-blood.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/artery-inner-lineage.json').read_text());native_outer,host_outer,host_blood,residual=v10.native_outer(host,blood,lineage,root);inside,ambiguous=ray_contains(fluid[0],blood.triangles);assert np.all(inside&~ambiguous)
    branch_outer=v10.previous.solid(ov,of,root);branch_blood=v10.previous.solid(iv,iff,root);joined=(host_outer+branch_outer)-(host_blood+branch_blood);jv,jf=v10.previous.arrays(joined,root);network=Mesh(jv,jf,FRAME);network.save(OUT/'artery-joined-wall.npz');bv,bf=v10.previous.arrays((branch_outer-host_blood)-branch_blood,root);branch=Mesh(bv,bf,FRAME);branch.save(OUT/'artery-branch-wall.npz')
    inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');cover=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz');witness,contact=pair_contacts(branch,inner,0.,max_pairs=250000,max_tests=2000000);atomic_npz(OUT/'innermost-witness.npz',**witness);tree=cKDTree(centers);outside=[]
    for i,solid in enumerate((branch_blood-branch_outer).decompose()):
        v,f=v10.previous.arrays(solid,root);mesh=Mesh(v,f,FRAME);mesh.save(OUT/('blood-outside-'+str(i)+'.npz'));point=mesh.vertices.mean(0);outside.append(dict(index=i,volume_m3=float(solid.volume())/1e9,nearestStation=int(tree.query(point)[1]),center_m=point.tolist()))
    lateral=[r for r in outside if r['nearestStation']<len(centers)-2 and r['volume_m3']>1e-23];top=topology(network);overlap=(host_blood^branch_blood).decompose();blood_top=topology(Mesh(*v10.previous.arrays(host_blood+branch_blood,root),FRAME))
    passed=not lateral and contact['crossing_face_pairs']==0 and top['closed_oriented_manifold'] and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==1 and len(overlap)==1 and blood_top['components'][0]['closed_orientable_genus']==0
    base.update(lateralBloodEscapeComponents=lateral,allBloodOutsideOuterComponents=outside,innermostCrossingPairs=contact['crossing_face_pairs'],joinedTopology=top,hostDonorBloodOverlapComponents=len(overlap),composedBloodTopology=blood_top,proximalBloodDiskInsideActualParent=True,actualParentBloodCommunication=True,allRequiredGeometryGatesPassed=passed,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'gates.json',base)
    reviews=[]
    if passed:
        for name,point,half in [('tip',centers[-1],.008),('junction',port_center,.009)]:
            meshes=[]
            for label,mesh,color in [('arterial-wall',network,[.82,.13,.1,1]),('GB-outer',cover,[.35,.58,.16,1]),('GB-inner',inner,[.93,.65,.25,1])]:
                tri=mesh.triangles;take=np.all(tri.max(1)>=point-half,axis=1)&np.all(tri.min(1)<=point+half,axis=1)
                if not np.any(take):continue
                cropped=Mesh(mesh.vertices,mesh.faces[take],FRAME,source=dict(diagnosticUncappedCrop=True));file=OUT/(name+'-'+label+'.npz');cropped.save(file);meshes.append(dict(mesh=str(file),name=label+' actual source crop',color=color))
            review_spec=dict(frame=FRAME,meshes=meshes,views=[dict(name=name+'-surface',direction=[0,0,1]),dict(name=name+'-cut',direction=[0,0,1],section=dict(normal=[0,0,1],point=point.tolist(),clear_outer=True))],resolution=800,wall_seconds=40);atomic_json(OUT/(name+'-review-spec.json'),review_spec);review=prepare(None,review_spec);atomic_json(OUT/(name+'-review.json'),review);reviews.append(review['request_id'])
    print(json.dumps(dict(output=str(OUT),continuousPrismPassed=True,lateralBloodEscapeComponents=len(lateral),GBLumenCrossings=contact['crossing_face_pairs'],joinedWallGenus=[c['closed_orientable_genus']for c in top['components']],actualBloodOverlapComponents=len(overlap),passed=passed,reviewIds=reviews,elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
