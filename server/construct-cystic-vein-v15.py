"""One source-specific cystic-vein course and actual right-portal junction."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.audit import topology
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
from geometry_tools.review import prepare
sp=importlib.util.spec_from_file_location('v14',Path(__file__).with_name('construct-admissible-cystic-departure-v14.py'));v14=importlib.util.module_from_spec(sp);sp.loader.exec_module(v14);v10=v14.v10
rp=importlib.util.spec_from_file_location('feasibility',Path(__file__).with_name('cystic-departure-feasibility-v13.py'));feas=importlib.util.module_from_spec(rp);rp.loader.exec_module(feas)
WORK=v10.WORK;OUT=WORK/'cystic-vein-v15';FRAME=v10.FRAME

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();audit=json.loads((WORK/'source-audit.json').read_text());profile=np.array(audit['rows'][1]['boundaryProfiles'][0]['coordinates_m']);donor=Mesh.load(WORK/'donor-1.npz');source_centers,radii,evidence=v10.previous.source_sections(donor,profile);scale=json.loads((WORK/'first-similarity/fit.json').read_text())['metrics']['scale']
    old=np.load(WORK/'surface-courses-v4/vein-section-fields.npz');curve=old['centers'].copy();root=curve[0].copy();host=Mesh.load(WORK/'fixed-target-21.npz');blood=Mesh.load(WORK/'native-parent-envelopes-v9/vein-blood.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/vein-inner-lineage.json').read_text());native_outer,host_outer,host_blood,residual=v10.native_outer(host,blood,lineage,root)
    inside,ambiguous=ray_contains(root[None],blood.triangles);assert inside[0] and not ambiguous[0];clearance=BVH(blood).closest_point(root)['distance_m'];connector_radius=min(.00035,clearance*.45);assert connector_radius>.00008
    direction=curve[24]-root;direction/=np.linalg.norm(direction);source_first_frame=v14.frame(profile,source_centers[1]-source_centers[0]);source_n=source_first_frame[:,2];cross=np.cross(source_n,direction);s=np.linalg.norm(cross);c=source_n@direction
    if s<1e-12:rotation=np.eye(3)
    else:
        axis=cross/s;K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]);rotation=np.eye(3)+K*s+K@K*(1-c)
    local_profile=(profile-profile.mean(0))*scale;posed_local=local_profile@rotation.T;blood_exit=feas.first_ray_exit(blood,root,direction);outer_exit=feas.first_ray_exit(native_outer,root,direction)
    # One conservative placement, not a sequence of displacement trials.
    axial_support=float(abs(posed_local@direction).max());distance=outer_exit['distance_m']+axial_support+.0015;port_center=root+distance*direction;posed=posed_local+port_center;inside,ambiguous=ray_contains(posed,native_outer.triangles);near,_=BVH(native_outer).closest_points(posed);placement_clearance=float(near['distance_m'].min())
    assert np.all(~inside&~ambiguous) and placement_clearance>.0001,'Single analytically placed original port has insufficient native-parent clearance'
    fixed=13;end=curve[fixed];end_tangent=curve[fixed+1]-curve[fixed-1];end_tangent/=np.linalg.norm(end_tangent);length=np.linalg.norm(end-port_center);path=CubicHermiteSpline([0,1],np.array([port_center,end]),np.array([direction,end_tangent])*length)
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(curve[:fixed+1],axis=0),axis=1))];parameter=arc/arc[-1];curve[:fixed+1]=path(parameter);tangent=np.gradient(curve,axis=0);tangent[:fixed+1]=path(parameter,1);tangent/=np.linalg.norm(tangent,axis=1)[:,None];frames=rotation_minimizing_frames(curve,tangent,initial_reference=posed_local[0])
    source_arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(source_centers,axis=0),axis=1))];source_t=source_arc/source_arc[-1];target_t=np.linspace(0,1,len(curve));local_sections=[]
    for j,raw in enumerate(evidence['sectionProfiles_m']):
        p=np.array(raw);forward=source_centers[min(j+1,len(source_centers)-1)]-source_centers[max(j-1,0)];basis=v14.frame(p,forward);q=(p-p.mean(0))@basis
        if j:
            q=v10.resample_profile(q);choices=[]
            for reverse in [False,True]:
                qq=q[::-1]if reverse else q
                for shift in range(18):
                    candidate=np.roll(qq,shift,axis=0);choices.append((float(np.sum((candidate/radii[j]-local_sections[-1]/radii[j-1])**2)),candidate))
            q=min(choices,key=lambda x:x[0])[1]
        q*=np.sqrt(evidence['sectionAreas_m2'][j]/projected_area(q));local_sections.append(q)
    local_sections=np.array(local_sections);newlocal=np.empty((len(curve),18,3))
    for k in range(18):
        for d in range(3):newlocal[:,k,d]=np.interp(target_t,source_t,local_sections[:,k,d])*scale
    expected_area=np.interp(target_t,source_t,evidence['sectionAreas_m2'])*scale**2;expected_diameter=np.interp(target_t,source_t,evidence['sectionDiameters_m'])*scale
    for j in range(len(newlocal)):newlocal[j]*=np.sqrt(expected_area[j]/projected_area(newlocal[j]))
    sections=np.einsum('sij,snj->sni',frames,newlocal)+curve[:,None,:];sections[0]=posed
    area=np.array([projected_area(p)for p in sections]);diameter=np.linalg.norm(sections[:,:,None]-sections[:,None,:],axis=3).max((1,2));area_error=float(abs(area/expected_area-1).max());diameter_error=float(abs(diameter/expected_diameter-1).max());assert area_error<=.15 and diameter_error<=.15
    fraction=np.linspace(0,1,7,endpoint=False);connector_centers=root+fraction[:,None]*(port_center-root);factor=connector_radius/np.linalg.norm(posed_local,axis=1).max();blend=fraction*fraction*(3-2*fraction);f=factor+(1-factor)*blend;connector=connector_centers[:,None]+f[:,None,None]*posed_local
    outer=np.r_[connector,sections];centers=np.r_[connector_centers,curve];fluid=centers[:,None]+.82*(outer-centers[:,None]);inside,ambiguous=ray_contains(fluid[0],blood.triangles);assert np.all(inside&~ambiguous);direction_end=centers[-1]-centers[-2];direction_end/=np.linalg.norm(direction_end);fluid[-1]+=direction_end*.001
    ov,of,oc=v10.swept_profiles(outer);iv,iff,ic=v10.swept_profiles(fluid);side=of[:-36];radial=ov[side].mean(1)-np.repeat(centers,36,axis=0)[:len(side)];fn=np.cross(ov[side[:,1]]-ov[side[:,0]],ov[side[:,2]]-ov[side[:,0]])
    if np.median(np.sum(radial*fn,axis=1))<0:side=side[:,::-1]
    nominal_inner=centers[:,None]+.82*(outer-centers[:,None]);prism=StationPrismConstraints(np.array([outer.reshape(-1,3),nominal_inner.reshape(-1,3)]),side,np.repeat(np.arange(len(outer)),18),np.ones(len(outer),bool),allow_source_faults=True);bad=np.flatnonzero(~prism.sourceHealthy[0]);invalid=np.unique(np.repeat(np.arange(len(outer)),18)[side[bad]].min(1)).tolist()
    atomic_npz(OUT/'sections.npz',outer_sections=outer,fluid_sections=fluid,centers=centers,original_source_profile=profile,original_profile_new_pose=posed,source_section_local_shapes=local_sections,source_centers=source_centers)
    base=dict(schema='cystic-vein-individual-source-transport/v15',sourceGeometrySha256=donor.geometry_sha,sourceSectionEvidence=evidence,sourceScale=scale,maximumAreaRelativeError=area_error,maximumDiameterRelativeError=diameter_error,calibreTolerance=.15,continuousPrismInvalidIntervals=invalid,maximumContinuousPrismDeterminant=float(prism.sourceMaximum.max()),continuousPrismPassed=not invalid,original18ControlProfileProperPosePreserved=True,nativeParentRootBloodClearance_m=clearance,originalPortNativeParentClearance_m=placement_clearance,proximalBloodDiskInsideNativeParent=True,parentBloodClassifierReconstructionResidual_m3=residual,organsAndSkeletonUnchanged=True,installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began)
    if invalid:
        base['stoppedBeforeCSG']=True;atomic_json(OUT/'gates.json',base);print(json.dumps(dict(output=str(OUT),invalidPrismIntervals=invalid,sourceDiameterError=diameter_error,stoppedBeforeCSG=True,elapsedSeconds=base['elapsedSeconds'])));return
    branch_outer=v10.previous.solid(ov,of,root);branch_blood=v10.previous.solid(iv,iff,root);joined=(host_outer+branch_outer)-(host_blood+branch_blood);jv,jf=v10.previous.arrays(joined,root);network=Mesh(jv,jf,FRAME);network.save(OUT/'vein-joined-wall.npz');bv,bf=v10.previous.arrays((branch_outer-host_blood)-branch_blood,root);branch=Mesh(bv,bf,FRAME);branch.save(OUT/'vein-branch-wall.npz')
    inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');cover=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz');witness,contacts=pair_contacts(branch,inner,0.,max_pairs=250000,max_tests=2000000);atomic_npz(OUT/'innermost-witness.npz',**witness);tree=cKDTree(centers);outside=[]
    for i,solid in enumerate((branch_blood-branch_outer).decompose()):
        v,f=v10.previous.arrays(solid,root);mesh=Mesh(v,f,FRAME);mesh.save(OUT/('blood-outside-'+str(i)+'.npz'));point=mesh.vertices.mean(0);outside.append(dict(index=i,volume_m3=float(solid.volume())/1e9,nearestStation=int(tree.query(point)[1]),center_m=point.tolist()))
    lateral=[r for r in outside if r['nearestStation']<len(centers)-2 and r['volume_m3']>1e-23];top=topology(network);overlap=(host_blood^branch_blood).decompose();blood_top=topology(Mesh(*v10.previous.arrays(host_blood+branch_blood,root),FRAME));_,wall_contact=pair_contacts(branch,cover,0.,max_pairs=250000,max_tests=2000000)
    passed=not lateral and not contacts['crossing_face_pairs'] and top['closed_oriented_manifold'] and len(top['components'])==1 and top['components'][0]['closed_orientable_genus']==1 and len(overlap)==1 and blood_top['components'][0]['closed_orientable_genus']==0
    base.update(innermostCrossingPairs=contacts['crossing_face_pairs'],recipientOuterWallContactPairs=wall_contact['crossing_face_pairs'],lateralBloodEscapeComponents=lateral,allBloodOutsideOuterComponents=outside,joinedTopology=top,hostDonorBloodOverlapComponents=len(overlap),composedBloodTopology=blood_top,allRequiredGeometryGatesPassed=passed,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'gates.json',base)
    reviews=[]
    if passed:
        for name,point,half in [('vein-tip',centers[-1],.008),('portal-junction',port_center,.009)]:
            entries=[]
            for label,mesh,color in [('venous-wall',network,[.13,.25,.78,1]),('GB-outer',cover,[.35,.58,.16,1]),('GB-inner',inner,[.93,.65,.25,1])]:
                tri=mesh.triangles;selected=np.all(tri.max(1)>=point-half,axis=1)&np.all(tri.min(1)<=point+half,axis=1)
                if not np.any(selected):continue
                cropped=Mesh(mesh.vertices,mesh.faces[selected],FRAME,source=dict(diagnosticUncappedCrop=True));file=OUT/(name+'-'+label+'.npz');cropped.save(file);entries.append(dict(mesh=str(file),name=label+' actual crop',color=color))
            review_spec=dict(frame=FRAME,meshes=entries,views=[dict(name=name+'-surface',direction=[0,0,1]),dict(name=name+'-cut',direction=[0,0,1],section=dict(normal=[0,0,1],point=point.tolist(),clear_outer=True))],resolution=800,wall_seconds=40);atomic_json(OUT/(name+'-review-spec.json'),review_spec);review=prepare(None,review_spec);atomic_json(OUT/(name+'-review.json'),review);reviews.append(review['request_id'])
    print(json.dumps(dict(output=str(OUT),continuousPrismsPassed=True,sourceDiameterError=diameter_error,lateralEscapes=len(lateral),GBLumenCrossings=contacts['crossing_face_pairs'],recipientWallContacts=wall_contact['crossing_face_pairs'],wallGenus=[c['closed_orientable_genus']for c in top['components']],actualBloodOverlapComponents=len(overlap),passed=passed,reviewIds=reviews,elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
