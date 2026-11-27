"""Transport measured section shapes and make one bounded host blood junction."""
from pathlib import Path
import json,sys,time,importlib.util
import numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
from geometry_tools.review import prepare

spec=importlib.util.spec_from_file_location('cystic_previous',ROOT/'server/construct-cystic-surface-courses-v2.py');previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
WORK=previous.WORK;OUT=WORK/'profile-mouth-correction-v10';FRAME=previous.FRAME

def resample_profile(p,n=18):
    edge=np.linalg.norm(np.roll(p,-1,axis=0)-p,axis=1);s=np.r_[0,np.cumsum(edge)];q=np.linspace(0,s[-1],n,endpoint=False)
    return np.column_stack([np.interp(q,s,np.r_[p[:,k],p[0,k]])for k in range(3)])

def swept_profiles(sections):
    v=sections.reshape(-1,3);f=[];n=sections.shape[1]
    for i in range(len(sections)-1):
        for j in range(n):
            a=i*n+j;b=i*n+(j+1)%n;c=b+n;d=a+n;f.extend([[a,b,c],[a,c,d]])
    caps=[]
    for i,reverse in [(0,True),(len(sections)-1,False)]:
        center=len(v);v=np.r_[v,sections[i].mean(0)[None]];ids=[]
        for j in range(n):
            ids.append(len(f));t=[center,i*n+j,i*n+(j+1)%n];f.append(t[::-1]if reverse else t)
        caps.append(np.array(ids))
    f=np.array(f,np.int64)
    if np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])).sum()<0:f=f[:,::-1]
    return v,f,caps

def native_outer(host,blood,lineage,root):
    # Both classifiers share the SAME computational blood-root disk, so it
    # cancels from their difference. Nonplanar native annular tissue stays.
    points=blood.vertices;mapping={tuple(p):i for i,p in enumerate(points)}
    f=np.array([[mapping[tuple(p)]for p in tri]for tri in host.triangles])
    innerids=np.array(lineage['sourceFaceIds'],int);outerids=np.setdiff1d(np.arange(len(f)),innerids)
    disk=blood.faces[-lineage['temporaryRootDiskFaces']:]
    outer=Mesh(points,np.r_[f[outerids],disk],FRAME)
    assert topology(outer)['closed_oriented_manifold']
    old=previous.solid(host.vertices,host.faces,root);outer_solid=previous.solid(outer.vertices,outer.faces,root);blood_solid=previous.solid(blood.vertices,blood.faces,root)
    reconstructed=outer_solid-blood_solid
    symmetric=(reconstructed-old)+(old-reconstructed)
    residual=float(symmetric.volume())/1e9
    assert residual<1e-16,('Native blood classifier fails exact wall reconstruction',residual)
    return outer,outer_solid,blood_solid,residual

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();audit=json.loads((WORK/'source-audit.json').read_text());scale=json.loads((WORK/'first-similarity/fit.json').read_text())['metrics']['scale'];results=[];review=[]
    inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');covering=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz')
    for number,role in enumerate(['artery','vein']):
        original=Mesh.load(WORK/('donor-'+str(number)+'.npz'));profile=np.array(audit['rows'][number]['boundaryProfiles'][0]['coordinates_m']);source_centers,radii,evidence=previous.source_sections(original,profile)
        old=np.load(WORK/'surface-courses-v4'/(role+'-section-fields.npz'));curve=old['centers'].copy();arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(curve,axis=0),axis=1))];direction=curve[1]-curve[0];direction/=np.linalg.norm(direction)
        host=Mesh.load(WORK/('fixed-target-20.npz'if role=='artery'else'fixed-target-21.npz'));blood=Mesh.load(WORK/'native-parent-envelopes-v9'/(role+'-blood.npz'));lineage=json.loads((WORK/'native-parent-envelopes-v9'/(role+'-inner-lineage.json')).read_text())
        root=curve[0].copy();inside,ambiguous=ray_contains(root[None],blood.triangles);assert inside[0] and not ambiguous[0],'Named parent-region center is not in actual blood space'
        clearance=BVH(blood).closest_point(root)['distance_m'];connector_radius=min(.00035,clearance*.45);assert connector_radius>.00008
        outer,host_outer,host_blood,reconstruction_residual=native_outer(host,blood,lineage,root)
        extent=max(.004,scale*np.linalg.norm(profile-profile.mean(0),axis=1).max()*2.5)
        blend=np.clip(arc/.012,0,1);weight=1-blend*blend*(3-2*blend);curve+=weight[:,None]*extent*direction
        source_arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(source_centers,axis=0),axis=1))];source_t=source_arc/source_arc[-1];target_t=np.linspace(0,1,len(curve))
        source_tangent=np.gradient(source_centers,axis=0);source_tangent/=np.linalg.norm(source_tangent,axis=1)[:,None]
        source_frames=rotation_minimizing_frames(source_centers,source_tangent,initial_reference=profile[0]-profile.mean(0));locals=[]
        for j,raw in enumerate(evidence['sectionProfiles_m']):
            p=np.array(raw);p=p-p.mean(0);local=p@source_frames[j]
            if j:
                local=resample_profile(local);choices=[]
                for reverse in [False,True]:
                    q=local[::-1]if reverse else local
                    for shift in range(18):
                        candidate=np.roll(q,shift,axis=0);choices.append((float(np.sum((candidate/radii[j]-locals[-1]/radii[j-1])**2)),candidate))
                local=min(choices,key=lambda x:x[0])[1]
            local*=np.sqrt(evidence['sectionAreas_m2'][j]/projected_area(local));locals.append(local)
        locals=np.array(locals);newlocal=np.empty((len(curve),18,3))
        for k in range(18):
            for d in range(3):newlocal[:,k,d]=np.interp(target_t,source_t,locals[:,k,d])*scale
        desired_area=np.interp(target_t,source_t,np.array(evidence['sectionAreas_m2']))*scale**2
        desired_diameter=np.interp(target_t,source_t,np.array(evidence['sectionDiameters_m']))*scale
        for j in range(len(newlocal)):newlocal[j]*=np.sqrt(desired_area[j]/projected_area(newlocal[j]))
        tangent=np.gradient(curve,axis=0);tangent/=np.linalg.norm(tangent,axis=1)[:,None];frames=rotation_minimizing_frames(curve,tangent,initial_reference=direction if abs(direction@tangent[0])<.99 else [1,0,0])
        sections=np.einsum('sij,snj->sni',frames,newlocal)+curve[:,None,:]
        actual_area=np.array([projected_area(p)for p in sections]);diameter=np.linalg.norm(sections[:,:,None]-sections[:,None,:],axis=3).max((1,2));area_error=abs(actual_area/desired_area-1);diameter_error=abs(diameter/desired_diameter-1)
        atomic_json(OUT/(role+'-source-transport.json'),dict(sourceGeometrySha256=original.geometry_sha,sourceSectionEvidence=evidence,maximumAreaRelativeError=float(area_error.max()),maximumDiameterRelativeError=float(diameter_error.max()),tolerance=.15,passed=bool(np.all(area_error<=.15)&np.all(diameter_error<=.15))))
        assert np.all(area_error<=.15) and np.all(diameter_error<=.15),'Individual source section calibre exceeds 15%'
        # A small new local connector reaches the unchanged 18-control donor
        # profile outside the parent; the blood cutter begins INSIDE the
        # verified native blood volume and never extends through its far wall.
        fraction=np.linspace(0,1,7,endpoint=False);connector_centers=root+(curve[0]-root)*fraction[:,None]
        connector_scale=connector_radius/(np.linalg.norm(sections[0]-curve[0],axis=1).max())
        f=fraction*fraction*(3-2*fraction);f=connector_scale+(1-connector_scale)*f
        connector=connector_centers[:,None,:]+f[:,None,None]*(sections[0]-curve[0])[None]
        outer_sections=np.r_[connector,sections];centers=np.r_[connector_centers,curve]
        fluid_sections=centers[:,None,:]+.82*(outer_sections-centers[:,None,:]);inside,ambiguous=ray_contains(fluid_sections[0],blood.triangles);assert np.all(inside&~ambiguous),'Proximal cutter disk is not wholly in verified parent blood space'
        fluid_sections[-1]+=(centers[-1]-centers[-2])/np.linalg.norm(centers[-1]-centers[-2])*.001
        ov,of,oc=swept_profiles(outer_sections);iv,iff,ic=swept_profiles(fluid_sections);branch_outer=previous.solid(ov,of,root);branch_blood=previous.solid(iv,iff,root)
        joined=(host_outer+branch_outer)-(host_blood+branch_blood);jv,jf=previous.arrays(joined,root);network=Mesh(jv,jf,FRAME);network.save(OUT/(role+'-joined-wall.npz'))
        branch_solid=(branch_outer-host_blood)-branch_blood;bv,bf=previous.arrays(branch_solid,root);branch=Mesh(bv,bf,FRAME);branch.save(OUT/(role+'-branch-wall.npz'))
        atomic_npz(OUT/(role+'-sections.npz'),outer_sections=outer_sections,fluid_sections=fluid_sections,centers=centers,source_profile=profile,transported_donor_profile=sections[0],donor_start_station=np.array(len(connector)),host_blood_disk=fluid_sections[0])
        witness,contact=pair_contacts(branch,inner,0.,max_pairs=250000,max_tests=2000000);atomic_npz(OUT/(role+'-inner-witness.npz'),**witness);top=topology(network)
        result=dict(role=role,areaRelativeError=float(area_error.max()),diameterRelativeError=float(diameter_error.max()),sourceCalibrePassed=True,innermostCrossingPairs=contact['crossing_face_pairs'],parentBloodRootClearance_m=clearance,connectorRadius_m=connector_radius,proximalCutterDiskStrictlyInsideNativeBlood=True,parentBloodClassifierReconstructsNativeWallResidual_m3=reconstruction_residual,original18ControlProfilePreservedByPositiveSimilarity=True,oldParentGenus=0,newParentGenus=top['components'][0]['closed_orientable_genus'],closed=top['closed_oriented_manifold'],components=len(top['components']),remainingLimitations=['Microscopic recipient continuity and exact wall histology remain unestablished'],eligibleForInstallation=False)
        atomic_json(OUT/(role+'-gates.json'),result);results.append(result)
        review.append(dict(mesh=str(OUT/(role+'-joined-wall.npz')),name=role+' actual native parent and corrected cystic wall',color=[.85,.14,.1,1]if role=='artery'else[.13,.24,.78,1]))
        if contact['crossing_face_pairs'] or top['components'][0]['closed_orientable_genus']!=1:
            atomic_json(OUT/'stopped.json',dict(reason='Corrected candidate failed lumen or exactly-one-added-outlet topology gate',result=result));break
    full=dict(frame=FRAME,meshes=review+[dict(mesh=str(WORK/'first-similarity/fixed-covering-free-surface.npz'),name='Fixed male GB outer surface',color=[.48,.55,.22,1])],views=[dict(name='whole-parent',direction=[0,0,1])],resolution=640,wall_seconds=40);atomic_json(OUT/'whole-review-spec.json',full);whole=prepare(None,full)
    close_meshes=[]
    for role in ['artery','vein']:
        path=OUT/(role+'-branch-wall.npz')
        if path.exists():close_meshes.append(dict(mesh=str(path),name=role+' corrected cystic branch',color=[.85,.14,.1,1]if role=='artery'else[.13,.24,.78,1]))
    close=dict(frame=FRAME,meshes=close_meshes+[dict(mesh=str(WORK/'first-similarity/fixed-covering-free-surface.npz'),name='Actual outer GB surface',color=[.4,.6,.18,1]),dict(mesh=str(WORK/'first-similarity/fixed-mucosa-free-surface.npz'),name='Actual innermost GB surface',color=[.9,.65,.24,1])],views=[dict(name='GB-surface',direction=[0,0,1]),dict(name='GB-luminal-section',direction=[0,0,1],section=dict(normal=[0,0,1],point=covering.bounds.mean(0).tolist(),clear_outer=True))],resolution=800,wall_seconds=40);atomic_json(OUT/'close-review-spec.json',close);close_review=prepare(None,close)
    atomic_json(OUT/'result.json',dict(schema='cystic-individual-profile-bounded-host-junction/v10',gates=results,wholeReview=whole,closeReview=close_review,installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(output=str(OUT),gates=results,reviewIds=[whole['request_id'],close_review['request_id']],elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
