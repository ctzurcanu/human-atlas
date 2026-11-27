"""One local paired-section correction between fixed source-profile stations."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.review import prepare
sp=importlib.util.spec_from_file_location('v10',Path(__file__).with_name('correct-cystic-profiles-and-host-mouths-v10.py'));v10=importlib.util.module_from_spec(sp);sp.loader.exec_module(v10)
WORK=v10.WORK;OUT=WORK/'station13-rounded-v12';FRAME=v10.FRAME

def rotate_normal(normal,target):
    cross=np.cross(normal,target);s=np.linalg.norm(cross);c=float(normal@target)
    if s<1e-12:return np.eye(3)
    cross/=s;K=np.array([[0,-cross[2],cross[1]],[cross[2],0,-cross[0]],[-cross[1],cross[0],0]])
    return np.eye(3)+K*s+K@K*(1-c)

def crop(mesh,center,half,name):
    center=np.array(center);tri=mesh.triangles;take=np.all(tri.max(1)>=center-half,axis=1)&np.all(tri.min(1)<=center+half,axis=1);assert np.any(take)
    result=Mesh(mesh.vertices,mesh.faces[take],FRAME,source=dict(diagnosticCropOnly=True,uncappedCrop=True,sourceGeometrySha256=mesh.geometry_sha));path=OUT/(name+'.npz');result.save(path);return path

def main():
    start=time.monotonic();assert not OUT.exists();OUT.mkdir();z=np.load(WORK/'profile-mouth-correction-v10/artery-sections.npz');old=z['outer_sections'];old_centers=z['centers'];outer=old.copy();centers=old_centers.copy();fluid=z['fluid_sections'].copy();low,high=7,20
    normal=np.linalg.svd(old[low]-old_centers[low],full_matrices=False)[2][-1];travel=old_centers[low+1]-old_centers[low]
    if normal@travel<0:normal=-normal
    end_direction=old_centers[high+1]-old_centers[high-1];end_direction/=np.linalg.norm(end_direction);length=np.linalg.norm(np.diff(old_centers[low:high+1],axis=0),axis=1).sum()
    path=CubicHermiteSpline([0,1],old_centers[[low,high]],np.array([normal,end_direction])*length)
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(old_centers[low:high+1],axis=0),axis=1))];t=s/s[-1];new_centers=path(t);tangent=path(t,1);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    for k in range(low+1,high):
        local=old[k]-old_centers[k];basis=np.linalg.svd(local,full_matrices=False)[2];n=basis[-1]
        if n@(old_centers[k+1]-old_centers[k-1])<0:n=-n
        rotation=rotate_normal(n,tangent[k-low]);centers[k]=new_centers[k-low];outer[k]=local@rotation.T+centers[k];fluid[k]=(z['fluid_sections'][k]-old_centers[k])@rotation.T+centers[k]
    protected=np.r_[np.arange(low+1),np.arange(high,len(outer))];assert np.array_equal(outer[protected],old[protected]);assert np.array_equal(fluid[protected],z['fluid_sections'][protected])
    area=lambda a:np.array([projected_area(p)for p in a]);diameter=lambda a:np.linalg.norm(a[:,:,None]-a[:,None,:],axis=3).max((1,2))
    area_error=float(abs(area(outer)/area(old)-1).max());diameter_error=float(abs(diameter(outer)/diameter(old)-1).max());assert area_error<1e-10 and diameter_error<1e-10
    nominal_inner=centers[:,None]+.82*(outer-centers[:,None]);ov,of,oc=v10.swept_profiles(outer);iv,iff,ic=v10.swept_profiles(fluid);side=of[:-36];radial=ov[side].mean(1)-np.repeat(centers,36,axis=0)[:len(side)];normal_faces=np.cross(ov[side[:,1]]-ov[side[:,0]],ov[side[:,2]]-ov[side[:,0]])
    if np.median(np.sum(radial*normal_faces,axis=1))<0:side=side[:,::-1]
    prism=StationPrismConstraints(np.array([outer.reshape(-1,3),nominal_inner.reshape(-1,3)]),side,np.repeat(np.arange(len(outer)),18),np.ones(len(outer),bool),allow_source_faults=True)
    bad=np.flatnonzero(~prism.sourceHealthy[0]);invalid=np.unique(np.repeat(np.arange(len(outer)),18)[side[bad]].min(1)).tolist()
    atomic_npz(OUT/'sections.npz',outer_sections=outer,fluid_sections=fluid,centers=centers,protected_station_ids=protected,source_profile=z['source_profile'],transported_donor_profile=outer[low])
    root=centers[0];host=Mesh.load(WORK/'fixed-target-20.npz');blood=Mesh.load(WORK/'native-parent-envelopes-v9/artery-blood.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/artery-inner-lineage.json').read_text());native_outer,host_outer,host_blood,residual=v10.native_outer(host,blood,lineage,root)
    branch_outer=v10.previous.solid(ov,of,root);branch_blood=v10.previous.solid(iv,iff,root);joined=(host_outer+branch_outer)-(host_blood+branch_blood);jv,jf=v10.previous.arrays(joined,root);network=Mesh(jv,jf,FRAME);network.save(OUT/'artery-joined-wall.npz')
    bv,bf=v10.previous.arrays((branch_outer-host_blood)-branch_blood,root);branch=Mesh(bv,bf,FRAME);branch.save(OUT/'artery-branch-wall.npz');inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');cover=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz');witness,contact=pair_contacts(branch,inner,0.,max_pairs=250000,max_tests=2000000);atomic_npz(OUT/'innermost-witness.npz',**witness)
    tree=cKDTree(centers);outside=[]
    for i,solid in enumerate((branch_blood-branch_outer).decompose()):
        v,f=v10.previous.arrays(solid,root);mesh=Mesh(v,f,FRAME);mesh.save(OUT/('blood-outside-'+str(i)+'.npz'));point=mesh.vertices.mean(0);outside.append(dict(index=i,volume_m3=float(solid.volume())/1e9,nearestStation=int(tree.query(point)[1]),center_m=point.tolist(),bounds_m=mesh.bounds.tolist()))
    top=topology(network);lateral=[r for r in outside if r['nearestStation']<len(centers)-2 and r['volume_m3']>1e-23]
    result=dict(schema='cystic-station13-local-rounded-transport/v12',movedStationIds=list(range(low+1,high)),fixedStationIds=protected.tolist(),fixedOriginal18ControlPortLiteral=bool(np.array_equal(outer[7],old[7])),fixedParentCutterAndTerminalLiteral=True,maximumExternalAreaRelativeChange=area_error,maximumExternalDiameterRelativeChange=diameter_error,maximumSourceDiameterError=.03756634960524918,sourceDiameterTolerance=.15,continuousPrismInvalidIntervals=invalid,continuousPrismPassed=not invalid,maximumContinuousPrismDeterminant=float(prism.sourceMaximum.max()),lateralBloodEscapeComponents=lateral,allBloodOutsideOuterComponents=outside,innermostCrossingPairs=contact['crossing_face_pairs'],joinedTopology=top,parentBloodClassifierReconstructionResidual_m3=residual,geometryAcceptance=False,installed=False,eligibleForInstallation=False,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'gates.json',result)
    reviews=[]
    for name,point,half in [('GB-tip',centers[-1],.008),('parent-junction',centers[7],.009),('side-opening-site',old_centers[13],.006)]:
        rows=[]
        for label,mesh,color in [('arterial-wall',network,[.82,.13,.1,1]),('GB-outer',cover,[.35,.58,.16,1]),('GB-inner',inner,[.93,.65,.25,1])]:
            tri=mesh.triangles
            if not np.any(np.all(tri.max(1)>=point-half,axis=1)&np.all(tri.min(1)<=point+half,axis=1)):continue
            path=crop(mesh,point,half,name+'-'+label);rows.append(dict(mesh=str(path),name=label+' actual uncapped local crop',color=color))
        spec=dict(frame=FRAME,meshes=rows,views=[dict(name=name+'-surface',direction=[0,0,1]),dict(name=name+'-section',direction=[0,0,1],section=dict(normal=[1,0,0]if name=='parent-junction'else[0,0,1],point=point.tolist(),clear_outer=True))],resolution=800,wall_seconds=40);atomic_json(OUT/(name+'-review-spec.json'),spec);review=prepare(None,spec);atomic_json(OUT/(name+'-review.json'),review);reviews.append(review['request_id'])
    print(json.dumps(dict(output=str(OUT),continuousPrismInvalidIntervals=invalid,lateralBloodEscapeComponents=len(lateral),innermostCrossingPairs=contact['crossing_face_pairs'],wallGenus=[c['closed_orientable_genus']for c in top['components']],areaRelativeChange=area_error,diameterRelativeChange=diameter_error,reviewIds=reviews,elapsedSeconds=time.monotonic()-start)))

if __name__=='__main__':main()
