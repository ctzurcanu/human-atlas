"""Analytic local departure/prism and actual-parent clearance; no new mesh."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_prism_constraints import StationPrismConstraints
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
sp=importlib.util.spec_from_file_location('v10',Path(__file__).with_name('correct-cystic-profiles-and-host-mouths-v10.py'));v10=importlib.util.module_from_spec(sp);sp.loader.exec_module(v10)
WORK=v10.WORK;OUT=WORK/'departure-feasibility-v13'

def first_ray_exit(mesh,p,direction):
    tri=mesh.triangles-p;e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];h=np.cross(direction,e2);den=np.einsum('ij,ij->i',e1,h);valid=abs(den)>1e-20
    inverse=np.divide(1,den,out=np.zeros_like(den),where=valid);s=-tri[:,0];u=np.einsum('ij,ij->i',s,h)*inverse;q=np.cross(s,e1);v=(q@direction)*inverse;t=np.einsum('ij,ij->i',e2,q)*inverse
    selected=valid&(u>=-1e-12)&(v>=-1e-12)&(u+v<=1+1e-12)&(t>1e-10);ids=np.flatnonzero(selected);assert len(ids)
    i=int(ids[np.argmin(t[ids])]);return dict(distance_m=float(t[i]),point_m=(p+t[i]*direction).tolist(),sourceTriangle=i)

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();z=np.load(WORK/'profile-mouth-correction-v10/artery-sections.npz');ring=z['outer_sections'][7];center=ring.mean(0);local=ring-center;area_normal=np.cross(local,np.roll(local,-1,axis=0)).sum(0);normal=area_normal/np.linalg.norm(area_normal);forward=z['centers'][20]-center
    if normal@forward<0:normal=-normal
    length=.001;sections=np.array([ring,ring+length*normal]);inner=sections.mean(1)[:,None]+.82*(sections-sections.mean(1)[:,None]);v,f,c=v10.swept_profiles(sections);faces=f[:-36]
    radial=v[faces].mean(1)-np.repeat(sections.mean(1),36,axis=0)[:len(faces)];fn=np.cross(v[faces[:,1]]-v[faces[:,0]],v[faces[:,2]]-v[faces[:,0]])
    if np.median(np.sum(radial*fn,axis=1))<0:faces=faces[:,::-1]
    gate=StationPrismConstraints(np.array([sections.reshape(-1,3),inner.reshape(-1,3)]),faces,np.repeat([0,1],18),np.array([True,False]),allow_source_faults=True)
    a0=gate.g0;a1=gate.g1;bernstein=np.stack([a0,a0+a1*.5,a0+a1],axis=0).reshape(-1,3);dot=bernstein@normal;perpendicular=np.linalg.norm(bernstein-dot[:,None]*normal,axis=1);angle=np.arctan2(-dot,perpendicular);half_angle=float(angle.min())
    # Check the omitted collapsed-prism constants: exact construction has
    # no volume at zero station spacing. Bounds retain arithmetic residuals.
    zero0=gate.c0-np.einsum('kfvi,i->kfv',gate.g0,length*normal);zero1=gate.c1-np.einsum('kfvi,i->kfv',gate.g1,length*normal)
    residual=max(float(abs(zero0).max()),float(abs(zero1).max()),float(abs(gate.c2).max()));reserve=-float(dot.max())*length
    certified_angle=max(0.,half_angle-np.deg2rad(.1));cone_reserve=-float((dot*np.cos(certified_angle)+perpendicular*np.sin(certified_angle)).max())*length-residual
    analytic_prism=bool(np.all(gate.sourceHealthy) and certified_angle>0 and cone_reserve>0 and residual<reserve*1e-10)
    root=z['centers'][0];blood=Mesh.load(WORK/'native-parent-envelopes-v9/artery-blood.npz');host=Mesh.load(WORK/'fixed-target-20.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/artery-inner-lineage.json').read_text());outer,_,_,reconstruction=v10.native_outer(host,blood,lineage,root)
    direction=z['centers'][24]-root;direction/=np.linalg.norm(direction);cross=np.cross(normal,direction);s=np.linalg.norm(cross);cosine=normal@direction
    if s<1e-12:rotation=np.eye(3)
    else:
        axis=cross/s;K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]);rotation=np.eye(3)+K*s+K@K*(1-cosine)
    root_inside,ambiguous=ray_contains(root[None],blood.triangles);assert root_inside[0] and not ambiguous[0];root_clearance=BVH(blood).closest_point(root)['distance_m'];blood_exit=first_ray_exit(blood,root,direction);outer_exit=first_ray_exit(outer,root,direction);assert outer_exit['distance_m']>blood_exit['distance_m']
    placements=[];outer_bvh=BVH(outer)
    for extra in np.linspace(.001,.012,12):
        distance=outer_exit['distance_m']+float(extra);posed=local@rotation.T+root+distance*direction;inside,ambiguous=ray_contains(posed,outer.triangles);near,_=outer_bvh.closest_points(posed);minimum=float(near['distance_m'].min());row=dict(distanceFromBloodRoot_m=distance,all18ControlsOutsideNativeOuter=bool(np.all(~inside&~ambiguous)),minimumParentOuterClearance_m=minimum)
        placements.append(row)
        if row['all18ControlsOutsideNativeOuter'] and minimum>=.0001:break
    clear=bool(placements[-1]['all18ControlsOutsideNativeOuter'] and placements[-1]['minimumParentOuterClearance_m']>=.0001)
    result=dict(schema='cystic-port-departure-feasibility/v13',scope='Analytic local departure and parent-clearance feasibility only; no full connector constructed',literal18ControlProfileShapePreserved=True,profileAreaNormal=normal.tolist(),departureCone=dict(halfAngleLowerBoundDegrees=float(np.degrees(certified_angle)),criticalZeroReserveHalfAngleDegrees=float(np.degrees(half_angle)),allAzimuthsAndContinuousWallDepthViaBernsteinBounds=True,nominalOneMillimetreExtrusionPrismPassed=bool(np.all(gate.sourceHealthy)),maximumPrismDeterminant=float(gate.sourceMaximum.max()),collapsedPrismArithmeticResidual=residual,negativeDeterminantReserve=reserve,certifiedConeNegativeDeterminantReserve=cone_reserve,analyticLocalDepartureFeasible=analytic_prism),nativeParent=dict(rootInsideVerifiedBlood=True,bloodSpaceRootClearance_m=root_clearance,smallConnectorRadius_m=.00035,smallConnectorClearanceMargin_m=root_clearance-.00035,bloodExit=blood_exit,outerExit=outer_exit,wallThicknessAlongDeparture_m=outer_exit['distance_m']-blood_exit['distance_m'],exactWallReconstructionResidual_m3=reconstruction),admissibleProposedWorldDirection=direction.tolist(),wholeProfileProperRotation=rotation.tolist(),externalProfilePlacementTrials=placements,externalProfilePlacementFeasible=clear,feasibleLocalPrerequisites=analytic_prism and clear,fullConstructedConnectorAdmissibilityEstablished=False,station7IsNotAnatomicalFixedEndpoint=True,geometryIteration=False,newMeshesWritten=0,organsAndSkeletonUnchanged=True,installed=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'feasibility.json',result);print(json.dumps(dict(output=str(OUT/'feasibility.json'),analyticLocalDepartureFeasible=analytic_prism,coneHalfAngleDegrees=float(np.degrees(half_angle)),externalProfilePlacementFeasible=clear,parentRootClearance_mm=root_clearance*1000,actualParentWallThicknessAlongDeparture_mm=(outer_exit['distance_m']-blood_exit['distance_m'])*1000,newMeshesWritten=0,elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
