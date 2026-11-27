"""Bounded local endpoint feasibility on the same actual gallbladder region."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.thyroid_integer_volume_field import ray_contains
sp=importlib.util.spec_from_file_location('v14',ROOT/'server/construct-admissible-cystic-departure-v14.py');v14=importlib.util.module_from_spec(sp);sp.loader.exec_module(v14)
WORK=v14.WORK;OUT=WORK/'vein-local-entry-v23';CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a'


def rotation(a,b):
    cross=np.cross(a,b);s=np.linalg.norm(cross);c=a@b
    if s<1e-14:return np.eye(3)
    x,y,z=cross/s;K=np.array([[0,-z,y],[z,0,-x],[-y,x,0]]);return np.eye(3)+K*s+K@K*(1-c)


def main():
    began=time.monotonic();assert digest_file(CURRENT)==EXPECTED;assert not OUT.exists();OUT.mkdir();a=json.loads(CURRENT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame=v14.FRAME;artery=reader.read_part(a,parts['TA98CYSTICBLOOD:male:artery'],frame=frame,catalogue_path=CURRENT);artery_bvh=BVH(artery);inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz');inner_bvh=BVH(inner);cover=Mesh.load(WORK/'first-similarity/fixed-covering-free-surface.npz');cover_bvh=BVH(cover)
    z=np.load(WORK/'cystic-vein-v16/sections.npz');center=z['centers'][-1];profile=z['outer_sections'][-1];local=profile-center;near=cover_bvh.closest_point(center);q=np.array(near['point']);t=cover.triangles[near['face_id']];normal=np.cross(t[1]-t[0],t[2]-t[0]);normal/=np.linalg.norm(normal);depth=float((center-q)@normal);assert depth>0
    # Bind the recipient region through its actual free outer faces.
    found=[]
    for p in a['parts']:
        if not p['id'].startswith('TA98PART:ZA:Gallbladder:GALLWALL:covering:')or p.get('suppressed'):continue
        m=reader.read_part(a,p,frame=frame,catalogue_path=CURRENT);ids=[i for g in p['tissueDrawGroups']if not g.get('neighbor')for i in range(g['start']//3,(g['start']+g['count'])//3)];m=Mesh(m.vertices,m.faces[np.array(ids)],frame);hit=BVH(m).closest_point(q);found.append((hit['distance_m'],p['id'],m))
    distance,recipient_id,recipient=min(found,key=lambda x:x[0]);assert distance<1e-10;rb=BVH(recipient)
    tangent=z['centers'][-1]-z['centers'][-2];x=tangent-normal*(tangent@normal);x/=np.linalg.norm(x);y=np.cross(normal,x)
    # Discrete search is explicit: report its resolution rather than claiming
    # an unproved global continuous optimum.
    angles=np.arange(64)*2*np.pi/64;radii=np.arange(1,31)*.0001;offsets=np.array([r*(np.cos(t)*x+np.sin(t)*y)for r in radii for t in angles]);projected,_=rb.closest_points(q+offsets);points=projected['points'];faces=recipient.triangles[projected['face_ids']];normals=np.cross(faces[:,1]-faces[:,0],faces[:,2]-faces[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None];centers=points+depth*normals;shift=np.linalg.norm(centers-center,axis=1);order=np.argsort(shift,kind='stable');accepted=None;tested=0;failures=dict(outside3mm=0,centerTooClose=0,artery=0,inner=0,depth=0,noWallEntry=0)
    # Actual installed artery outer-envelope classifier, quantized in its
    # unchanged bank's coordinates. This disk is construction data only.
    az=np.load(WORK/'admissible-departure-v14/sections.npz');ov,of,_=v14.v10.swept_profiles(az['outer_sections'][7:]);origin=np.array(parts['TA98CYSTICBLOOD:male:artery']['geometryOrigin']);ov=(ov-origin).astype('<f4').astype(float)+origin;outer_artery=Mesh(ov,of,frame)
    disk_faces=np.array([[18,i,(i+1)%18]for i in range(18)]);clearance=.00015
    for i in order:
        if time.monotonic()-began>48:break
        if shift[i]>.003:failures['outside3mm']+=1;continue
        newcenter=centers[i];closest=artery_bvh.closest_point(newcenter)
        if closest['distance_m']<clearance:failures['centerTooClose']+=1;continue
        R=rotation(normal,normals[i]);newprofile=local@R.T+newcenter;disk=Mesh(np.r_[newprofile,newcenter[None]],disk_faces,frame);inside,amb=ray_contains(newcenter[None],outer_artery.triangles)
        if inside[0]or amb[0]:failures['artery']+=1;continue
        _,contact=pair_contacts(disk,artery,clearance,second_bvh=artery_bvh,max_pairs=250000,max_tests=2000000);assert contact['complete_within_pair_budget']
        if contact['returned_triangle_pairs']:failures['artery']+=1;continue
        _,lumen=pair_contacts(disk,inner,0.,second_bvh=inner_bvh,max_pairs=250000,max_tests=2000000);assert lumen['complete_within_pair_budget']
        if lumen['crossing_face_pairs']or lumen['coplanar_overlap_face_pairs']:failures['inner']+=1;continue
        measured=cover_bvh.closest_point(newcenter);actual_depth=measured['distance_m'];depth_error=abs(actual_depth-depth)
        if depth_error>.00002:failures['depth']+=1;continue
        _,entry=pair_contacts(disk,cover,0.,second_bvh=cover_bvh,max_pairs=250000,max_tests=2000000);assert entry['complete_within_pair_budget']
        if not entry['crossing_face_pairs']:failures['noWallEntry']+=1;continue
        accepted=dict(index=int(i),newCenter_m=newcenter.tolist(),shift_m=(newcenter-center).tolist(),shiftLength_m=float(shift[i]),rotation=R.tolist(),targetNormal=normals[i].tolist(),recipientOuterPoint_m=points[i].tolist(),sourceRecipientRegion=recipient_id,sourceRecipientTriangle=int(projected['face_ids'][i]),sourceCenterSignedOuterOffset_m=depth,measuredNewOuterOffset_m=actual_depth,outerOffsetError_m=depth_error,arteryClearanceThreshold_m=clearance,arteryWithinClearancePairs=contact['returned_triangle_pairs'],GBInnermostCrossings=lumen['crossing_face_pairs'],recipientWallEntryPairs=entry['crossing_face_pairs']);atomic_npz(OUT/'selected-endpoint.npz',original_outer_profile=profile,outer_profile=newprofile,center=newcenter,rotation=R,normal=normals[i],old_outer_normal=normal,shift=newcenter-center);disk.save(OUT/'selected-endpoint-construction-disk.npz');break
    proof=dict(schema='cystic-vein-local-entry-feasibility/v23',sourceCatalogueSha256=EXPECTED,stage=2,maximumStages=3,arteryGeometrySha256=artery.geometry_sha,GBRecipientGeometrySha256=recipient.geometry_sha,radialGridResolution_m=.0001,angularGridResolution_degrees=5.625,totalCandidates=len(order),testedCandidates=sum(failures.values())+(accepted is not None),failureCounts=failures,smallestFeasibleTestedShift=accepted,globalContinuousMinimumNotProved=True,maximumShift_m=.003,fullSectionShapeAndCalibrePreservedByProperTransform=True,nativePortalRootAndAllOrgansFixed=True,feasible=accepted is not None,routeConstructed=False,installed=False,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'feasibility.json',proof);print(json.dumps(dict(output=str(OUT/'feasibility.json'),feasible=accepted is not None,selectedShift_mm=None if accepted is None else accepted['shiftLength_m']*1000,failureCounts=failures,elapsedSeconds=proof['elapsedSeconds'])))


if __name__=='__main__':main()
