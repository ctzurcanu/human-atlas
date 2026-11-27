"""Actual ITA mouths, source calibre and sternum/costal recipient witnesses.

This is a landmark feasibility stage, not a selected arterial takeoff.
"""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,projected_area
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.distance import barycentric

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
OUT=BASE/'ita-actual-interface-landmarks-v3'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def segment_hits(start,stop,triangles):
    n=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);ds=np.einsum('ij,ij->i',n,start-triangles[:,0]);dt=np.einsum('ij,ij->i',n,stop-triangles[:,0]);ids=np.flatnonzero(ds*dt<0)
    if not len(ids):return []
    t=-ds[ids]/(dt[ids]-ds[ids]);p=start+t[:,None]*(stop-start);w,valid=barycentric(p,triangles[ids]);keep=valid&np.all(w>=0,axis=1)
    return [dict(sourceFace=int(i),segmentParameter=float(mu),point_m=x,barycentric=ww)for i,mu,x,ww in zip(ids[keep],t[keep],p[keep],w[keep])]

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    a=json.loads(CAT.read_bytes());ps={p['id']:p for p in a['parts']};runtime=json.loads((BASE/'actual-runtime-context-v2.json').read_bytes());available={p['id']for p in runtime['rows']if p['normalAvailable']}
    reader=AtlasReader();meshes={};trees={};bindings={};rows=[]
    def get(pid):
        if pid not in meshes:
            p=ps[pid];m=reader.read_part(a,p,frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED);path=OUT/(identity(pid)[:20]+'.npz');m.save(path);meshes[pid]=m
            bindings[pid]=dict(path=str(path),sha256=digest_file(path),recordSha256=identity(p),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'],normalRuntimeAvailable=pid in available)
        return meshes[pid]
    def tree(pid):
        if pid not in trees:trees[pid]=BVH(get(pid))
        return trees[pid]
    sternum=['ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6']
    for side in ['l','r']:
        pid='ZA:Internal thoracic artery.'+side;m=get(pid);v,first,inv=np.unique(m.vertices,axis=0,return_index=True,return_inverse=True);f=inv[m.faces];course=recover_station_rings(v,f);assert len(course)==1
        rr=list(course[0]['rings'])
        if v[rr[-1]].mean(0)[1]>v[rr[0]].mean(0)[1]:rr=rr[::-1]
        directed=next(q for q in boundary_loops(f)if set(q)==set(rr[0]));rr[0]=directed
        rings,phase=align_station_correspondence(v,f,rr);centers=np.array([v[q].mean(0)for q in rings]);area=np.array([projected_area(v[q])for q in rings]);diam=np.array([np.linalg.norm(v[q][:,None]-v[q][None],axis=2).max()for q in rings])
        stationRows=[];cartilage=['ZA:Costal cartilage of '+n+' rib.'+side for n in ['first','second','third','fourth','fifth','sixth','seventh']]
        for j in range(0,len(rings),4):
            c=centers[j];best=min([(tree(s).closest_point(c),s)for s in sternum],key=lambda q:q[0]['distance_m']);w,sid=best
            endpoint=np.array(w['point']);hits=[]
            for cid in cartilage:
                q=segment_hits(c,endpoint,get(cid).triangles)
                if q:hits.append(dict(partId=cid,hits=q))
            stationRows.append(dict(sourceStation=j,center_m=c,sourceRingVertexIds=first[rings[j]],sourceRingArea_m2=area[j],sourceDiameter_m=diam[j],
                closestActualSternalSourcePart=sid,sternalTriangleWitness=w,actualCostalCartilageSegmentHits=hits,
                sourceSectionPlaneNormal=None,notChosenBranchOrigin=True,periostealRecipientLayerExplicitlyModeled=False))
        mouths=[]
        sub='ZA:Left subclavian artery'if side=='l'else'ZA:Right subclavian artery'
        for label,station,recipients in [('proximal',0,[sub]),('distal',len(rings)-1,['ZA:Musculophrenic artery.'+side,'ZA:Superior epigastric artery.'+side])]:
            q=v[rings[station]];rmesh=Mesh(q,np.c_[np.arange(len(q)),np.roll(np.arange(len(q)),-1),np.arange(len(q))],FRAME);rrtree=BVH(rmesh);recipientRows=[]
            for rid in recipients:
                rm=get(rid);rc=rrtree.closest_mesh(tree(rid));pv,ix=np.unique(rm.vertices,axis=0,return_inverse=True);rf=ix[rm.faces];matches=[]
                for loop in boundary_loops(rf):
                    target=pv[loop]
                    if len(target)!=len(q):continue
                    variants=[]
                    for reverse in [False,True]:
                        for shift in range(len(loop)):
                            z=np.roll(target[::-1]if reverse else target,shift,axis=0);dx=np.linalg.norm(q-z,axis=1);variants.append((float(dx.max()),reverse,shift,float(np.sqrt((dx*dx).mean()))))
                    maximum,reverse,shift,rms=min(variants)
                    matches.append(dict(recipientPhysicalVertexIds=loop,recipientBoundaryPoints_m=target,completeMaximumMismatch_m=maximum,completeRmsMismatch_m=rms,
                        reverse=reverse,shift=shift,centerDistance_m=float(np.linalg.norm(q.mean(0)-target.mean(0))),exactCoordinateMatchedAnnulus=maximum==0))
                detail,gate=pair_contacts(m,rm,0.,first_bvh=tree(pid),second_bvh=tree(rid),max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget']
                recipientRows.append(dict(recipientPartId=rid,closestActualSourceSurface=rc,recipientMouths=matches,wholeSourceContact=gate,sourceFacePairs=detail['pairs'],notAFiniteAnnulus=True))
            mouths.append(dict(kind=label,station=station,sourceRingCompactVertexIds=first[rings[station]],sourceRingPoints_m=q,sourceRingArea_m2=area[station],sourceDiameter_m=diam[station],recipientRows=recipientRows))
        atomic_npz(OUT/(side+'-ita-source-stations.npz'),physical_vertices=v,physical_faces=f,source_to_physical=inv,rings=np.array(rings),centers=centers,areas=area,diameters=diam,source_faces=m.faces)
        # Complete actual source-stem contacts only. This does not pre-clear
        # an authored branch or accept vessel penetration of other tissues.
        contactRows=[]
        for rid in sternum+cartilage+['ZA:Transversus thoracis muscle.'+side,'ZA:Internal intercostal muscles.'+side,'ZA:Innermost intercostal muscles.'+side]:
            detail,gate=pair_contacts(m,get(rid),0.,first_bvh=tree(pid),second_bvh=tree(rid),max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget']
            contactRows.append(dict(partId=rid,gate=gate,sourceFacePairs=detail['pairs']))
        rows.append(dict(side=side,parentPartId=pid,parentTerm='A12.2.08.029',sourceStations=len(rings),sourceProfileAreaRange_m2=[float(area.min()),float(area.max())],sourceDiameterRange_m=[float(diam.min()),float(diam.max())],
            unbranchedLiteralSourceTrunk=True,sourceTakeoffSocketsAvailable=0,actualOriginalMouths=mouths,sternalRecipientStationWitnesses=stationRows,wholeStemContextContacts=contactRows,
            noActualAnatomicalBranchTakeoffDefinedYet=True))
    result=dict(schema='thoracic-parietal-ita-actual-source-interface-landmarks/v3',sourceCatalogueSha256=EXPECTED,sourceBindings=bindings,rows=rows,
        currentRuntimeContextSha256=digest_file(BASE/'actual-runtime-context-v2.json'),noGeometryOrSourceTransform=True,anatomicalAcceptance=False,
        conclusion='Named native ITA trunks are unbranched 109-station outer source cylinders. Native junctions and recipient patches are measured; no missing branch exists inside the trunk and no specific new takeoff follows merely from a nearest bone witness.',
        elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)));assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],
        sides=[dict(side=q['side'],stations=q['sourceStations'],diameterRange_mm=[1000*x for x in q['sourceDiameterRange_m']],sternalWitnesses=len(q['sternalRecipientStationWitnesses']),
            stemContacts=[dict(id=z['partId'],crossings=z['gate']['crossing_face_pairs'])for z in q['wholeStemContextContacts']if z['gate']['crossing_face_pairs']])for q in rows])))

if __name__=='__main__':main()
