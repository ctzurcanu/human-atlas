"""Bind full authority, fixed current context, source stations and mouths.

This stage records constraints before any anatomical movement. It emits only
literal source arrays and independently measured interface diagnostics.
"""
from pathlib import Path
import json,sys,time,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
RUNTIME=BASE/'course-contract-v21/runtime-context.json'
OUT=BASE/'course-contract-v23'
FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def polygon_area(v):
    q=v-v.mean(0);return np.linalg.norm(np.cross(q,np.roll(q,-1,axis=0)).sum(0))/2

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False)
    runtime=json.loads(RUNTIME.read_bytes());sha=digest_file(CAT);assert sha==runtime['sourceCatalogueSha256']
    a=json.loads(CAT.read_bytes());parts={p['id']:p for p in a['parts']};reader=AtlasReader();bindings={};meshes={}
    proposal=json.loads((BASE/'boundary-review-v18/proposal.json').read_bytes());summary=json.loads((BASE/'summary-v17.json').read_bytes())
    def get(pid):
        if pid not in meshes:
            m=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
            path=OUT/(identity(pid)[:20]+'.npz');m.save(path);meshes[pid]=m
            bindings[pid]=dict(mesh=str(path),meshSha256=digest_file(path),recordSha256=identity(parts[pid]),
                geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'])
        return meshes[pid]
    rows=[]
    for side in proposal['rows']:
        for p in side['marginProposals']:
            pid=p['sourcePartId'];m=get(pid);recipient=get(p['fullMouthReplay']['recipientPartId']);physical,inverse=np.unique(m.vertices,axis=0,return_inverse=True);f=inverse[m.faces]
            components=recover_station_rings(physical,f);start=np.array(p['fullMouthReplay']['firstSourceRim']['center_m'])
            component=min(components,key=lambda c:min(np.linalg.norm(physical[q].mean(0)-start)for q in [c['rings'][0],c['rings'][-1]]))
            rings=list(component['rings'])
            if np.linalg.norm(physical[rings[-1]].mean(0)-start)<np.linalg.norm(physical[rings[0]].mean(0)-start):rings=rings[::-1]
            directed=next(q for q in boundary_loops(f[component['face_ids']])if set(q)==set(rings[0]));rings[0]=directed
            aligned,phase=align_station_correspondence(physical,f[component['face_ids']],rings);centers=np.array([physical[q].mean(0)for q in aligned])
            assert np.max(abs(centers-np.array(p['sourceCourseRingsFromArch'])))<1e-15
            # Source strip phase is recovered by literal connectivity. The
            # arch ring candidate is selected by the already measured mouth,
            # never by a generic nearest-side-wall query.
            archphysical,archinverse=np.unique(recipient.vertices,axis=0,return_inverse=True)
            archloop=np.array(p['fullMouthReplay']['secondSourceRim']['vertexIds'],int);v0=physical[aligned[0]];y=archphysical[archloop]
            assert np.array_equal(y,np.array(p['fullMouthReplay']['secondSourceRim']['points_m']))
            actualarchloop=next(q for q in boundary_loops(archinverse[recipient.faces])if set(q)==set(archloop))
            archloop=actualarchloop
            choices=[]
            for reverse in [True]:
                for shift in range(len(y)):
                    order=np.roll(archloop[::-1]if reverse else archloop,shift);yy=archphysical[order]
                    x=v0-v0.mean(0);z=yy-yy.mean(0);u,s,vt=np.linalg.svd(x.T@z);R=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
                    replay=x@R.T+yy.mean(0);res=np.linalg.norm(replay-yy,axis=1)
                    choices.append((float(res.max()),float(np.linalg.norm(replay-v0,axis=1).max()),reverse,shift,order,R,yy))
            # All fits within 50 nm share the original Float32 profile up to
            # its actual encoding error. Choose least movement, avoiding an
            # arbitrary twelvefold phase change from last-bit SVD residuals.
            equivalent=[q for q in choices if q[0]<=5e-8];assert equivalent
            equivalent.sort(key=lambda q:(q[1],q[0],q[3]));error,motion,rev,shift,order,R,target=equivalent[0]
            archcompact=np.array([np.flatnonzero(archinverse==i)[0]for i in order])
            area0=polygon_area(v0);area1=polygon_area(target);angle=float(np.arccos(np.clip((np.trace(R)-1)/2,-1,1)))
            npz=OUT/(side['side']+'-'+p['kind']+'-source-stations.npz')
            atomic_npz(npz,source_vertices=m.vertices,source_faces=m.faces,physical_vertices=physical,physical_faces=f,
                source_to_physical=inverse,rings=np.array(aligned),component_source_faces=component['face_ids'],centers=centers,
                first_ring_target=target,arch_order=order,arch_compact_order=archcompact,arch_fit_rotation=R,
                source_profile_area=np.array([polygon_area(physical[q])for q in aligned]))
            rows.append(dict(side=side['side'],kind=p['kind'],term=p['term'],sourcePartId=pid,recipientPartId=p['fullMouthReplay']['recipientPartId'],
                stations=str(npz),stationsSha256=digest_file(npz),sourceRingCount=len(aligned),controlsPerRing=len(aligned[0]),
                footBoundaryStation=p['coarseBoundaryRing'],exactSourceFaceOwnership=p['proposedLiteralFootSourceFaceIds'],
                sourceMouthInterface=dict(sourceRingPhysicalVertexIds=aligned[0],archPhysicalVertexIds=order,archCompactVertexIds=archcompact,
                    properFitMatrix=R,properFitResidualMaximum_m=error,sourceArea_m2=area0,recipientArea_m2=area1,
                    relativeAreaDifference=float(abs(area1/area0-1)),rotationRadians=angle,translation_m=target.mean(0)-v0.mean(0),
                    directedBoundaryWindingOpposed=True,firstFitSourcePointsReplayedLiterally=True,
                    maximumEndpointControlMotion_m=motion,positionCopyRequiresOnlyStorageRoundingAfterProperSectionFit=True,
                    finiteWallOrBloodLumenEstablished=False,originalMismatch_m=p['fullMouthReplay']['maxVertexDifference_m']),
                fixedSourceIntervals=f'All source faces and controls beyond station {p["coarseBoundaryRing"]}; all other source components; native recipient arch fixed.',
                plannedMouthRepair='Copy this complete existing arch boundary, with proper rotation tapered through source rings and no scaling or cap. Preserve source strip topology and surface normals through the identical section rotations.',
                plannedCourseRepair='Lateral only: one measured smooth course displacement into the posterior/lateral native soft-tissue corridor, endpoints and downstream source rings fixed, rigid source-section transport only. No automatic radius reduction.',
                sourceProfileAreaRange_m2=[float(min(polygon_area(physical[q])for q in aligned)),float(max(polygon_area(physical[q])for q in aligned))],
                maximumMotion_m=.010,anatomicalBoundaryProvisional=True))
    for context in runtime['rows']:
        for p in context['rows']:
            if p['retained']:get(p['id'])
    inventory=[]
    for row in summary['fullApplicableLowerLimbInventory']:
        auth=row['authority'];inventory.append(dict(id=auth['id'],parent_id=auth['parent_id'],name_en=auth['name_en'],name_la=auth['name_la'],
            classification='whole-organ',entity='macroscopic named vessel/network or collective vessel container; no new finite wall claimed',
            existingSources=row['directPartIds'],activeExactSources=row['activeExactPartIds'],representationStatus=row['representationStatus'],
            missingRequirements=row['missingRequirements']))
    contract=dict(schema='foot-marginal-course-whole-interface-contract/v23',sourceCatalogueSha256=sha,
        authoritySha256=summary['authoritySha256'],censusSha256=summary['censusSha256'],fullApplicableAuthorityInventory=inventory,
        authorityParentCaveat=proposal['authorityParentCaveat'],primaryReferences=proposal['primaryReferences'][:2],
        fixedCurrentContext=str(RUNTIME),fixedCurrentContextSha256=digest_file(RUNTIME),
        supersedesInvalidPhysicalVertexIndexingContract=str(BASE/'course-contract-v21/contract.json'),
        sourceBindings=bindings,regions=rows,sourceTopologyDefectsRetained=summary['sourceMeshes'],
        sourceSurfaceOnly=True,noBloodLumenOrFiniteWallClaim=True,noGeometryMovementYet=True,
        explicitConstraints=['Fixed complete skeleton, heart, muscles, tendons, nerve courses, arteries, skin and all retained current physical source context.',
            'Fixed recipient dorsal arches. No lumen plugs or duplicate faces; original face, UV and material ownership retained.',
            'Exact complete mouth-coordinate interface after repair, with opposed face winding; global outer-surface self-intersections and all surrounding current contacts checked.',
            'Every station receives one proper rigid transform, same for all source normal vectors; calibrated source polygon area and diameter remain within explicit Float32 storage roundoff.'],
        anatomicalAcceptance=False,installed=False,elapsedSeconds=time.monotonic()-began)
    # Deliberately avoid huge repeated source component arrays in the contract.
    contract['sourceTopologyDefectsRetained']=[dict(partId=x['partId'],topology=x['topology'])for x in summary['sourceMeshes']if 'saphenous'in x['partId']]
    atomic_json(OUT/'contract.json',json.loads(json.dumps(contract,default=clean)))
    report=ROOT/'reports/mesh-audit/REVIEW.md';heading='## Male foot marginal veins: corrected physical mouth indexing contract v23'
    text=report.read_text();assert heading not in text
    lines=['',heading,'',f'Current source: `{sha}`. Authority: `{summary["authoritySha256"]}`. Full current lower-limb artery/vein inventory is retained below before reconstruction. The two marginal terms remain uninstalled. Macroscopic named vessel/network terms are whole-organ or collective containers; their existing source surfaces do not establish finite wall layers, openings or microscopic tissue. No new point, lumen or wall claims are made.','',
        'Both marginal-vein terms retain their canonical TA98 parent under small saphenous vein. Physical medial continuation into the great saphenous vein is a separate anatomical relation. The actual source portions contain 360 medial / 1,728 lateral triangles per side. Source boundary stations 15 / 72 remain provisional malleolar display boundaries. All source profiles, source UVs and original faces must be retained, with proper section transport only.','',
        'The inherited lateral foot source crosses actual calcaneus at 158 triangle pairs per side (stations 14–29 and 35–44). Its original dorsal-arch mouth mismatch is 27.674 µm. The medial arch mouth is visibly oblique, with 1.04444 mm maximum coordinate mismatch. Neither is accepted as a patent tube. The whole source GSV has 221 tibial crossing pairs per side outside the proposed medial foot face interval; these remain explicit source defects.','',
        'Before candidate construction, the actual normal loader retains 105 nearby physical/context records per side across skeletal, muscular, attachments, integumentary, connective, arterial, fascia, nervous and venous systems; 13 alternatives per side are excluded only by executing existing loader rules. Context, recipient arches, skeleton and installed heart remain fixed. Full source-bound records and exact source station ownership are in local `course-contract-v23/contract.json`. Earlier v22 proper mouth-fit rows indexed physical loop IDs as attributed vertex IDs and are rejected; v23 independently replays all saved physical loop coordinates and requires opposed directed boundary winding. No geometry has moved in this contract stage.','',
        '| TA98 | Parent | Name | Current source controls | Status |','|---|---|---|---|---|']
    for r in inventory:
        sources=', '.join(r['existingSources'])or'—';lines.append('| '+ ' | '.join([r['id'],r['parent_id']or'—',r['name_en'].replace('|','/'),sources.replace('|','/'),r['representationStatus']])+' |')
    report.write_text(text+'\n'.join(lines)+'\n')
    assert digest_file(CAT)==sha
    print(json.dumps(dict(contract=str(OUT/'contract.json'),inventoryTerms=len(inventory),fixedSourceMeshes=len(bindings),elapsedSeconds=contract['elapsedSeconds'],
        mouths=[dict(side=r['side'],kind=r['kind'],maxFitResidual_um=1e6*r['sourceMouthInterface']['properFitResidualMaximum_m'],areaRelative=r['sourceMouthInterface']['relativeAreaDifference'],rotationDegrees=180*r['sourceMouthInterface']['rotationRadians']/np.pi)for r in rows])))

if __name__=='__main__':main()
