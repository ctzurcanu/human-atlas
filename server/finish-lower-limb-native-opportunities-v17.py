"""Finish the source-only marginal-vein feasibility summary, without geometry writes."""
from pathlib import Path
import sys, time, json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
import importlib.util
helper=importlib.util.spec_from_file_location('lower_limb_v16',Path(__file__).with_name('audit-lower-limb-native-opportunities-v16.py'))
module=importlib.util.module_from_spec(helper); helper.loader.exec_module(module)
physical,course,clean,CAT,OUT,FRAME=(getattr(module,k) for k in ['physical','course','clean','CAT','OUT','FRAME'])

def rim_pair(first,second):
    a=np.asarray(first['points_m']); b=np.asarray(second['points_m']); variants=[]
    for reverse in [False,True]:
        v=b[::-1] if reverse else b
        for shift in range(len(v)):
            d=np.linalg.norm(a-np.roll(v,shift,axis=0),axis=1)
            variants.append((float(d.max()),reverse,shift,float(np.sqrt((d*d).mean()))))
    maximum,reverse,shift,rms=min(variants)
    def area(q): return .5*np.linalg.norm(np.cross(q,np.roll(q,-1,axis=0)).sum(0))
    return dict(correspondence='complete cyclic 12-vertex original source boundary, no fitting',
                reverse=reverse,cyclicShift=shift,maxVertexDifference_m=maximum,rmsVertexDifference_m=rms,
                firstArea_m2=area(a),secondArea_m2=area(b),firstSourceRim=first,secondSourceRim=second,
                literalSharedCoordinateAnnulus=bool(maximum==0.),pairedInnerOuterWall=False,
                patencyEstablished=False)

def main():
    began=time.monotonic(); before=digest_file(OUT); r=json.loads(OUT.read_text()); assert digest_file(CAT)==r['sourceCatalogueSha256']
    a=json.loads(CAT.read_text()); parts={p['id']:p for p in a['parts']}; sources={p['partId']:p for p in r['sourceMeshes']}; reader=AtlasReader()
    rims=[]; extra=[]
    for side in ['l','r']:
        great=sources['ZA:Great saphenous vein.'+side]; small=sources['ZA:Small saphenous vein.'+side]; arch=sources['ZA:Dorsal venous arch of foot.'+side]
        for donor,term in [(great,'A12.3.11.021'),(small,'A12.3.11.020')]:
            best=min(((np.linalg.norm(np.array(q['center_m'])-p['center_m']),q,p) for q in donor['exactCoordinateBoundaryLoops']for p in arch['exactCoordinateBoundaryLoops']),key=lambda x:x[0])
            rr=rim_pair(best[1],best[2]); rr.update(term=term,side=side,sourcePartId=donor['partId'],recipientPartId=arch['partId']); rims.append(rr)
        pid=great['partId']; m=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=r['sourceCatalogueSha256'])
        sourceFaces=np.array(great['components'][0]['sourceFaceIds']); cm=Mesh(m.vertices,m.faces[sourceFaces],FRAME)
        pm=physical(cm); loops=boundary_loops(pm.faces); cr=course(pm,loops)
        great['components'][0].update(exactCoordinateTopology=topology(cm),sourceCourse=cr,
            sourceOwnership='Only literal component 0 is the principal ankle-to-groin course. Component 1 is a separate leg course and must not be silently included in a marginal-vein control.')
        bid='ZA:Tibia.'+side+':Bone'; b=reader.read_part(a,parts[bid],frame=FRAME,catalogue_path=CAT,catalogue_sha=r['sourceCatalogueSha256'])
        d,g=pair_contacts(m,b,0.,max_pairs=250000,max_tests=2000000); assert g['complete_within_pair_budget']
        extra.append(dict(first=pid,second=bid,actualContactGate=g,sourceFacePairs=d['pairs'],
            firstSourceTriangles_m=m.triangles[np.unique(d['pairs'][:,0])] if len(d['pairs']) else [],
            note='Whole current native source screen. Crossings are inherited and are not resolved by adding a reference.'))
    r.update(schema='lower-limb-native-opportunities/v17',priorSummarySha256=before,
        sourceMouthReplay=rims,additionalNativePackingChecks=extra,
        opportunityFamily=dict(terms=[dict(id='A12.3.11.020',name='lateral marginal vein',
            physicalSourceBySide=['ZA:Small saphenous vein.l','ZA:Small saphenous vein.r']),
            dict(id='A12.3.11.021',name='medial marginal vein',
            physicalSourceBySide=['ZA:Great saphenous vein.l#component-0','ZA:Great saphenous vein.r#component-0'])],
            representationIfReviewed='Non-owning bounded anatomical course controls referencing existing source faces; parent whole vessel retains its physical ownership. A partition would instead need complete disjoint ancestry and a whole-assembly control.',
            sourceAnatomyEvidence='Actual source courses terminate on the lateral and medial original dorsal arch rims. This supports a foot continuation; it does not independently identify the proximal named transition at each malleolus.',
            exactBoundaryStatus='HELD: no independent current medial/lateral malleolus reference or source-face region establishes the ankle transition. No horizontal plane, nearest whole-bone triangle or arbitrary axial fraction is substituted.',
            proposedNextStep='One source-only local ankle/foot anatomical review, identifying the actual medial/lateral malleolar transition on literal source courses. Keep all organs and bones fixed. Only after that define bounded source-face intervals and reciprocal parent controls.',
            flowInterfaceChecksRequired=['Full original dorsal-arch mouth correspondence and normals, including mismatched original coordinates.',
                'Confirm single source wall is only an outer vascular representation; no double-wall or lumen claim.',
                'Retain literal existing proximal saphenous course, side and all other physical records.',
                'Review actual distal source placement against ankle bones; inherited great-saphenous/tibia contacts need anatomical screening.']),
        rejectedAlternatives=[dict(sourcePartIds=['O3M:Dorsal venous netwo.l','O3M:Dorsal venous netwo.r'],
            reason='Actual bounds place this truncated-name network at hand height, not on the foot; its nearest actual foot-arch distance is 824.568 mm. No foot identity is assigned.'),
            dict(sourcePartIds=['FJ2703','FJ2704','FJ2705','FJ2727','FJ2728','FJ2729'],
            reason='Original BP3D left/right marginal-vein FMA4708/FMA4716 belong to cardiac veins (A12.3.01), not medial/lateral foot marginal veins. Excluded by original identity, not a loose marginal-vein name match.')],
        readyNewTerms=[],identificationStatus='Two anatomically plausible literal-course opportunities, zero ready controls until ankle boundaries and source packing are reviewed.',
        finalStageElapsedSeconds=time.monotonic()-began,noGeometryWrites=True,noLiveChanges=True)
    output=OUT.with_name('summary-v17.json')
    assert digest_file(CAT)==r['sourceCatalogueSha256']; atomic_json(output,json.loads(json.dumps(r,default=clean)))
    print(json.dumps(dict(output=str(output),readyNewTerms=[],opportunityTerms=['A12.3.11.020','A12.3.11.021'],
        sourceRimMatches=[dict(term=q['term'],side=q['side'],maxVertexDifference_nm=1e9*q['maxVertexDifference_m'],rms_nm=1e9*q['rmsVertexDifference_m'])for q in rims],
        inheritedGreatSaphenousTibiaChecks=[dict(side=q['first'][-1],crossingPairs=q['actualContactGate']['crossing_face_pairs'])for q in extra],
        elapsedSeconds=r['finalStageElapsedSeconds'])))

if __name__=='__main__': main()
