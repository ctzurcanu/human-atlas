"""Selective licensed source retrieval; no registration or anatomical edits."""
from pathlib import Path
import sys,time,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts

CAT=ROOT/'public/models/atlas.json'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1/bp3d-anterior-circle-source-v6'
FRAME='source:BodyParts3D-4.0:packed-native-metres-y-up'

def clean(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)

def main():
    started=time.monotonic();OUT.mkdir(exist_ok=False);sha=digest_file(CAT);atlas=json.loads(CAT.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    wanted=['FJ1654','FJ1654M','FJ1655','FJ1682','FJ1682M'];meshes={};rows=[]
    for pid in wanted:
        p=parts[pid];m=reader.read_part(atlas,p,frame=FRAME,catalogue_path=CAT,catalogue_sha=sha);meshes[pid]=m
        rows.append(dict(id=pid,originalName=p['name'],originalFmaId=p['conceptId'],record=p,frame=FRAME,sourceBankSha256=m.source['bank_sha256'],geometrySha256=m.geometry_sha,attributeSha256=m.attributes_sha,dimensions_m=np.ptp(m.vertices,axis=0).tolist(),actualBounds_m=m.bounds.tolist(),topology=topology(m),catalogueSourceTriangles=atlas.get('sourceTriangles'),optimizedCatalogue=atlas.get('optimized'),licence='CC BY 4.0; official BodyParts3D source attribution in ATTRIBUTION.md',notRegisteredToCurrentTa98Body=True))
        m.save(OUT/(pid+'.npz'))
    pairs=[]
    for aid,bid,role in [('FJ1654','FJ1655','right ACA–ACom'),('FJ1654M','FJ1655','left ACA–ACom'),('FJ1654','FJ1682','right ACA–ICA'),('FJ1654M','FJ1682M','left ACA–ICA')]:
        a,b=meshes[aid],meshes[bid];nearest=BVH(a).closest_mesh(BVH(b));af,bf=nearest['first_face_id'],nearest['second_face_id']
        aa={tuple(q)for q in a.vertices};bb={tuple(q)for q in b.vertices};shared=sorted(aa&bb)
        contacts,summary=pair_contacts(a,b,0.,max_pairs=250000,max_tests=2000000);assert summary['complete_within_pair_budget']
        pairs.append(dict(role=role,firstPartId=aid,secondPartId=bid,closest=nearest,firstTriangle_m=a.triangles[af],secondTriangle_m=b.triangles[bf],literalSharedCoordinates_m=shared,literalSharedCoordinatesCount=len(shared),actualZeroDistanceContactSummary=summary,actualZeroDistanceSourceFacePairs=contacts['pairs'],anatomicalInterfaceInferred=False))
    result=dict(schema='licensed-bp3d-anterior-circle-source-correspondence/v6',sourceCatalogue=str(CAT),sourceCatalogueSha256=sha,sourceFrame=FRAME,sourceProvenance=atlas['provenance'],frameDefinition='Existing packed native BodyParts3D source frame: original OBJ millimetres/Z-up converted by scripts/convert-anatomy.py to [x*.001, z*.001+.0781112, -y*.001-.1], then existing source catalogue optimization. No transform to the current TA98 anatomy is inferred or applied.',frameCodeSha256=digest_file(ROOT/'scripts/convert-anatomy.py'),sourceRecords=rows,pairs=pairs,onlyFiveSpecifiedRecordsDecoded=True,wholeGlbDecoded=False,noGeometryRegistrationOrLiveChanges=True,noA1A2SubdivisionClaimed=True,anteriorCircleCoherenceEstablished=False,elapsedSeconds=time.monotonic()-started)
    atomic_json(OUT/'proof.json',json.loads(json.dumps(result,default=clean)));assert digest_file(CAT)==sha
    print(json.dumps(dict(output=str(OUT/'proof.json'),sourceSha256=sha,measurements=[dict(pair=p['role'],gap_mm=1000*p['closest']['distance_m'],sharedCoordinates=p['literalSharedCoordinatesCount'],crossingPairs=p['actualZeroDistanceContactSummary']['crossing_face_pairs'],touchingPairs=p['actualZeroDistanceContactSummary']['touching_face_pairs'])for p in pairs],elapsedSeconds=result['elapsedSeconds']),default=clean))

if __name__=='__main__':main()
