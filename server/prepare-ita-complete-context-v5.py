"""Bind every current thoracic context source and its actual corridor faces.

The clipping is only an exact triangle-AABB broad phase, never an anatomical
exclusion. Full sources stay fixed; every retained triangle keeps its source ID.
"""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file, identity
from geometry_tools.model import AtlasReader, Mesh

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
FRAME='model:local-ta98:assembled-metres'

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--distal-closure',action='store_true'); args=parser.parse_args()
    OUT=BASE/('ita-distal-closure-context-v6'if args.distal_closure else'ita-complete-context-v5')
    began=time.monotonic(); OUT.mkdir(exist_ok=False); assert digest_file(CAT)==EXPECTED
    atlas=json.loads(CAT.read_bytes()); parts={p['id']:p for p in atlas['parts']}
    runtime=json.loads((BASE/'actual-runtime-context-v2.json').read_bytes())
    # Named, finite anatomical corridor: original bilateral whole ITA bounds,
    # plus 12 mm posterior and 3 mm on the other sides. Lateral parasternal
    # correction remains within the original proximal ITA's lateral extent.
    profiles=[np.load(BASE/'ita-actual-interface-landmarks-v3'/(s+'-ita-source-stations.npz'))['physical_vertices'] for s in ['l','r']]
    points=np.vstack(profiles); low=points.min(0)-[.003,.003,.012]; high=points.max(0)+.003
    if args.distal_closure:
        # The native ITA outlet lies between the literal sternal surfaces.
        # Include the proximal native terminal recipients down to the lower
        # limit already covered by the actual 656-record runtime census.
        low[1]=1.285; low[0]=-.065; high[0]=.065; high[2]=.135
    reader=AtlasReader(); rows=[]
    for record in runtime['rows']:
        if not record['normalAvailable']:continue
        pid=record['id']; p=parts[pid]
        m=reader.read_part(atlas,p,frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED)
        tri=m.triangles
        keep=np.flatnonzero(np.all(tri.max(1)>=low,axis=1)&np.all(tri.min(1)<=high,axis=1))
        row=dict(partId=pid,system=record['system'],recordSha256=identity(p),geometrySha256=m.geometry_sha,
            attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'],normalRuntimeAvailable=True,
            wholeSourceTriangles=len(m.faces),sourceFaceIds=keep.tolist(),retainedTriangles=len(keep),wholeSourceBounds_m=m.bounds.tolist())
        if len(keep):
            path=OUT/(identity(pid)[:20]+'.npz'); clipped=Mesh(m.vertices,m.faces[keep],m.frame,m.attributes,dict(m.source,actualSourceFaceIds=keep.tolist(),broadPhaseOnly=True))
            clipped.save(path); row.update(mesh=str(path),sha256=digest_file(path))
        rows.append(row)
        if len(rows)%100==0:atomic_json(OUT/('progress-'+str(len(rows))+'.json'),dict(processed=len(rows),elapsedSeconds=time.monotonic()-began))
    assert len(rows)==656
    result=dict(schema='ita-complete-current-context-corridor/v5',sourceCatalogueSha256=EXPECTED,
        runtimeContextSha256=digest_file(BASE/'actual-runtime-context-v2.json'),conservativeCorridorBounds_m=[low.tolist(),high.tolist()],
        all656RuntimeSourcesDecoded=True,allSystemsIncluded=True,noSourceMovedOrHidden=True,rows=rows,
        candidateMustRemainInsideBoundedCorridor=True,triangleSelectionRule='all original triangles whose individual closed AABB intersects the declared corridor',
        wholeTissueContainmentNotInferredFromClippedFaces=True,elapsedSeconds=time.monotonic()-began)
    if args.distal_closure:
        result.update(schema='ita-adjoining-native-distal-closure-context/v6',
            movingAssembly=['ZA:'+name+'.'+side for side in ['l','r']for name in ['Internal thoracic artery','Musculophrenic artery','Superior epigastric artery']],
            authorizedConstruction='One coupled proper-section C2 candidate, no repeated offset trials',
            sourceCoordinateAxes=dict(superior='+Y',anterior='+Z',patientLeft='+X'),
            originalSternalCrossingStationIntervals=[[66,77]],originalFourthCartilageCrossingStationIntervals=[[79,82]],
            distalControlReleasedBecause='Every native ITA terminal-ring point is bracketed by literal posterior/anterior sternal source faces. This is direct open-source slice enclosure evidence, not a closed-solid assertion.',
            wholeOrganNamedTissueLandmarks=['ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6']+['ZA:Costal cartilage of '+n+' rib.'+s for s in ['l','r']for n in ['first','second','third','fourth','fifth','sixth','seventh']],
            plannedCourse='About 10 mm lateral to measured sternal border, posterior to the complete costal-cartilage neighborhood; one constrained C2 source-ring field, without radius scaling.',
            protectedControls='ITA stations 0 through 12 and both terminal recipients beyond 45 mm physical source arclength from their ITA-adjacent source mouths; disconnected superior-epigastric source component stays literal.',
            knownUnrepairedInterfaces='Source subclavian and terminal casts have no literal matching annuli. The common distal source-mouth relative geometry will be transported coherently, not relabeled a completed lumen.',
            fixedOwners='All other actual anatomy, including bones/cartilages, skin, muscles, nerves, diaphragm, installed heart and original coronaries',
            sourceRaysSha256=digest_file(BASE/'ita-actual-interface-landmarks-v3/independent-parasternal-rays.json'))
    atomic_json(OUT/'evidence.json',result); assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),allSources=len(rows),actualCorridorSources=sum(bool(r['retainedTriangles'])for r in rows),elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__': main()
