"""Geometric source-mouth dependencies, distinct from clinical flow acceptance."""
import hashlib
import json
from pathlib import Path
import numpy as np
from .common import ROOT, digest_file
from .model import AtlasReader
from .algorithms.tube_continuity import boundary_loops
from .algorithms.cached_source_mouth import SourceMouthChecks


def check_source_mouth(first, second, first_ids, second_ids, *, tolerance_m=1e-6, area_tolerance=.01):
    """Check complete external boundaries, all correspondences and opposed winding."""
    if first.frame != second.frame:
        raise ValueError('Source mouths must use the same explicit coordinate frame')
    rings, normals = [], []
    for mesh, ids in [(first,first_ids),(second,second_ids)]:
        ids=np.asarray(ids,int)
        if ids.ndim!=1 or len(ids)<3 or np.any(ids<0) or np.any(ids>=len(mesh.vertices)):
            raise ValueError('Every declared mouth control must exist')
        unique,inv=np.unique(mesh.vertices,axis=0,return_inverse=True)
        target=inv[ids]
        actual=[loop for loop in boundary_loops(inv[mesh.faces])if set(loop)==set(target)]
        if len(np.unique(target))!=len(target)or len(actual)!=1:
            raise ValueError('Declared controls are not one complete external source mouth')
        ring=unique[actual[0]];center=ring.mean(0)
        normal=np.cross(ring-center,np.roll(ring,-1,axis=0)-center).sum(0)
        if np.linalg.norm(normal)<=1e-20:
            raise ValueError('Zero-area source mouth')
        normals.append(normal);rings.append(mesh.vertices[ids])
    curve_proof = None
    if rings[0].shape!=rings[1].shape:
        from .algorithms.closed_curve_correspondence import closed_curve_correspondence
        # Use actual source boundary winding, rather than a supplied point
        # permutation. Every corner of each complete cycle participates.
        actual_rings=[]
        for mesh,ids in [(first,first_ids),(second,second_ids)]:
            unique,inv=np.unique(mesh.vertices,axis=0,return_inverse=True)
            target=inv[np.asarray(ids,int)]
            loop=next(x for x in boundary_loops(inv[mesh.faces])if set(x)==set(target))
            actual_rings.append(unique[loop])
        curve_proof=closed_curve_correspondence(*actual_rings,opposite=True)
        gap=curve_proof['maximumGap_m']
    else:
        gap=float(np.linalg.norm(rings[0]-rings[1],axis=1).max())
    area_error=float(abs(np.linalg.norm(normals[1])/np.linalg.norm(normals[0])-1))
    cosine=float(normals[0]@normals[1]/(np.linalg.norm(normals[0])*np.linalg.norm(normals[1])))
    result=dict(passed=bool(gap<=tolerance_m and area_error<=area_tolerance and cosine<-.99),
                maximumGap_m=gap,relativeAreaError=area_error,normalCosine=cosine,
                completeSourceBoundaries=True,geometricCorrespondenceOnly=True,
                wholeLumenPatencyEstablished=False,anatomicalAttachmentEstablished=False)
    if curve_proof is not None:result['curveCorrespondence']=curve_proof
    return result


def source_port_interfaces(atlas, ledger_path=None, *, root=ROOT):
    root=Path(root)
    path=Path(ledger_path)if ledger_path else root/'.local-models/geometry-workbench/biodigital-male-circulation-v1/source-port-dependencies.json'
    if not path.exists():
        return [],dict(available=False,verified=0,pending=0)
    raw=path.read_bytes();ledger=json.loads(raw);ledger_sha=hashlib.sha256(raw).hexdigest();parts={p['id']:p for p in atlas['parts']}
    reader=AtlasReader(root);loaded={};links=[];verified=pending=0;checks=SourceMouthChecks()
    for item in ledger['interfaces']:
        if item['first']not in parts or item['second']not in parts:
            continue
        for pid in [item['first'],item['second']]:
            if pid not in loaded:
                loaded[pid]=reader.read_part(atlas,parts[pid],frame='model:local-ta98:assembled-metres')
        try:
            proof=checks.check(loaded[item['first']],loaded[item['second']],item['firstRing'],item['secondRing'])
        except ValueError as exc:
            proof=dict(passed=False,reason=str(exc),geometricCorrespondenceOnly=True)
        verified+=int(proof['passed']);pending+=int(not proof['passed'])
        links.append(dict(**{'from':item['first'],'to':item['second']},kind='source-mouth-correspondence',
                           coMove=True,verificationPending=not proof['passed'],sourcePortProof=proof,
                           firstRing=item['firstRing'],secondRing=item['secondRing'],
                           sourceLedger=str(path),sourceLedgerSha256=ledger_sha,
                           operation='Keep both complete matched source mouths in the fitting scope. Replay these controls after deformation; geometry-only correspondence does not prove whole flow or clinical correctness.'))
    if digest_file(path)!=ledger_sha:
        raise ValueError('Source-mouth ledger changed during verification')
    return links,dict(available=True,path=str(path),sha256=ledger_sha,verified=verified,pending=pending,
                       anatomicalAttachmentCompleteness=False,wholeLumenPatencyEstablished=False)


def current_ownership_forest(catalogue_path, *, root=ROOT):
    """Select a hash-bound augmented forest without rewriting archived evidence."""
    path=Path(catalogue_path)
    pointer=Path(root)/'.local-models/geometry-workbench/biodigital-male-circulation-v1/current-ownership-forest.json'
    if pointer.exists():
        value=json.loads(pointer.read_text())
        if value['sourceCatalogueSha256']==digest_file(path):
            forest=Path(value['forest'])
            if digest_file(forest)!=value['forestSha256']:
                raise ValueError('Current ownership-forest pointer does not match its immutable file')
            return forest
    return path.parent/'ownership/complete-forest.json'
