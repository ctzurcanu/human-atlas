"""Rebase only four reviewed vascular carriers and their existing surface views.

No anatomical coordinates, faces, normal words, organ bulk, cameras or public
assets are edited. This emits an uninstalled candidate requiring fresh audits.
"""
from pathlib import Path
import copy, json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.constraints import check_constraint

WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-hepatic-junction-rebase-v1'
CURRENT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
SAVED = ROOT / '.local-models/geometry-workbench/liver-completion/final-connections-v260-male/atlas.json'
EXPECTED = '1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39'
SAVED_SHA = '61fe994ba0182be90ac4c6104d24739b8f4045f6de0721ecd908210945cc50ba'
PHYSICAL = ['ZA:Coeliac trunk', 'TA98HEPATIC:male:artery',
            'ZA:Inferior vena cava (abdominal part)', 'TA98HEPATIC:male:hepatic_vein']
GEOMETRY = ['positions', 'normals', 'normalSigns', 'uvs', 'colors', 'indices',
            'vertexCount', 'indexCount', 'triangles', 'geometryOrigin',
            'geometryScale', 'geometryInstance', 'sourceOffset']


def load(p): return json.loads(Path(p).read_bytes())
def require(v, message):
    if not v: raise ValueError(message)


def main():
    began = time.monotonic()
    inventory_path = ROOT / 'reports/ta98-circulation-inventory.json'
    authority_report = ROOT / 'TA98-CIRCULATION.md'
    inventory = load(inventory_path)
    require(digest_file(CURRENT) == EXPECTED == inventory['source_catalogue_sha256'],
            'Current source or complete authoritative inventory changed')
    require(digest_file(SAVED) == SAVED_SHA, 'Saved source changed')
    readiness = load(WORK / 'source-readiness.json')
    require(readiness['source_sha256'] == EXPECTED and readiness['saved_sha256'] == SAVED_SHA
            and all(p['nonpositive_effective_normal_corners'] == 0 for p in readiness['parts']),
            'Actual source/normal field readiness proof is absent')
    terms = {r['id']: r for r in inventory['terms']}
    roots = ['A12.2.12.012', 'A12.2.12.015', 'A12.3.09.005', 'A12.3.09.001']
    subsets = {}
    for term, expected in zip(roots, [41, 25, 4, 82]):
        selected = {term}
        while True:
            expanded = selected | {r['id'] for r in terms.values() if r['parent_id'] in selected}
            if expanded == selected: break
            selected = expanded
        require(len(selected) == expected, 'Incomplete authoritative subtree: ' + term)
        subsets[term] = [terms[t] for t in sorted(selected)]
    before = load(CURRENT); saved = load(SAVED); after = copy.deepcopy(before)
    old = {p['id']: p for p in before['parts']}; source = {p['id']: p for p in saved['parts']}
    parts = {p['id']: p for p in after['parts']}
    views = sorted(p['id'] for p in before['parts'] if p.get('hepaticRegionOf') in PHYSICAL)
    require(len(views) == 14, 'Expected fourteen existing dependent surface views')
    changed = PHYSICAL + views
    chunks = {}
    for pid in changed:
        record, incoming = parts[pid], source[pid]
        ci = incoming['chunk']
        if ci not in chunks:
            chunks[ci] = len(after['chunks'])
            after['chunks'].append(copy.deepcopy(saved['chunks'][ci]))
        record['chunk'] = chunks[ci]
        for key in GEOMETRY:
            record.pop(key, None)
            if key in incoming: record[key] = copy.deepcopy(incoming[key])
        if pid in PHYSICAL:
            for key in ['pairedTissue', 'sectionSolid', 'tissueDrawGroups']:
                if key in incoming: record[key] = copy.deepcopy(incoming[key])
        if 'hepaticAttributeSignature' in incoming:
            record['hepaticAttributeSignature'] = incoming['hepaticAttributeSignature']
    reader = AtlasReader(); actual = {}; rows = []
    for pid in changed:
        record = parts[pid]
        mesh = reader.read_part(after, record, frame='model:local-ta98:assembled-metres', catalogue_path=CURRENT)
        actual[pid] = mesh
        record['bounds'] = mesh.bounds.tolist()
        if pid in PHYSICAL:
            path = WORK / 'meshes' / (pid.replace(':', '_').replace('/', '_') + '.npz')
            path.parent.mkdir(parents=True, exist_ok=True)
            mesh.save(path)
            rows.append(dict(id=pid, mesh=str(path), sha256=digest_file(path),
                             geometry_sha256=mesh.geometry_sha, attributes_sha256=mesh.attributes_sha,
                             bank=str(reader.asset(after['chunks'][record['chunk']]['url'])),
                             bank_sha256=reader.file_sha(reader.asset(after['chunks'][record['chunk']]['url']))))
    interfaces = []
    for first, second in [(PHYSICAL[0], PHYSICAL[1]), (PHYSICAL[2], PHYSICAL[3])]:
        parameters = {}
        for role, pid, peer in [('first', first, second), ('second', second, first)]:
            faces = [i for g in parts[pid].get('tissueDrawGroups', []) if g.get('neighbor') == peer
                     for i in range(g['start'] // 3, (g['start'] + g['count']) // 3)]
            require(len(faces) == len(set(faces)) == 24, 'Incomplete actual annular seam: ' + pid)
            parameters[role + '_faces'] = faces
        result = check_constraint(dict(kind='shared-interface', parameters=parameters), actual[first], actual[second])
        require(result['passed'], 'Actual emitted seam is not literally reciprocal')
        interfaces.append(dict(kind='shared-interface', first=first, second=second,
                               parameters=parameters, actual_check=result))
    # Compare literal world triangles: view index compaction is independent.
    def triangle_key(t):
        return min(tuple(np.roll(t, k, axis=0).ravel()) for k in range(3))
    host_triangles = {host: {triangle_key(t) for t in actual[host].triangles}
                      for host in {parts[pid]['hepaticRegionOf'] for pid in views}}
    for pid in views:
        host = parts[pid]['hepaticRegionOf']
        require(all(triangle_key(t) in host_triangles[host] for t in actual[pid].triangles),
                'Surface view is not a literal subset of its actual host: ' + pid)
        require(parts[pid].get('pairedTissue') == old[pid].get('pairedTissue')
                and parts[pid].get('sectionSolid') == old[pid].get('sectionSolid'),
                'A surface reference acquired false finite-tissue semantics')
    require(all(parts[pid] == record for pid, record in old.items() if pid not in changed),
            'Unrelated source anatomy changed')
    require(after['concepts'] == before['concepts'] and after['materials'] == before['materials']
            and after['defaultHidden'] == before['defaultHidden']
            and after['chunks'][:len(before['chunks'])] == before['chunks'],
            'Unrelated concepts/materials/visibility or source bank prefix changed')
    marker = dict(schema='vascular-annular-junction-rebase/v1', sex='male',
                  sourceCatalogueSha256=EXPECTED, savedCatalogueSha256=SAVED_SHA,
                  partIds=PHYSICAL, dependentSurfaceViewIds=views, interfaces=interfaces,
                  anatomicalAcceptance=False, installationReady=False,
                  currentContextReviewPending=True, otherTributaryPatencyPending=True,
                  heartGeometryUnchanged=True, liverBulkAndClassifiersUnchanged=True,
                  stomachAndSkeletonUnchanged=True,
                  limitations='Two saved annular vessel junctions only. Current organ containment, other coeliac/IVC tributary openings and heart-facing IVC continuity require fresh review. Arterial origin repacking changes controls by at most5.59nm and retessellates18old/12new intrahepatic facets; no cavity preservation or whole-network patency claim is inferred.')
    after['vascularJunctionIntegration'] = marker
    atomic_json(WORK / 'applicable-terms.json', dict(schema='vascular-junction-authority-subsets/v1',
                source_sha256=EXPECTED, full_inventory=str(inventory_path),
                full_inventory_sha256=digest_file(inventory_path),
                checked_in_report=str(authority_report), checked_in_report_sha256=digest_file(authority_report),
                roots=roots, terms_by_root=subsets))
    atomic_json(WORK / 'candidate-atlas.json', after)
    atomic_json(WORK / 'packet.json', dict(schema='vascular-junction-rebase-packet/v1',
                candidate=str(WORK / 'candidate-atlas.json'), candidate_sha256=digest_file(WORK / 'candidate-atlas.json'),
                source_catalogue=str(CURRENT), source_sha256=EXPECTED,
                saved_catalogue=str(SAVED), saved_sha256=SAVED_SHA,
                physical_parts=rows, dependent_surface_views=views, changed_parts=changed,
                interfaces=interfaces, whole_network_patency_accepted=False,
                context_audit_pending=True, installed=False, elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(candidate=str(WORK / 'candidate-atlas.json'),
                          sha256=digest_file(WORK / 'candidate-atlas.json'), physical_parts=4,
                          dependent_views=14, exact_annular_interfaces=2,
                          context_audit_pending=True, installed=False, elapsed_seconds=time.monotonic()-began)), flush=True)


if __name__ == '__main__': main()
