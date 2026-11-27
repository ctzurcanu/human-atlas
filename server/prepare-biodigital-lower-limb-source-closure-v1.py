#!/usr/bin/env python3
"""Bounded source hierarchy/face ownership inventory; no spatial queries."""
from pathlib import Path
import hashlib
import json
import time


def main():
    started = time.monotonic()
    base = Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
    source_path = base / 'source-v2/instances.json'
    ownership_path = base / 'body-envelope-fit-v2/candidate-v1/instance-ownership.json'
    manifest_path = base / 'source-v2/manifest.json'
    paths = [source_path, ownership_path, manifest_path]
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    bindings = {str(p): sha(p) for p in paths}
    source = json.loads(source_path.read_text())
    instances, nodes = source['instances'], source['nodes']
    by_id = {x['id']: x for x in instances}
    node_by_id = {x['id']: x for x in nodes}
    candidate = {x['sourceInstanceId']: x for x in json.loads(ownership_path.read_text())['instances']}
    groups, focused_ids = [], set()
    routes = ('popliteal', 'genicular', 'femoral', 'tibia', 'fibula', 'peroneal', 'malleolar', 'calcaneal', 'dorsalis pedis', 'dorsal artery of foot', 'plantar', 'circumflex iliac', 'obturator', 'external iliac')
    major_exact = {f'{side} {route}' for side in ('right', 'left') for route in ('femoral artery', 'femoral vein', 'popliteal artery', 'popliteal vein', 'anterior tibial artery', 'anterior tibial veins', 'posterior tibial artery', 'posterior tibial veins', 'fibular artery', 'fibular veins', 'external iliac artery', 'external iliac vein')}
    for side in ('right', 'left'):
        for family in ('Arteries', 'Veins'):
            group = next(n for n in nodes if n['name'] == f'{family} of {side} lower limb')
            members = [x for x in instances if group['id'] in x['lineage']]
            groups.append({'side': side, 'system': family.lower(), 'sourceGroupId': group['id'], 'parentId': group['parentId'], 'sourceInstanceCount': len(members), 'sourceInstanceIds': [x['id'] for x in members], 'candidateInstanceCount': sum(x['id'] in candidate for x in members)})
    focused = []
    for x in instances:
        name = x['name'].lower()
        if not any(side in name for side in ('right', 'left')) or not any(term in name for term in routes):
            continue
        row = candidate[x['id']]
        focused_ids.add(x['id'])
        regions = []
        if any(term in name for term in ('femoral', 'iliac', 'obturator', 'acetabular')):
            regions.append('hip-thigh')
        if any(term in name for term in ('popliteal', 'genicular', 'tibial recurrent', 'circumflex fibular')) or name in major_exact and any(term in name for term in ('femoral', 'tibial', 'fibular')):
            regions.append('knee-route-and-continuation')
        if any(term in name for term in ('tibia', 'fibula', 'peroneal', 'malleolar', 'calcaneal', 'dorsalis pedis', 'plantar', 'dorsal artery of foot')):
            regions.append('ankle-foot-route-and-continuation')
        if 'nutrient' in name:
            role = 'named nutrient: possible intended bone entry; exact cortical opening and communicating interface unverified'
        elif name in major_exact:
            role = 'major named trunk: prioritize actual solid-bone crossing witnesses; no nutrient exemption'
        elif any(term in name for term in ('perforator', 'perforating')):
            role = 'named perforating route: no automatic permission to cross cortical bone; target tissue/opening unverified'
        elif 'middle genicular' in name or 'acetabular' in name:
            role = 'named joint-directed branch: joint/capsular relation is plausible; a bone-solid crossing is not accepted from the name'
        else:
            role = 'named regional branch: preserve with common parent spatial deformation; route target unverified'
        focused.append({'sourceInstanceId': x['id'], 'name': x['name'], 'parentId': x['parentId'], 'ancestorIds': x['lineage'][:-1], 'sourceGeometryId': x['geometryId'], 'bloodSystem': x['bloodSystem'], 'vascularRole': x['vascularRole'], 'candidateMesh': row['candidateMesh'], 'candidateFaceStart': row['candidateFaceStart'], 'candidateFaceCount': row['candidateFaceCount'], 'candidateVertexStart': row['candidateVertexStart'], 'candidateVertexCount': row['candidateVertexCount'], 'sourceInstanceIndex': row['instanceIndex'], 'regions': regions, 'qualification': role})
    knee_families = []
    for side in ('right', 'left'):
        relevant = [x for x in focused if side in x['name'].lower() and ('knee-route-and-continuation' in x['regions'])]
        knee_families.append({'side': side, 'focusedSourceInstanceIds': [x['sourceInstanceId'] for x in relevant], 'geometryLeafIds': [x['sourceInstanceId'] for x in relevant if not node_by_id[x['sourceInstanceId']].get('children')], 'poplitealBranchGroups': [{'id': n['id'], 'parentId': n['parentId'], 'geometryDescendantIds': [x['id'] for x in instances if n['id'] in x['lineage']]} for n in nodes if n['name'] in (f'Branches of {side} popliteal artery', f'Branches of {side} popliteal vein')]})
    output = base / 'lower-limb-joints-v1'
    output.mkdir(exist_ok=True)
    report = {'schema': 'biodigital-lower-limb-source-closure/v1', 'inputBindings': bindings, 'candidateStage': 'body-envelope-fit-v2/candidate-v1', 'sideSemanticSubtrees': groups, 'focusedRouteCount': len(focused), 'focusedRoutes': focused, 'kneeFamilies': knee_families, 'hierarchyQualification': ['Source nodes are semantic ownership groups; geometry leaves and branch-group siblings do not establish a communicating vascular graph.', 'Popliteal artery and its Branches grouping are siblings under each side lower-limb arterial group. Popliteal vein is in Primary deep veins; its branch grouping is a sibling semantic subtree under Deep veins.', 'No geometric mouths, foramina, exact lumen communication, cortical opening or route clearance are inferred from ancestry or bone centroids.', 'A correction must act as one common spatial field on all 1724 vascular instances so dependent arterial and venous geometry, siblings, original banks and source ancestry remain together. Focused IDs are evidence/query ownership, not a per-vessel displacement mask.'], 'crossingRolePolicy': {'majorTrunks': 'Do not excuse strict cortical-bone crossings as nutrient routes.', 'nutrientNames': 'Keep named nutrient intersections separately qualified until the actual matching bone opening/entry patch is verified; do not mark wrong by default.', 'perforatingNames': 'Perforating/perforator does not establish a bony foramen; do not grant automatic bone-crossing exemption.', 'jointDirectedBranches': 'Middle genicular/acetabular names can target capsular or joint structures; no solid bone crossing is accepted on label alone.'}, 'geometryChanged': False, 'newSpatialQueries': False, 'anatomicalAcceptance': False, 'elapsedSeconds': time.monotonic() - started}
    (output / 'source-closure.json').write_text(json.dumps(report, indent=2) + '\n')
    assert all(sha(p) == value for p, value in ((Path(p), value) for p, value in bindings.items()))
    print(json.dumps({'path': str(output / 'source-closure.json'), 'focusedRouteCount': len(focused), 'kneeCounts': {x['side']: len(x['focusedSourceInstanceIds']) for x in knee_families}, 'semanticSubtreeCounts': {x['sourceGroupId']: x['sourceInstanceCount'] for x in groups}, 'elapsedSeconds': time.monotonic() - started}))


if __name__ == '__main__':
    main()
