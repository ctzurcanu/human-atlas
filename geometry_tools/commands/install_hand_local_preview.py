"""Install the completed imported hand construction in the ordinary local atlas.

This reuses the already reviewed source bank, never an unfinished fitting job.
It does not claim that the remaining wrist/body placement defects are repaired.
Public delivery and the opt-in whole-circulation catalogue remain unchanged.
"""
import argparse
import copy
import json
import time
from pathlib import Path

from geometry_tools.common import ROOT
from geometry_tools.course_io import sha, write_json
from geometry_tools.install import transaction
from geometry_tools.model import AtlasReader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    started = time.monotonic()
    work = ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
    base_path = ROOT/'.local-models/organ-models/spinal-cord/male-body.json'
    source_path = work/'atlas-cardiovascular-review-v37/atlas.json'
    source_bank = work/'atlas-cardiovascular-review-v34/geometry.bin'
    route = ROOT/'app/local-ta98-layer-catalogue.ts'
    output = work/'hand-palms-local-install-v968'
    if output.exists():
        raise ValueError('Install stage already exists; inspect its saved result')
    base, source = [json.loads(p.read_text()) for p in [base_path, source_path]]
    if sha(source_path) != '3af34239777882340c3512cbe8bbacb2d2de00bd23813d6af6ce0a0742986617':
        raise ValueError('Reviewed source catalogue changed')
    if sha(source_bank) != '7a1c5ecc032b044a564b368e8ad440d99f389bbfd6ce8837035b001803f6759a':
        raise ValueError('Reviewed arterial source bank changed')
    reviewed = json.loads((work/'atlas-cardiovascular-review-v34/visual-review.json').read_text())
    if not reviewed['accepted_for_local_preview'] or reviewed['anatomical_acceptance']:
        raise ValueError('Expected the explicitly limited completed source review')
    arterial = json.loads((work/'atlas-cardiovascular-review-v34/emitted-hand-report.json').read_text())
    if not arterial['complete'] or arterial['skeletonChanged'] or len(arterial['pairs']) != 90:
        raise ValueError('Completed source interface construction required')
    selected = {r['partId'] for r in arterial['rows']}
    for concept in source['concepts']:
        if concept['id'].endswith(('-veins_of_left_hand_and_wrist_ID', '-veins_of_right_hand_and_wrist_ID')):
            selected.update(concept['elements'])
    chosen = [p for p in source['parts'] if p['id'] in selected]
    if len(chosen) != 52 or len({p['id'] for p in chosen}) != 52:
        raise ValueError('Expected 34 arterial and 18 venous original hand source meshes')
    if selected.intersection(p['id'] for p in base['parts']):
        raise ValueError('Ordinary atlas already contains imported hand identities')
    after = copy.deepcopy(base)
    chunks = {}
    for part in chosen:
        if part['chunk'] not in chunks:
            chunks[part['chunk']] = len(after['chunks'])
            after['chunks'].append(copy.deepcopy(source['chunks'][part['chunk']]))
        p = copy.deepcopy(part)
        p['chunk'] = chunks[part['chunk']]
        after['parts'].append(p)
        material = p.get('material')
        if material:
            if material in after['materials'] and after['materials'][material] != source['materials'][material]:
                raise ValueError('Material identity collision')
            after['materials'][material] = copy.deepcopy(source['materials'][material])
    old_concepts = {c['id'] for c in base['concepts']}
    new_concepts = []
    for concept in source['concepts']:
        elements = [i for i in concept['elements'] if i in selected]
        if elements:
            if concept['id'] in old_concepts:
                raise ValueError('Source concept conflicts with ordinary atlas')
            new_concepts.append(dict(copy.deepcopy(concept), elements=elements))
    concepts = old_concepts | {c['id'] for c in new_concepts}
    for c in new_concepts:
        if c.get('partOf') not in concepts:
            c.pop('partOf', None)
    after['concepts'].extend(new_concepts)
    bindings = [copy.deepcopy(t) for t in source['ta98TermBindings']
                if t.get('representation') == 'source-identity' and set(t['partIds']).issubset(selected)]
    after.setdefault('ta98TermBindings', []).extend(bindings)
    # Explicit equivalent old hand vessels, already retired in the full review.
    # Unique palmar carpal branches and native tissues are retained.
    stems = ['Common palmar digital arteries', 'Deep palmar arch', 'Deep venous palmar arch',
             'Dorsal carpal anastomosis', 'Dorsal carpal branch of ulnar artery',
             'Dorsal digital arteries of hand', 'Dorsal digital veins of hand', 'Dorsal metacarpal arteries',
             'Palmar digital veins', 'Palmar metacarpal arteries', 'Proper palmar digital arteries',
             'Radial artery', 'Superficial palmar arch', 'Superficial venous palmar arch', 'Ulnar artery',
             'Anterior interosseous artery', 'Posterior interosseous artery']
    stems3 = ['Deep venous palmar ar', 'Dorsal carpal arch', 'Dorsal carpal network',
              'Dorsal digital arteries', 'Dorsal digital vei', 'Palmar metacarpal vei',
              'Palmar venous netwo', 'Princeps pollicis arteries', 'Radialis indicis',
              'Superficial palmar venous ar', 'Ulnar artery (dorsal carpal br)']
    retired = {prefix+stem+'.'+side for prefix, names in [('ZA:', stems), ('O3M:', stems3)]
               for stem in names for side in ['l', 'r']}
    retired &= set(source['defaultHidden']) & {p['id'] for p in base['parts']}
    hidden = set(base.get('defaultHidden', [])) | retired
    after['defaultHidden'] = sorted(hidden)
    after['partsCount'] = len(after['parts'])
    after['triangles'] = sum(p['indexCount']//3 for p in after['parts'])
    after['handLocalInstallation'] = dict(
        ordinaryLocalAtlas=True, importedParts=sorted(selected), sourceCatalogue=str(source_path),
        sourceCatalogueSha256=sha(source_path), sourceArterialBankSha256=sha(source_bank),
        sourceConstructionComplete=True, arterialInterfaces=90, allImportedDetailRetained=True,
        preferredOver=sorted(retired), skeletonChanged=False, anatomicalAcceptance=False,
        remainingLimits=['Wrist capsule, bone and skin clearance remain unfinished.',
                         'Venous joining and forearm-to-native-vessel interfaces are not certified.',
                         'Existing source identity controls do not establish finite vessel walls.'])
    output.mkdir()
    catalogue = output/'atlas.json'
    write_json(catalogue, after)
    # Independently decode every emitted record through its actual new chunk map.
    reader = AtlasReader()
    for source_part, installed_part in zip(chosen, after['parts'][len(base['parts']):]):
        original = reader.read_part(source, source_part, frame='model:local-ta98:assembled-metres', catalogue_path=source_path)
        emitted = reader.read_part(after, installed_part, frame=original.frame, catalogue_path=catalogue)
        if original.geometry_sha != emitted.geometry_sha or original.attributes.keys() != emitted.attributes.keys():
            raise ValueError('Source arrays changed during installation')
        import numpy as np
        if any(not np.array_equal(original.attributes[k], emitted.attributes[k]) for k in original.attributes):
            raise ValueError('Source normal, UV or material attribute changed')
        if time.monotonic()-started > 18:
            raise TimeoutError('18-second native source-replay stage')
    if after['parts'][:len(base['parts'])] != base['parts']:
        raise ValueError('Existing anatomy or skeleton changed')
    old_route = route.read_bytes()
    needle = b"`/local-models/organ-models/spinal-cord/${model==='local-ta98'?'male':'female'}-body.json`:`/local-models/ta98-review/${folder}/${model==='local-ta98'?'male':'female'}/atlas.json`"
    replacement = b"(model==='local-ta98'?'/local-models/geometry-workbench/biodigital-male-circulation-v1/hand-palms-local-install-v968/atlas.json':'/local-models/organ-models/spinal-cord/female-body.json'):`/local-models/ta98-review/${folder}/${model==='local-ta98'?'male':'female'}/atlas.json`"
    if old_route.count(needle) != 1:
        raise ValueError('Unexpected ordinary local routing source')
    result = dict(complete=True, installed=False, importedMeshes=52, arteries=34, veins=18,
                  distinctTA98Terms=len({t['term'] for t in bindings}), sourceIdentityBindings=len(bindings),
                  retiredAlternatives=sorted(retired), baseCatalogue=str(base_path), baseCatalogueSha256=sha(base_path),
                  sourceCatalogueSha256=sha(source_path), catalogueSha256=sha(catalogue),
                  allImportedBuffersReplayed=True, existingAnatomyUnchanged=True, skeletonChanged=False,
                  anatomicalAcceptance=False, sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    if args.apply:
        result['transaction'] = transaction({route: old_route.replace(needle, replacement)},
                                            output/'rollback-route', expected={str(route): sha(route)})
        result['installed'] = True
    write_json(output/'install-report.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
