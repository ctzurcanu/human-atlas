"""Prepare one guarded installation of the complete female liver transfer."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
base=ROOT/'.local-models/geometry-workbench/liver-completion';candidate=base/'female-transfer-v206';inputs=base/'female-transfer-v202';stage=base/'female-liver-install-v207';stage.mkdir(exist_ok=False)
a=json.loads((candidate/'atlas.json').read_text());plan=json.loads((inputs/'plan.json').read_text());current=Path(plan['targetCatalogue']);before=json.loads(current.read_text());assert digest_file(current)==plan['targetCatalogueSha256']
for p in [candidate/'independent-transfer-checks.json',candidate/'fixed-context-fit.json',inputs/'clearance-solver.json']:assert json.loads(p.read_text())['passed'],str(p)
old={p['id']:p for p in before['parts']};changed=[p['id']for p in a['parts']if p!=old.get(p['id'])];ids=a['liverTransferIntegration']['partIds'];parts={p['id']:p for p in a['parts']}
closed=[pid for pid in ids if parts[pid].get('pairedTissue')and parts[pid].get('sectionSolid')];assert len(closed)>=21
evidence=[inputs/'plan.json',inputs/'fit-cleared-v204.json',inputs/'clearance-solver.json',candidate/'transfer-evidence.json',candidate/'independent-transfer-checks.json',candidate/'fixed-context-fit.json',base/'completion-v155/coverage.json']
spec=dict(model='local-female-ta98',proposed_catalogue=str(candidate/'atlas.json'),expected_catalogue_sha256=digest_file(current),changed_parts=changed,inventory=str(inputs/'inventory.json'),require_closed=False,closed_parts=closed,max_pairs=3000000,max_tests=24000000,
 assembly=dict(whole_organ='LOCAL:female:liver',part_ids=ids,fixed=['female skeleton','female diaphragm','female skin','all unrelated female anatomy'],
  adjoining_interfaces=[dict(kind='shared-tissue-interface',scope='28 reciprocal tissue interfaces',evidence=str(candidate/'independent-transfer-checks.json')),
   dict(kind='disjoint',scope='fixed diaphragm and 22 nearby skeletal targets',evidence=str(candidate/'fixed-context-fit.json')),
   dict(kind='disjoint',scope='six distinct channel-wall pairs',evidence=str(candidate/'independent-transfer-checks.json')),
   dict(kind='open-conduit',scope='source internal hepatic networks and transferred proximal parents',note='One positive uniform transform preserves the reviewed source assembly.'),
   dict(kind='adjacent-unjoined',scope='native distal female vessels and bile trunk',note='Proximity and physiological relationships do not establish target wall/lumen continuity.')],
  remaining_limits=['Same six unrepresented TA98 terms as the male; 57 organ-scale and six illustrative microscopic representations.','Regional surface, functional and branch boundaries remain provisional.','Target distal vascular/biliary junctions are not certified by this transfer.']),evidence_files=list(map(str,evidence)))
atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec));print(json.dumps(dict(stage=str(stage),changed=len(changed),newParts=len(ids),closed=len(closed))))
