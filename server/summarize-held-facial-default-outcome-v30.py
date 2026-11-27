"""Aggregate exact default-scope passes and unknowns; no new anatomy."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file


def main():
    D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'default-runtime-context-v29';assert not(O/'held-default-outcome-v30.json').exists()
    load=lambda p:json.loads(p.read_bytes());binding=load(D/'source-binding-v19/evidence.json');atlas=load(ROOT/binding['snapshot']);parts={p['id']:p for p in atlas['parts']};v17=load(D/'feasibility-v2/evidence.json');v20=load(D/'native-classifiers-v20/evidence.json');v21=load(D/'conservative-envelopes-v21/evidence.json');v25=load(D/'registered-context-v25/evidence.json');v29=load(O/'evidence.json');runtime=load(D/'nearby-unresolved-participation-frame-v25.json');hairpath=D/'scalp-hair-instanced-bounds-v29.json';assert hairpath.exists();hair=load(hairpath)
    old={r['source_object_id']:r for r in v20['rows']};envelopes={r['source_object_id']:r for r in v21['rows']};allunknown=v17['unverified_native_occupied_volumes']+v25['unverified_native_occupied_volumes']+v29['new_unverified_physical_volumes'];branches={}
    for name in ['ascending-palatine','tonsillar','glandular-0','glandular-1']:
        uncertain={};resolved=[];sheets=[]
        for row in allunknown:
            if row['branch']!=name:continue
            part=parts[row['part_id']];key=part.get('sourceId')or part['id'];earlier=old.get(key);env=envelopes.get(key)
            if earlier and earlier['classification']=='closed-source-whole'and any(x['branch']==name and x['empty']for x in earlier['overlaps']):resolved.append(dict(source_object_id=key,method='Exact literal sibling-material closed occupied classifier',proof='native-classifiers-v20/evidence.json'));continue
            if env and any(x['branch']==name and x['certificate_pass']for x in env['results']):resolved.append(dict(source_object_id=key,method='Whole emitted conduit/full source-point strict convex separation certificate',proof='conservative-envelopes-v21/evidence.json'));continue
            if earlier and earlier['classification']in{'actual-open-integumentary-surface','two-triangle-nonvolume-native-reference'}:
                sheets.append(dict(source_object_id=key,part_id=row['part_id'],classification=earlier['classification'],actual_triangle_tests_not_a_complete_volume_or_lumen_proof=True));continue
            entry=uncertain.setdefault(key,dict(source_object_id=key,part_ids=[],reason='Actual rendered/source triangles tested; occupied physical source volume remains unresolved for this branch.',closed_whole_classifier_available_for_additional_pair_check=bool(earlier and earlier['classification']=='closed-source-whole')))
            if row['part_id']not in entry['part_ids']:entry['part_ids'].append(row['part_id'])
        branches[name]=dict(known_actual_crossings=[r for r in v29['known_actual_crossings']if r['branch']==name],known_occupied_overlaps=[r for r in v29['known_occupied_overlaps']if r['branch']==name],remaining_unique_physical_volume_unknowns=list(uncertain.values()),open_skin_or_helper_surface_limits=sheets,previously_resolved_unknowns=resolved,oral_pharyngeal_complete_fluid_exclusion_unresolved=True,geometrically_complete_or_accepted=False)
    files=['source-binding-v19/evidence.json','nearby-derivative-context-inventory-v24.json','registered-context-v25/evidence.json','nearby-unresolved-participation-frame-v25.json','scalp-hair-instanced-bounds-v29.json','default-runtime-context-v29/evidence.json','default-runtime-context-v29/unbuilt-tonsillar-glandular-source-binding.json']
    result=dict(schema='held-facial-current-default-outcome/v30',source_catalogue_sha256=binding['current_catalogue_sha256'],registry_revision=runtime['registry_revision'],candidate_geometry_sha256=v29['candidate_geometry_sha256'],source_evidence=[dict(path=str((D/p).relative_to(ROOT)),sha256=digest_file(D/p))for p in files],per_branch=branches,scalp_hair_actual_instance_bound_evidence=hair,default_alternatives_excluded_only_with_loader_proof=v29['excluded_default_alternative_ids'],oral_airspace_qualification='No complete native oral/pharyngeal closed-fluid exclusion has been established by the saved metadata/geometry evidence; this is not a proof that no such representation could exist elsewhere.',held_v17_and_v25_v28_failures_preserved=True,no_new_route_or_narrower_geometry_constructed=True,no_packing_or_installation=True,anatomical_acceptance=False)
    atomic_json(O/'held-default-outcome-v30.json',result);print(json.dumps(dict(source=binding['current_catalogue_sha256'],per_branch={k:dict(known_crossings=len(v['known_actual_crossings']),physical_unknowns=[x['source_object_id']for x in v['remaining_unique_physical_volume_unknowns']])for k,v in branches.items()},hair_evidence_sha256=digest_file(hairpath))))


if __name__=='__main__':main()
