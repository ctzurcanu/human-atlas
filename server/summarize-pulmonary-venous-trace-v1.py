"""Record source-bound venous review groups without accepting anatomy."""
import json
import math
from pathlib import Path

from geometry_tools.common import digest_file

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    directory=BASE/'venous-trace-evidence-v1'
    evidence=json.loads((directory/'evidence.json').read_text())
    curves=json.loads((BASE/'native-curves-v3/evidence.json').read_text())['curves']
    by_index={c['index']:c for c in curves}
    contract=json.loads((BASE/'family-contract-v2.json').read_text())
    groups=[
        dict(key='right-upper-posterior-intralobar-course',curves=[124,125,127,128,129],
             candidate_terms=['A12.3.02.110'],
             finding='Posterior source c2/c3 follow BII; c5 is the shared central ascending course with c6/c7 tributaries approaching BI. Actual source c2 root matches c5 station48; c6/c7 roots match c5 stations12/24.',
             qualification='Strong literal source drainage-course group; full intralobar ownership requires fissure/lobe surface review.'),
        dict(key='right-upper-posterior-infralobar-course',curves=[126],candidate_terms=['A12.3.02.109'],
             finding='Posterior source c4 enters the exact common root of c1/c5, courses inferiorly toward the BVI/S2 lower edge, and remains distinct from the BII-following c2/c3 tree.',
             qualification='Candidate only; an inferior direction is insufficient to establish infralobar/fissural ownership. c0/c1 are additional anterior/fissural review candidates, not yet assigned.'),
        dict(key='left-superior-segment-tree',curves=[158,160,161,162,163],candidate_terms=['A12.3.02.403'],
             finding='Superior-left source c0 with c2/c3/c4/c5 is a literal connected course tree following BVI (median centre-to-bronchus surfaces 5.4–12.3 mm).',
             qualification='Strong S6-course source group; intrasegmental subdivision still requires tributary territory review.'),
        dict(key='left-anteromedial-basal-tree-misfiled-under-superior-vein',curves=[159,164,165,166],candidate_terms=['A12.3.02.407','A12.3.02.408','A12.3.02.409'],
             finding='Superior-left source c1 with c6/c7/c8 follows actual anteromedial basal/BVIII guides (medians2.3–5.1mm) and shares a literal root centre with the S6 tree, not its distal course. c6/c7 root exactly at c1 station36, c8 at c1 station60.',
             qualification='Strong V8/anteromedial-basal source candidate, not name-only assignment. Exact anterior/medial drainage and intra/intersegmental ownership require rendered bronchovascular review. Do not force these four components into V6 merely because the source object is called Superior vein.'),
        dict(key='right-superior-segment-upper-tree',curves=[168,169,170,171],candidate_terms=['A12.3.02.203'],
             finding='Superior-right c1 plus c2/c3/c4 follow BVI (median7.0–13.6mm), with c2 root at c1 station12 and c3 root at c1 station36; c4 enters c3 station24.',
             qualification='Strong S6-course tree; root stations after72 are shared collector candidates, not automatically intrasegmental.'),
        dict(key='right-superior-segment-other-root-tributaries',curves=[167,172,173],candidate_terms=['A12.3.02.204'],
             finding='Superior-right c0 sweeps laterally across BVI/basal course; c5/c6 enter the common c1 root but descend toward basal territory.',
             qualification='Review candidates only: basal tributaries cannot be called V6 intersegmental solely from distance or common root.'),
        dict(key='left-apicoposterior-own-course-branches',curves=[52,53,54,56,57,59],candidate_terms=['A12.3.02.303'],
             finding='Apicoposterior-left c4/c5/c6 and c8/c9/c11 tributaries follow the actual combined BI+BII guide and root into main c1 at literal stations36 or72.',
             qualification='S1+2-course branch group; intrasegmental identity remains provisional.'),
        dict(key='left-apicoposterior-anterior-facing-branches',curves=[50,51,55,58],candidate_terms=['A12.3.02.304'],
             finding='Apicoposterior-left c2/c3 approach BIII; c7/c10 follow the lateral course nearer BI+BII and BIII than other lobar guides. Their actual parent main c1 roots are at stations60 and48.',
             qualification='Candidate S1+2/S3 drainage interface; no plane classifier or accepted ownership.'),
        dict(key='right-apical-distal-versus-anterior-side-branches',curves=[43,44,45,46,47],candidate_terms=['A12.3.02.103','A12.3.02.104'],
             finding='Apical-right c1 is the long root-to-apex trunk. c2 enters at c1 station48 and carries c0 at its station12; c3/c4 enter c1 at stations72/84.',
             qualification='This supplies an explicit disjoint branch topology for review; V1a/V1b identities need actual S1/S3 drainage validation, not a c2=intersegmental assumption.'),
        dict(key='right-anterior-two-root-sharing-branches',curves=[26,27],candidate_terms=['A12.3.02.106','A12.3.02.107'],
             finding='Anterior-right c0/c1 share literal larger-radius endpoint coordinates within0.02micrometres; both follow BIII more closely than BI/BII.',
             qualification='Two real separate source branches, but neither intra/intersegmental classification is established.'),
        dict(key='left-anterior-lingular-interface-review',curves=[23,24,25],candidate_terms=['A12.3.02.306','A12.3.02.307'],
             finding='Anterior-left c1 roots into c0 station24. c1/c2 follow BIV more closely than BIII; c0 is the mediastinal root-bearing course. c2 has no submicron root match to another source curve.',
             qualification='Potential V3a+b upper/lingular interface, but source drainage and c2 continuity remain unqualified. Do not force V3 intersegmental identity from nearest territory.'),
    ]
    for group in groups:
        group['source_components']=[dict(curve=i,part_id=by_index[i]['part_id'],component=by_index[i]['component'],geometry_sha256=by_index[i]['geometry_sha256'])for i in group['curves']]
        group['accepted_anatomical_assignment']=False
    nodes={r['curve']:r['root_candidate_m']for r in evidence['literal_root_candidates']}
    def gap(i,j):return math.dist(nodes[i],nodes[j])*1000
    collectors=[
        dict(term='A12.3.02.111',name='middle lobe vein',source_curve_roots=[89,104,105],
             finding='Lateral c0 and medial c2 share a literal node; medial c1 roots4.68mm away. No separate literal shared collector stem was recovered.',
             root_gaps_mm={'lateral0_to_medial2':gap(89,105),'lateral0_to_medial1':gap(89,104)},new_geometry_required=True),
        dict(term='A12.3.02.205',name='right common basal vein',source_curve_roots=[66,135,139],
             finding='Inferior basal c2 and the two superior-basal source trees end at distinct roots. No literal common basal stem was recovered.',
             root_gaps_mm={'inferior2_to_superior0':gap(66,135),'inferior2_to_superior4':gap(66,139)},new_geometry_required=True),
        dict(term='A12.3.02.405',name='left common basal vein',source_curve_roots=[60,131],
             finding='Inferior basal c0 and superior basal c0 share a literal root centre, but their differently sized source mouths are disconnected casts, with no common collector stem.',
             root_gaps_mm={'inferior0_to_superior0':gap(60,131)},new_geometry_required=True),
    ]
    missing=[r for r in contract['inventory']if r['id'].startswith('A12.3.02.')and r['representation_status']=='missing-direct-representation']
    inventory=[]
    for term in missing:
        matching=[g['key']for g in groups if term['id']in g['candidate_terms']]
        collector=next((c for c in collectors if c['term']==term['id']),None)
        inventory.append(dict(id=term['id'],name=term['name_en'],parent_id=term['parent_id'],source_review_group_keys=matching,
                              collector_construction=collector,installed=False,accepted=False,
                              status='source-course review group'if matching else'collector construction required'if collector else'unresolved source course',
                              reason=None if matching or collector else'No branch under the proposed source was proven to follow BVIII; construct or trace the correct anterior-basal drainage before partitioning.'))
    result=dict(schema='pulmonary-venous-source-review-groups/v1',source_catalogue_sha256=evidence['source_catalogue_sha256'],
                source_evidence_sha256=digest_file(directory/'evidence.json'),groups=groups,collectors=collectors,missing23_terms=inventory,
                root_cast_interface_gap_summary=[dict(id=r['id'],minimum_source_root_center_to_cast_surface_mm=min(r['larger_end_to_cast_distance_mm']),
                                                     shared_annulus_or_lumen_verified=False)for r in evidence['root_casts']],
                primary_sources=[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC5035324/',use='Right upper vein course and drainage vary; V1b/V2a/V2c are intersegmental, V1a/V2b/V3a intersubsegmental; V2t/V3b/c are surface/fissural.'),
                                 dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC9735388/',use='Actual reconstructed V3a+b and V6b+c course/drainage may lie within upper/S6 territory rather than exactly on a segment plane; do not classify by nearest plane.')],
                geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False)
    path=directory/'source-review-groups.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),sha256=digest_file(path),groups=len(groups),missing_terms=len(inventory),accepted=0)))


if __name__=='__main__':main()
