"""Independent v3 cell replay, exact v2 gate reuse, and new cortical containment."""
import importlib.util,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'
OLD=BASE/'new-A-arterial-walls-v2';NEW=BASE/'new-A-arterial-walls-v3';DEST=NEW/'independent-context-v2'

def main():
    began=time.monotonic();assert not DEST.exists();reader=AtlasReader(ROOT);metas={p:json.loads((p/'finite-cells/tissue.json').read_text())for p in [OLD,NEW]};meshes={}
    for directory,meta in metas.items():
        bank=Path(meta['bank']);assert digest_file(bank)==meta['bank_sha256'];atlas={'chunks':[{'url':'/local-models/'+str(bank.relative_to(ROOT/'.local-models'))}]}
        for cell,row in meta['parts'].items():
            m=reader.read_part(atlas,dict(row,chunk=0),frame=meta['frame']);expected=meta['checks']['cells'][cell];assert(m.geometry_sha,m.attributes_sha)==(expected['geometry_sha256'],expected['attributes_sha256']);meshes[(directory,cell)]=m
    ev=json.loads((BASE/'evidence.json').read_text());assert digest_file(ev['catalogue'])==ev['catalogue_sha256'];bindings={r['part_id']:r for r in ev['source_bindings']};reused=[]
    for cell in ['interlobar','arcuate']:
        a,b=meshes[(OLD,cell)],meshes[(NEW,cell)];assert a.frame==b.frame and np.array_equal(a.triangles,b.triangles)
        assert set(a.attributes)==set(b.attributes)
        for attr in a.attributes:assert np.array_equal(a.attributes[attr][a.faces],b.attributes[attr][b.faces])
        prior=OLD/'independent-context-v2'/cell/'evidence.json';x=json.loads(prior.read_text());assert x['passed']and x['all_context_pairs_complete']and x['source_catalogue_sha256']==ev['catalogue_sha256']
        assert(x['candidate_geometry_sha256'],x['candidate_attributes_sha256'])==(b.geometry_sha,b.attributes_sha)
        for context in x['contexts']:
            assert digest_file(context['pairs'])==context['pairs_sha256'];source=context['context_source']
            if context['context']=='actual-urine-facing-lining':
                assert digest_file(source['path'])==source['sha256']and digest_file(source['ownership_archive'])==source['ownership_archive_sha256']
                assert digest_file(bindings['RENALDETAIL:male:left:collecting-wall:mucosa']['mesh'])==source['source_mesh_sha256']
            else:assert digest_file(source['mesh'])==source['mesh_sha256']
        out=DEST/cell;out.mkdir(parents=True);y=dict(x,candidate_bank=metas[NEW]['bank'],candidate_bank_sha256=metas[NEW]['bank_sha256'],candidate_tissue_manifest_sha256=digest_file(NEW/'finite-cells/tissue.json'),
            reuse=dict(previous_evidence=str(prior),previous_evidence_sha256=digest_file(prior),every_actual_world_triangle_literal=True,every_actual_corner_attribute_literal=True,all_context_sources_and_witness_archives_hash_replayed=True),
            current_v3_source_bound_reuse=True)
        atomic_json(out/'evidence.json',y);reused.append(cell)
    # Reuse the independent checker implementation, with the corrected bank as
    # its sole candidate input. It performs all26 actual cortical pair gates.
    spec=importlib.util.spec_from_file_location('renal_independent_context',ROOT/'server/check-new-renal-arterial-context-v1.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);helper.CAND=NEW;helper.main('cortical-radiate')
    cortical=meshes[(NEW,'cortical-radiate')];cortical_tree=BVH(cortical);carrier_path=BASE/'interlobar-routing-grid-v1/registered-native-kidney-outer-carrier.npz';qualification=BASE/'interlobar-routing-grid-v1/capsule-carrier-qualification-v1.json';q=json.loads(qualification.read_text());assert digest_file(carrier_path)==q['outer_mesh']['sha256'];carrier=Mesh.load(carrier_path)
    # Rebind the actual native carrier directly to the current source record.
    current=json.loads(Path(ev['catalogue']).read_text());record=next(p for p in current['parts']if p['id']=='ZA:Kidney.l');actual=reader.read_part(current,record,frame=carrier.frame,catalogue_path=ev['catalogue'],catalogue_sha=ev['catalogue_sha256'])
    assert(actual.geometry_sha,actual.attributes_sha)==(carrier.geometry_sha,carrier.attributes_sha);assert np.array_equal(actual.triangles,carrier.triangles)
    carrier_tree=BVH(carrier);r,s=pair_contacts(cortical,carrier,0.,first_bvh=cortical_tree,second_bvh=carrier_tree,max_pairs=250000,max_tests=2000000);out=DEST/'cortical-radiate';p=out/'carrier-pairs.npz';atomic_npz(p,**r)
    first=cortical.vertices[cortical.faces[0,0]];w=helper.solid_winding(carrier,first);near=carrier_tree.closest_point(first);nearest=carrier_tree.closest_mesh(cortical_tree)
    passed=s['complete_within_pair_budget']and not(s['crossing_face_pairs']or s['coplanar_overlap_face_pairs']or s['touching_face_pairs'])and abs(w)>.5
    containment=dict(schema='independent-new-cortical-carrier-containment/v1',passed=passed,candidate_geometry_sha256=cortical.geometry_sha,candidate_attributes_sha256=cortical.attributes_sha,candidate_bank_sha256=metas[NEW]['bank_sha256'],
        actual_carrier=str(carrier_path),actual_carrier_sha256=digest_file(carrier_path),carrier_geometry_sha256=carrier.geometry_sha,carrier_attributes_sha256=carrier.attributes_sha,
        qualification_source=str(qualification),qualification_sha256=digest_file(qualification),carrier_direct_current_source_record_replayed=True,carrier_source_part='ZA:Kidney.l',carrier_role='Registered hidden native outer reference explicitly retained by finite capsule, not another installed kidney',
        source_catalogue_sha256=ev['catalogue_sha256'],complete_triangle_contact=s,contact_witnesses=str(p),contact_witnesses_sha256=digest_file(p),sample_point_m=first.tolist(),float64_solid_angle_winding=w,sample_point_inside=abs(w)>.5,sample_nearest_surface_m=near['distance_m'],
        exact_triangle_minimum_distance_m=float(nearest['distance_m']),nearest_triangle_witness={k:(v.tolist()if isinstance(v,np.ndarray)else v.item()if isinstance(v,np.generic)else v)for k,v in nearest.items()},
        inference='Both boundaries are connected closed source-checked cells; zero boundary intersections plus an internal source boundary point establishes this vessel wall inside the carrier volume.',
        closed_whole_urine_occupancy_remains_unknown=True,anatomical_acceptance=False,installed=False)
    atomic_json(out/'carrier-containment.json',containment)
    cells=[]
    for cell in ['interlobar','arcuate','cortical-radiate']:
        p=DEST/cell/'evidence.json';e=json.loads(p.read_text());cells.append(dict(cell=cell,evidence=str(p),sha256=digest_file(p),passed=e['passed'],literal_v2_reuse=cell in reused))
    all_pass=all(x['passed']for x in cells)and passed;result=dict(schema='independent-renal-A-wallsv3-context-summary/v1',candidate_bank_sha256=metas[NEW]['bank_sha256'],source_catalogue_sha256=ev['catalogue_sha256'],cells=cells,
        cortical_carrier_containment=str(out/'carrier-containment.json'),cortical_carrier_containment_sha256=digest_file(out/'carrier-containment.json'),passed=all_pass,complete_renal_context_pairs=78,additional_complete_carrier_pairs=1,
        reused_actual_renal_pairs=52,new_actual_renal_pairs=26,closed_whole_urine_occupancy_remains_unknown=True,geometry_or_master_changed=False,installed=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(ev['catalogue'])==ev['catalogue_sha256'];atomic_json(DEST/'summary.json',result);print(json.dumps(dict(output=str(DEST/'summary.json'),passed=all_pass,seconds=result['elapsed_seconds'],reused_cells=reused,carrier_containment=passed,carrier_minimum_distance_m=nearest['distance_m'])))
if __name__=='__main__':main()
