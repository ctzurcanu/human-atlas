"""Incrementally extend the saved native domain with all runtime-retained physical sources."""
from pathlib import Path
import json,sys,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';I=D/'global-palatine-graph-v33';O=D/'global-palatine-actual-radius-v39';assert not O.exists();O.mkdir()
    C=I/'construction-contract.json';c=json.loads(C.read_bytes());E=I/'domain-evidence.json';e=json.loads(E.read_bytes());F=I/'additional-participation-frame-v37.json';f=json.loads(F.read_bytes());source=ROOT/c['source_catalogue'];assert digest_file(source)==f['source_catalogue_sha256']==c['source_catalogue_sha256'];assert digest_file(I/'clearance-domain.npz')==e['domain_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    for path,sha in f['code_bindings'].items():assert digest_file(ROOT/path)==sha
    ids=f['retained_indexed_surface_or_tissue_ids'];assert len(ids)==80 and not set(ids)&{r['part_id']for r in e['source_queries']};status={r['id']:r for r in f['rows']};z=np.load(I/'clearance-domain.npz');q=z['points'];clear=z['clearance_lower_bound_m'].copy();limiting=z['limiting_source'].copy();radius=c['source_profile_radius_m'];margin=c['minimum_clearance_m']
    contract=dict(schema='global-palatine-actual-radius-clearance-contract/v39',source_catalogue=c['source_catalogue'],source_catalogue_sha256=c['source_catalogue_sha256'],registry_revision=f['registry_revision'],prior_domain=str((I/'clearance-domain.npz').relative_to(ROOT)),prior_domain_sha256=e['domain_sha256'],runtime_classification=str(F.relative_to(ROOT)),runtime_classification_sha256=digest_file(F),prior_contract_sha256=digest_file(C),
      authority_inventory=c['authority_inventory'],authority_inventory_sha256=c['authority_inventory_sha256'],unchanged_parent_station=8,unchanged_start_m=c['start_m'],unchanged_recipient=c['recipient'],unchanged_SG_SP_band_m=c['SG_SP_band_m'],unchanged_upper_margin_witnesses_m=c['superior_constrictor_upper_witnesses_m'],unchanged_other_route_hashes=c['unmodified_other_routes'],
      reused_literal_sources=99,new_retained_physical_tissue_sources=80,total_declared_sources=179,additional_source_ids=ids,radius_m=radius,margin_m=margin,
      node_rule='Node must retain actual full source radius plus 0.070 mm clearance. No extra grid-cell radius is imposed as anatomy.',
      edge_rule='Use endpoint 1-Lipschitz certificate where valid. Otherwise evaluate the complete finite straight segment against actual source triangles and valid closed-source occupancy, or use adaptive midpoint distance lower bounds. Interface-mask transitions require dedicated analytic source checks.',
      excluded_extra_reference_surfaces=[dict(part_id=r['id'],representation_class=r['representation_class'],basis=r['representation_class_basis'])for r in f['rows']if not r['normal_runtime_suppressed']and r['id']not in ids],
      exclusions_not_provider_based=True,open_physical_volume_unknowns_remain=True,complete_oral_fluid_exclusion_unresolved=True,no_emitted_geometry_pack_install=True)
    atomic_json(O/'construction-contract.json',contract)
    helper_spec=importlib.util.spec_from_file_location('palatine_v39_occupancy',ROOT/'server/check-facial-cervical-feasibility.py');h=importlib.util.module_from_spec(helper_spec);helper_spec.loader.exec_module(h);rows=[];complete=True;checkpoints=[]
    for i,pid in enumerate(ids):
        if time.monotonic()-began>48:complete=False;break
        s=status[pid];r=s['registry_geometry_binding'];m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);assert m.geometry_sha==r['geometry_sha256']and m.attributes_sha==r['attributes_sha256']and m.source['bank_sha256']==r['bank_sha256']
        # Indexed draw bounds ignore unused shared-bank controls but preserve
        # every source triangle. This is a true geometric AABB lower bound.
        used=m.vertices[np.unique(m.faces)];lo,hi=used.min(0),used.max(0);box=np.linalg.norm(np.maximum(lo-q,0)+np.maximum(q-hi,0),axis=1);take=np.flatnonzero(box<.0022);values=np.minimum(box,.010)
        classifier,mode=h.occupancy_classifier(m,parts[pid].get('system')in{'arterial','venous'});inside=np.zeros(len(take),bool)
        if len(take):
            d2,face,near=igl.point_mesh_squared_distance(q[take],m.vertices,m.faces);distance=np.sqrt(d2)
            if classifier is not None:inside=abs(igl.winding_number(classifier.vertices,classifier.faces,q[take]))>.99;distance[inside]*=-1
            values[take]=distance
        changed=values<clear;clear[changed]=values[changed];limiting[changed]=99+i
        rows.append(dict(part_id=pid,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,bank_sha256=m.source['bank_sha256'],frame=frame,representation_class=s['representation_class'],normal_visible_skin0=s['normal_default_visible_skin0'],normal_visible_skin1=s['normal_visible_skin1'],classifier=mode,draw_bounds_m=[lo.tolist(),hi.tolist()],actual_source_controls_unchanged=True,queried_points=len(take),inside_valid_classifier_grid_points=int(inside.sum()),new_limiting_grid_points=int(changed.sum()),radius_margin_violations=int((values<radius+margin).sum())))
        if (i+1)%20==0:
            path=O/('clearance-after-'+str(i+1)+'.npz');atomic_npz(path,clearance_lower_bound_m=clear,limiting_source=limiting);checkpoints.append(dict(completed_sources=i+1,path=str(path.relative_to(ROOT)),sha256=digest_file(path),elapsed_seconds=time.monotonic()-began))
    if len(rows)%20 or not checkpoints:
        path=O/('clearance-after-'+str(len(rows))+'.npz');atomic_npz(path,clearance_lower_bound_m=clear,limiting_source=limiting);checkpoints.append(dict(completed_sources=len(rows),path=str(path.relative_to(ROOT)),sha256=digest_file(path),elapsed_seconds=time.monotonic()-began))
    result=dict(schema='global-palatine-clearance-field-extension/v39',source_catalogue_sha256=c['source_catalogue_sha256'],contract_sha256=digest_file(O/'construction-contract.json'),prior_domain_sha256=e['domain_sha256'],reused_prior_sources=e['source_queries'],new_sources=rows,completed=complete,completed_added_sources=len(rows),total_sources=99+len(rows),grid_points=len(q),actual_radius_plus_margin_m=radius+margin,nodes_at_actual_clearance=int((clear>=radius+margin).sum()),old_extra_half_edge_guard_not_required_for_nodes=True,checkpoints=checkpoints,current_field=checkpoints[-1],open_occupied_volume_unknowns=[r['part_id']for r in rows if r['classifier']in ['source_is_not_a_closed_occupied_volume','invalid_native_conduit_topology','unresolved_native_conduit_boundary','temporary_conduit_classifier_is_not_closed']],no_route_solve_yet=True,no_C2_or_emitted_geometry_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'field-extension-evidence.json',result);assert digest_file(source)==c['source_catalogue_sha256'];print(json.dumps({k:result[k]for k in ['completed','completed_added_sources','total_sources','grid_points','nodes_at_actual_clearance','elapsed_seconds']}))

if __name__=='__main__':main()
