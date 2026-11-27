"""Actual emitted wall and occupied-wall checks for explicit registry omissions."""
from pathlib import Path
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts


def main():
    began=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'registered-context-v25';assert not O.exists();O.mkdir()
    inventory_path=D/'nearby-derivative-context-inventory-v24.json';inventory=json.loads(inventory_path.read_bytes());expected=inventory['source_catalogue_sha256'];source=ROOT/inventory['source_catalogue'];assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();wall=Mesh.load(D/'finite-wall.npz');data=np.load(D/'surfaces-and-ancestry.npz');ports=json.loads((D/'ports.json').read_bytes());names=ports['component_names'];origin=wall.bounds.mean(0)
    spec=importlib.util.spec_from_file_location('current_context_v25',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());parent,mode=helper.occupancy_classifier(Mesh.load(ROOT/plan['fixed_parent']['mesh']),True);assert parent is not None;parent_solid=helper.solid(parent,origin)
    components={};occupied={}
    for ci,name in enumerate(names[1:],1):
        f=np.concatenate([data['outer_faces'][data['outer_component']==ci],data['inner_faces'][data['inner_component']==ci],data[f'annulus_{ci+1}_faces']]);components[name]=Mesh(data['vertices'],f,wall.frame);z=np.load(D/(name+'-fields.npz'));outer=Mesh(z['fields'][0],z['faces'],wall.frame);classifier,_=helper.occupancy_classifier(outer,True);assert classifier is not None;occupied[name]=helper.solid(classifier,origin)-parent_solid
    rows=[];crossings=[];overlaps=[];unknown=[];cache={};complete=True
    for item in inventory['omitted_candidates']:
        if time.monotonic()-began>48:complete=False;break
        pid=item['part_id'];name=item['branch'];assert item['catalogue_record_matches_registry']and item['frame']==wall.frame and not item['suppressed']and not item['hidden']
        if pid not in cache:
            target=reader.read_part(atlas,parts[pid],frame=wall.frame,catalogue_path=source,catalogue_sha=expected);assert target.geometry_sha==item['geometry_sha256'];classifier,mode=helper.occupancy_classifier(target,parts[pid].get('system')in{'arterial','venous'});solid=helper.solid(classifier,origin)if classifier is not None else None;cache[pid]=(target,classifier,mode,solid)
        target,classifier,mode,solid=cache[pid];contact,summary=pair_contacts(components[name],target,0.,max_pairs=1000000,max_tests=12000000);assert summary['complete_within_pair_budget'],(name,pid)
        intentional=pid=='ZA:Facial artery.l'or(name=='tonsillar'and pid=='ZA:Palatine tonsil.l')
        row=dict(part_id=pid,branch=name,registry_mesh_id=item['registry_mesh_id'],registry_revision=inventory['registry_revision'],geometry_sha256=target.geometry_sha,source_bank_sha256=item['source_bank_sha256'],actual_emitted_outer_inner_and_terminal_annulus_tested=True,source_category=item['source_category'],summary=summary,intentional_source_parent_or_recipient_interface=intentional,occupied_volume_classifier=mode)
        if summary['crossing_face_pairs']:
            file=O/('contacts-'+str(len(rows))+'.npz');atomic_npz(file,**contact);row['contacts']=str(file.relative_to(ROOT))
            if not intentional:crossings.append(dict(branch=name,part_id=pid,crossing_face_pairs=summary['crossing_face_pairs']))
        if solid is not None:
            overlap=occupied[name]^solid;assert overlap.status()==md.Error.NoError;row.update(new_conduit_vs_actual_native_occupied_wall_mm3=float(overlap.volume()),occupied_intersection_empty=bool(overlap.is_empty()))
            if not intentional and not overlap.is_empty():overlaps.append(dict(branch=name,part_id=pid,occupied_overlap_mm3=float(overlap.volume())))
        elif not intentional:unknown.append(dict(branch=name,part_id=pid,reason=mode))
        rows.append(row);atomic_json(O/'progress.json',dict(rows=rows,completed_pairs=len(rows),expected_pairs=len(inventory['omitted_candidates']),elapsed_seconds=time.monotonic()-began),overwrite=True)
    result=dict(schema='explicit-current-registered-facial-context/v25',source_catalogue_sha256=expected,registry_revision=inventory['registry_revision'],inventory_sha256=digest_file(inventory_path),candidate_geometry_sha256=wall.geometry_sha,expected_pairs=len(inventory['omitted_candidates']),completed_pairs=len(rows),complete=complete,actual_crossings=crossings,actual_occupied_overlaps=overlaps,unverified_native_occupied_volumes=unknown,unregistered_or_unresolved_active_sources=inventory['unregistered_or_unresolved_candidates'],unregistered_source_policy='Frame or ownership unresolved; neither silently excluded nor counted as passed native context.',rows=rows,remaining_scope='Closed finite-vessel source walls are checked as their actual occupied wall. A complete foreign-vessel blood/lumen exclusion is not inferred from a tissue-wall-only classifier. Earlier native-source occupancy and full oral-lumen exclusions remain separate unresolved gates.',source_aware_gate_pass=False,installed=False,blender_queued=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);assert digest_file(source)==expected and Mesh.load(D/'finite-wall.npz').geometry_sha==wall.geometry_sha;print(json.dumps({k:v for k,v in result.items()if k not in{'rows','unregistered_or_unresolved_active_sources'}}))


if __name__=='__main__':main()
