"""Retain complete near-contact triangle witnesses for the five native fans."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.contact import pair_contacts
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],default='male');ap.add_argument('--stage',required=True);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/args.sex/'atlas.json';a=json.loads(source.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();rows=[];start=time.monotonic()
for i,pid in enumerate(a['finiteCardiacIntegration']['chordalPartIds']):
 if time.monotonic()-start>35:raise TimeoutError('Papillary-contact stage deadline')
 p=parts[pid];pm=parts[p['cardiacChordalPartition']['papillaryPartId']];first=reader.read_part(a,p,frame='model:'+args.sex+':assembled-metres');second=reader.read_part(a,pm,frame=first.frame);witness,summary=pair_contacts(first,second,.0002,max_pairs=500000,max_tests=4000000);assert summary['complete_within_pair_budget'];atomic_npz(out/f'pair-{i}.npz',**witness);minimum=float(witness['distance_m'].min())if len(witness['pairs'])else None;row=dict(chordae=pid,papillary=pm['id'],first_geometry_sha256=first.geometry_sha,second_geometry_sha256=second.geometry_sha,first_attributes_sha256=first.attributes_sha,second_attributes_sha256=second.attributes_sha,minimum_sampled_surface_distance_m=minimum,contact=summary,physical_shared_interface=False,anatomical_attachment_inferred=False);rows.append(row);print(json.dumps(dict(chordae=pid,minimum_surface_distance_m=minimum,crossings=summary['crossing_face_pairs'],coplanar=summary['coplanar_overlap_face_pairs'],touching=summary['touching_face_pairs'])),flush=True)
atomic_json(out/'contacts.json',dict(complete=True,catalogue_sha256=digest_file(source),sex=args.sex,rows=rows,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False))
