#!/usr/bin/env python3
"""Bounded fixed-context audit of the frozen new pancreatic wall face scope."""
import argparse,json,time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--audit-dir');args=p.parse_args();started=time.monotonic()
 B=ROOT/'.local-models/geometry-workbench/circulatory-completion/pancreatic-public-cohort-v1/finite-network-construction-v1';D=Path(args.audit_dir).resolve()if args.audit_dir else B/'current-context-audit-v1';scopepath=D/'scope-v2-before-tests.json';scope=json.loads(scopepath.read_text());A=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(A)==scope['catalogue_sha256'];atlas=json.loads(A.read_text());parts={p['id']:p for p in atlas['parts']};owners=[q for q in scope['owners']if q['stage']==args.stage];outpath=D/(args.stage+'-checkpoint.json')
 output=json.loads(outpath.read_text())if outpath.exists()else{'schema':'pancreatic-fixed-context-stage/v1','stage':args.stage,'scope_sha256':digest_file(scopepath),'catalogue_sha256':digest_file(A),'owner_ids_fixed_before_tests':[q['id']for q in owners],'results':{},'anatomical_acceptance':False}
 assert output['scope_sha256']==digest_file(scopepath)
 frame='model:local-ta98:assembled-metres';wall=Mesh.load(D/'new-wall-face-scope.npz',frame=frame);tree=BVH(wall);z=np.load(D/'new-wall-face-scope.npz');wholeids=z['whole_wall_face_ids'];occupied=Mesh.load(D/'new-occupied-envelope.npz',frame=frame);origin=occupied.vertices.mean(0);osolid,_=make_solid(occupied.vertices-origin,occupied.faces);reader=AtlasReader();positive=[]
 for q in owners:
  key=q['id']
  if key in output['results']:continue
  if time.monotonic()-started>31:break
  part=parts[key];m=reader.read_part(atlas,part,frame=frame,catalogue_path=A,catalogue_sha=scope['catalogue_sha256']);b=BVH(m);r,s=pair_contacts(wall,m,0.,first_bvh=tree,second_bvh=b,max_pairs=100000,max_tests=1000000);ix=np.flatnonzero(np.isin(r['relation'],['crossing','coplanar-overlap']));witness=[]
  for j in ix[:12]:witness.append({'new_scope_face':int(r['pairs'][j,0]),'whole_wall_face':int(wholeids[r['pairs'][j,0]]),'fixed_face':int(r['pairs'][j,1]),'point_on_new':r['first'][j],'point_on_fixed':r['second'][j],'relation':str(r['relation'][j])})
  row={'source':m.source,'geometry_sha256':m.geometry_sha,'category':q['category'],'triangle_contact_summary':s,'strict_witnesses':witness,'allowed_pancreatic_parenchyma':q['category']=='allowed-pancreatic-parenchyma'}
  if q['category']!='allowed-pancreatic-parenchyma':
   topo=topology(m);row['fixed_topology']=topo
   if topo['closed_oriented_manifold']:
    try:
     solid,_=make_solid(m.vertices-origin,m.faces);inter=osolid^solid;row['new_occupied_envelope_intersection']={'status':str(inter.status()),'volume_m3':float(inter.volume()),'target_signed_volume_m3':float(solid.volume()),'context_closed_oriented':True,'clinical_volume_interpretation_provisional':True,'source_tissueWallLayer':part.get('tissueWallLayer'),'source_pairedTissue':part.get('pairedTissue',False)}
    except Exception as error:row['new_occupied_envelope_intersection']={'unknown':True,'error_type':type(error).__name__,'error':str(error)}
   else:row['new_occupied_envelope_intersection']={'unknown':True,'reason':'Actual fixed source is not a closed oriented manifold; no generic bounding-box occupancy substituted.'}
  output['results'][key]=row;output['elapsed_seconds_current_stage']=time.monotonic()-started;output['complete']=len(output['results'])==len(owners);atomic_json(outpath,output,overwrite=outpath.exists())
  if len(ix):positive.append((key,s['crossing_face_pairs'],s['coplanar_overlap_face_pairs']))
  print(key,'cross',s['crossing_face_pairs'],'coplanar',s['coplanar_overlap_face_pairs'],'complete',s['complete_within_pair_budget'],flush=True)
 output['complete']=len(output['results'])==len(owners);output['remaining_owner_ids']=[q['id']for q in owners if q['id']not in output['results']];output['elapsed_seconds_current_stage']=time.monotonic()-started;atomic_json(outpath,output,overwrite=outpath.exists());print('STAGE',args.stage,'complete',output['complete'],'completed',len(output['results']),'of',len(owners),'positive',positive,flush=True)
if __name__=='__main__':main()
