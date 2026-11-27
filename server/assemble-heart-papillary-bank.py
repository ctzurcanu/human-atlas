"""Check the complete emitted AV apparatus, including contacts between fans."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.tissue import validate_cells
ap=argparse.ArgumentParser();ap.add_argument('--stages',nargs='+',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);start=time.monotonic();data=bytearray();parts={};rows=[]
for stage in args.stages:
 folder=h/stage;e=json.loads((folder/'emission.json').read_text());assert e['passed']and json.loads((folder/'checked-geometry.json').read_text())['passed']and digest_file(folder/'cells.bin')==e['bank_sha256']
 while len(data)%4:data.append(0)
 offset=len(data);data.extend((folder/'cells.bin').read_bytes());layout=copy.deepcopy(e['parts'])
 for role,p in layout.items():
  for key in ['positions','normals','normalSigns','uvs','indices']:
   if key in p:p[key]+=offset
  p['geometryOrigin']=e['geometry_origin'];p['chunk']=0;assert p['id']not in parts;parts[p['id']]=p
 rows.append(dict(stage=stage,source=e['source'],bank_sha256=e['bank_sha256'],parts=layout,source_meshes=e['source_meshes'],source_catalogue_sha256=e['source_catalogue_sha256'],maximum_control_rounding_m=e['maximum_control_rounding_m']))
atomic_bytes(out/'cells.bin',bytes(data));mini=dict(chunks=[dict(url='/local-models/'+str((out/'cells.bin').relative_to(ROOT/'.local-models')),bytes=len(data))]);reader=AtlasReader();meshes={pid:reader.read_part(mini,p,frame='model:local-ta98:assembled-metres')for pid,p in parts.items()};checks=validate_cells(meshes,max_pairs=1500000,max_tests=16000000,checkpoint_directory=out/'checkpoints');atomic_json(out/'checked-geometry.json',checks);result=dict(passed=checks['passed'],failures=checks['failures'],parts=parts,sources=rows,bank_sha256=digest_file(out/'cells.bin'),bytes=len(data),physical_parts=len(parts),elapsed_seconds=time.monotonic()-start,installed=False);atomic_json(out/'emission.json',result);print(json.dumps({k:v for k,v in result.items()if k not in ['parts','sources']}),flush=True)
