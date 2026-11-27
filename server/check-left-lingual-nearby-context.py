"""One bounded native-frame spatial query for the new lingual tree only."""
from pathlib import Path
import json,sys,time,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
args=argparse.ArgumentParser();args.add_argument('--version',default='v5');args.add_argument('--encoded',action='store_true');args.add_argument('--side',choices=['left','right'],default='left');args=args.parse_args()
T=time.monotonic();D=ROOT/('.local-models/geometry-workbench/circulatory-completion/male-'+args.side+'-lingual-finite-network-'+args.version);O=D/('nearby-context-emitted-v1'if args.encoded else'nearby-context-v1');O.mkdir(exist_ok=False);source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';sha=digest_file(source);a=json.loads(source.read_bytes());new=Mesh.load(D/('independent-check-emitted-v1'if args.encoded else'independent-check-v1')/'actual-new-arterial-wall.npz');used=new.vertices[np.unique(new.faces)];low,high=used.min(0)-.002,used.max(0)+.002;reader=AtlasReader();rows=[]
for part in a['parts']:
    pid=part['id'];text=(pid+' '+part.get('name','')).lower();bounds=part.get('bounds')
    # Only the installed source frame corroborated by the native ECA/tongue
    # contract. Unaligned supplementary source banks cannot certify packing.
    if not pid.startswith('ZA:')or part.get('suppressed')or bounds is None:continue
    if not any(s in text for s in ['nerve','mandible','tooth','teeth','gland','artery','vein','duct']):continue
    if any(s in text for s in ['node','triangle']):continue
    if not(np.all(np.array(bounds)[1]>=low)and np.all(np.array(bounds)[0]<=high)):continue
    if pid in ['ZA:External carotid artery.'+('l'if args.side=='left'else'r'),'ZA:Sublingual gland.'+('l'if args.side=='left'else'r')]:continue
    target=reader.read_part(a,part,frame=new.frame,catalogue_path=source,catalogue_sha=sha);contacts,summary=pair_contacts(new,target,0.,max_pairs=100000,max_tests=1000000);assert summary['complete_within_pair_budget'];name=str(len(rows));atomic_npz(O/(name+'-contacts.npz'),**contacts)
    if summary['crossing_face_pairs']:target.save(O/(name+'-fixed-source.npz'))
    rows.append(dict(part_id=pid,name=part.get('name'),record=part,geometry_sha256=target.geometry_sha,bank_sha256=target.source['bank_sha256'],summary=summary,contact_file=name+'-contacts.npz'))
assert digest_file(source)==sha
failures=[dict(part_id=r['part_id'],crossings=r['summary']['crossing_face_pairs'])for r in rows if r['summary']['crossing_face_pairs']]
report=dict(schema='left-lingual-native-nearby-context/v1',source_catalogue_sha256=sha,new_arterial_geometry_sha256=new.geometry_sha,expanded_bounds_m=[low.tolist(),high.tolist()],native_source_frame_only=True,candidates=len(rows),checked=rows,forbidden_crossings=failures,passed=not failures,recipient_tongue_muscle_and_sublingual_entry_checked_separately=True,elapsed_seconds=time.monotonic()-T)
atomic_json(O/'evidence.json',report);print(json.dumps(dict(candidates=len(rows),forbidden_crossings=failures,passed=not failures,elapsed_seconds=report['elapsed_seconds'])))
