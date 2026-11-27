"""Audit native coverings and proximal roots before reconstructing their joins."""
from pathlib import Path
import argparse,json,time
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import audit

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads(path.read_text());reader=AtlasReader();rows=[]
ids=['spinal_cord_dura_mater:959','cauda_equina_base:845','l_cervical_spinal_nerve_c1:892','r_cervical_spinal_nerve_c1:883','l_lumbar_plexus:858','l_sacral_plexus:859']
for tail in ids:
    pid='LOCAL:female:'+tail;p=next((p for p in a['parts']if p['id']==pid),None)
    if p is None:
        if 'cervical_spinal_nerve_c1'in tail:p=next(p for p in a['parts']if p['id'].startswith('LOCAL:female:'+tail.split(':')[0]+':'))
        else:raise ValueError('Missing declared source: '+pid)
    mesh=reader.read_part(a,p,frame='model:local-female-ta98:assembled-metres',catalogue_path=path)
    file=out/(tail.replace(':','-')+'.npz');mesh.save(file)
    report=audit(mesh,require_closed=False,max_pairs=200000,max_tests=1500000)
    atomic_json(file.with_suffix('.json'),dict(part_id=p['id'],report=report))
    row=dict(part_id=p['id'],name=p['name'],mesh=str(file),geometry_sha256=mesh.geometry_sha,
        bounds_m=mesh.bounds.tolist(),failures=report['failures'],topology=report['topology'],
        self_crossings=report['self_intersections']['defect_pairs'])
    rows.append(row)
    print(json.dumps(dict(part=p['name'],triangles=len(mesh.faces),open_edges=report['topology']['open_edges'],
        components=len(report['topology']['components']),self_crossings=row['self_crossings'],elapsed_seconds=time.monotonic()-start)),flush=True)
atomic_json(out/'sources.json',dict(source_catalogue=str(path),source_catalogue_sha256=digest_file(path),sources=rows,
    skeleton_changed=False,geometry_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start))
