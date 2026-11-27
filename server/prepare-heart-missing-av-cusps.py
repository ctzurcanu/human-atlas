"""Register the two licensed named anterior AV cusps; audit emitted arrays."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'anterior-valves-v52';out.mkdir(exist_ok=False);path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};registration=parts['TA98CARDIAC:FJ2418']['registration'];scale=np.array(registration['scale']);shift=np.array(registration['translation']);height=a.get('heightAdjustment');shift=shift*height['scale']+np.asarray(height['translationMeters'])-scale*np.asarray(height['translationMeters'])if height else shift;o=np.array([0.,1.375,.03125]);reader=AtlasReader();rows=[]
for pid in ['FJ2421','FJ2420']:
 m=reader.read_part(a,parts[pid],frame='model:local-ta98:assembled-metres',catalogue_path=path);v=((m.vertices*scale+shift-o).astype('float32').astype(float)+o);attrs={k:q for k,q in m.attributes.items()if k!='normals'}
 if 'normals'in m.attributes:
  n=m.attributes['normals'].astype(float)/scale;n/=np.linalg.norm(n,axis=1)[:,None];attrs['normals']=np.rint(n*32767).astype('<i2')
 candidate=Mesh(v,m.faces,m.frame,attrs);candidate.save(out/(pid+'.npz'));proof=audit(candidate,require_closed=True,max_pairs=350000,max_tests=3000000);atomic_json(out/(pid+'-audit.json'),proof);rows.append(dict(source_id=pid,source_geometry_sha256=m.geometry_sha,geometry_sha256=candidate.geometry_sha,passed=proof['passed'],failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs'],triangles=len(candidate.faces),bounds=candidate.bounds.tolist()));print(pid,rows[-1],flush=True)
atomic_json(out/'construction.json',dict(source_catalogue_sha256=digest_file(path),registration=registration,parts=rows,elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False))
