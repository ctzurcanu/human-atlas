"""Independent actual-bank containment, including triangle and edge witnesses."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.open3d_backend import library
from geometry_tools.hepatic_course_field import HepaticCourseField
base=ROOT/'.local-models/geometry-workbench/liver-completion';version=sys.argv[1];out=base/('hepatic-contained-banks-'+version);out.mkdir(exist_ok=True);a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader()
body=Mesh.load(base/'vessel-containment-v157/liver-envelope.npz');origin=body.vertices.mean(0);o3d=library();scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((body.vertices-origin).astype(np.float32)),o3d.core.Tensor(body.faces.astype(np.uint32)))
field=HepaticCourseField(a,base/('hepatic-contained-course-'+(sys.argv[2]if len(sys.argv)>2 else'v163'))/'field.npz',base/'separated-courses-v105/evidence.json');sources={};banks={}
for role in ['artery','portal','hepatic_vein','biliary']:
    chosen=dict(artery='v170',portal='v181',biliary='v181',hepatic_vein='v197'if version=='v199'else'v191')[role]if version in ['v193','v199']else('v170'if role=='artery'and version=='v181'else version)
    folder=base/('hepatic-contained-banks-'+chosen);ep=folder/(role+'-evidence.json');ep=ep if ep.exists()else folder/(role+'-pending.json');record=json.loads(ep.read_text());layout=record['layout'];raw=Path(record['bank']).read_bytes();assert digest_file(record['bank'])==record['bankSha256'];p=parts['TA98HEPATIC:male:'+role]
    v=np.ndarray((layout['vertexCount'],3),'<f4',buffer=raw,offset=layout['positions']).astype(float)+layout['geometryOrigin'];faces=np.ndarray((layout['indexCount']//3,3),'<u4',buffer=raw,offset=layout['indices']);old=reader.read_part(a,p,frame=body.frame)
    if(folder/(role+'-parents.npz')).exists():
        z=np.load(folder/(role+'-parents.npz'));parent=z['parent'];source=z['source_vertices']
    else:
        parent=np.arange(len(faces));original_count=p['vertexCount'];lookup={tuple(x):i for i,x in enumerate(v[:original_count])};source=np.array([old.vertices[lookup[tuple(x)]]for x in v])
    labels=np.empty(len(old.faces),object);labels[:]=None
    for g in p['tissueDrawGroups']:labels[g['start']//3:(g['start']+g['count'])//3]=g['replacement']
    sources[role]=dict(v=v,f=faces,source=source,owner=labels[parent]);banks[role]=record['bankSha256']
rows=[];fit_probes=[]
for old in json.loads((base/'vessel-containment-v157/evidence.json').read_text())['measuredIntrahepaticBranches']:
    role=parts[old['id']]['hepaticConnectionRole'];m=sources[role];f=m['f'][m['owner']==old['id']];used=np.unique(f);edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
    p=np.r_[m['v'][used],m['v'][f].mean(1),m['v'][edges].mean(1)];original=np.r_[m['source'][used],m['source'][f].mean(1),m['source'][edges].mean(1)]
    fit_probes.append(original)
    signed=scene.compute_signed_distance(o3d.core.Tensor((p-origin).astype(np.float32)),nthreads=1,nsamples=7).numpy().astype(float);weight=field.weight(original);protected=weight<=.15
    if version=='v199'and role=='hepatic_vein':
        distance=field.cava.compute_distance(o3d.core.Tensor((original-field.query_origin).astype(np.float32)),nthreads=1).numpy()
        protected|=(distance<=.012)&(np.linalg.norm(p-original,axis=1)<=3e-7)
    bad=(signed>3e-7)&~protected
    row=dict(term=old['term'],name=old['name'],samples=len(p),unprotectedOutside=int(bad.sum()),maximumUnprotectedOutsideMetres=float(signed[~protected].max(initial=0)),protectedEntryExitOutside=int(((signed>3e-7)&protected).sum()),originalMaximumOutsideMetres=old['maximum_outside_m'])
    rows.append(row);print(json.dumps(row),flush=True)
result=dict(passed=all(r['unprotectedOutside']==0 for r in rows),banks=banks,sourceCatalogueSha256=digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'),branches=rows,coordinateGuardMetres=3e-7,allUsedVerticesTriangleCentresAndEdgeMidpoints=True,excludedEntryExitProtection='Explicit hilar anchors; native caval vein pedicles within 12 mm of the actual caval wall are exempt only while their source coordinates remain unchanged within the storage guard.',cavalOutletReference='https://pubmed.ncbi.nlm.nih.gov/17957095/')
if not(out/'fit-extra-probes.npz').exists():atomic_npz(out/'fit-extra-probes.npz',points=np.unique(np.concatenate(fit_probes),axis=0))
if not(out/'containment.json').exists():atomic_json(out/'containment.json',result)
else:assert json.loads((out/'containment.json').read_text())==result,'Immutable containment evidence changed'
if not result['passed']:raise ValueError('Actual branch samples still protrude outside the protected entries/exits')
