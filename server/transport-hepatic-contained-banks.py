"""Transport hepatic walls while retaining fixed interfaces and face ownership.
Split sharp normal corners at identical coordinates, never move tissue to
correct shading. Exact-coordinate topology and full geometry are audited.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.hepatic_course_field import HepaticCourseField,HepaticVeinShear
from geometry_tools.audit import audit
from geometry_tools.algorithms.transport_slivers import regularize_transport_slivers
from geometry_tools.algorithms.adaptive_course_transport import refine_course
ap=argparse.ArgumentParser();ap.add_argument('--version',required=True);ap.add_argument('--role',required=True,choices=['artery','portal','hepatic_vein','biliary']);ap.add_argument('--stage',choices=['construct','audit'],default='construct');ap.add_argument('--deviation-m',type=float,default=1e-5);ap.add_argument('--refine-defects-from');ap.add_argument('--field-version',default='v163');ap.add_argument('--shear-field');args=ap.parse_args()
started=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/('hepatic-contained-banks-'+args.version);out.mkdir(exist_ok=True)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());reader=AtlasReader();parts={p['id']:p for p in a['parts']}
field_path=Path(args.shear_field)if args.shear_field else base/('hepatic-contained-course-'+args.field_version)/'field.npz';field=HepaticVeinShear(json.loads(field_path.read_text()))if args.shear_field else HepaticCourseField(a,field_path,base/'separated-courses-v105/evidence.json')
role=args.role;p=parts['TA98HEPATIC:male:'+role];oldpath=reader.asset(a['chunks'][p['chunk']]['url']);raw=oldpath.read_bytes();count=p['vertexCount'];origin=np.array(p['geometryOrigin'])
if args.stage=='audit':
    record=json.loads((out/(role+'-pending.json')).read_text());layout=record['layout'];raw=Path(record['bank']).read_bytes();n=layout['vertexCount'];local=np.ndarray((n,3),'<f4',buffer=raw,offset=layout['positions']);faces=np.ndarray((layout['indexCount']//3,3),'<u4',buffer=raw,offset=layout['indices']);words=np.ndarray((n,3),'<i2',buffer=raw,offset=layout['normals']);sign=np.ndarray((n,),'i1',buffer=raw,offset=layout['normalSigns']);mesh=Mesh(local.astype(float)+origin,faces,'model:local-ta98:assembled-metres',dict(normals=words,normalSigns=sign))
    assert digest_file(Path(record['bank']))==record['bankSha256']
    result=audit(mesh,require_closed=True,max_pairs=3000000,max_tests=24000000);atomic_json(out/(role+'-audit.json'),result)
    if not result['passed']:raise ValueError(role+': encoded mesh failed '+str(result['failures']))
    record.update(passed=True,geometrySha256=mesh.geometry_sha,auditSeconds=time.monotonic()-started);atomic_json(out/(role+'-evidence.json'),record);print(json.dumps(dict(role=role,triangles=len(faces),passed=True,seconds=record['auditSeconds'])));sys.exit(0)
oldlocal=np.ndarray((count,3),'<f4',buffer=raw,offset=p['positions']);v=oldlocal.astype(float)+origin
oldwords=np.ndarray((count,3),'<i2',buffer=raw,offset=p['normals']);sign=np.ndarray((count,),'i1',buffer=raw,offset=p['normalSigns']).copy()
faces=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices']).copy();source_faces=faces.copy()
source_count=count;source_local=oldlocal.copy();source_words=oldwords.copy()
fine=[]
if args.refine_defects_from:
    previous=base/('hepatic-contained-banks-'+args.refine_defects_from);bad=np.unique(np.asarray(json.loads((previous/(role+'-audit.json')).read_text())['self_intersections']['examples'],dtype=int));parents=np.load(previous/(role+'-parents.npz'))['parent'];fine=np.unique(parents[bad]).tolist()
pretriangulation=None
if role=='hepatic_vein'and fine:
    faces,pretriangulation=regularize_transport_slivers(v,v,faces,field.weight(v)==0,quality_floor=.01,plane_tolerance_m=1e-5,allowed_faces=fine)
v,faces,effective,parent,refinement=refine_course(v,faces,oldwords.astype(float)*sign[:,None],field,deviation_m=args.deviation_m,max_rounds=18,fine_parents=fine)
extra=len(v)-count
oldlocal=np.r_[oldlocal,(v[count:]-origin).astype('<f4')];unit=effective[count:]/np.linalg.norm(effective[count:],axis=1)[:,None]
oldwords=np.r_[oldwords,np.rint(unit*32767).astype('<i2')];sign=np.r_[sign,np.ones(extra,dtype='i1')];count=len(v)
unique,inv=np.unique(v,axis=0,return_inverse=True);mapped=field(unique);J=field.jacobian(unique);det=np.linalg.det(J);fixed=np.all(mapped==unique,axis=1);changed=~fixed[inv]
record=dict(sourceCatalogueSha256=digest_file(path),sourceBankSha256=digest_file(oldpath),field=str(field_path),fieldSha256=digest_file(field_path),minimumSampledJacobian=float(det.min()),maximumMovementMetres=float(np.linalg.norm(mapped-unique,axis=1).max()),fixedCoordinates=int(fixed.sum()),sourceFaceOrderLiteral=True,sourceFaceOwnershipLiteral=True,refinement=refinement,sourcePatchTriangulation=pretriangulation)
atomic_json(out/(role+'-field.json'),record)
if det.min()<=.05:raise ValueError(role+': course field folds or nearly collapses')
local=oldlocal.copy();local[changed]=(mapped[inv[changed]]-origin).astype('<f4');words=oldwords.copy()
normal=np.linalg.solve(J[inv[changed]].transpose(0,2,1),oldwords[changed].astype(float)[:,:,None])[:,:,0];normal/=np.linalg.norm(normal,axis=1)[:,None]
words[changed]=np.rint(np.clip(normal,-1,1)*32767).astype('<i2')
faces,slivers=regularize_transport_slivers(v,local.astype(float)+origin,faces,~changed)
source_faces=faces.copy();record.update(sliverTriangulation=slivers,sourceFaceOrderLiteral=not slivers['changes'])
_,weld=np.unique(local,axis=0,return_inverse=True);wf=weld[faces];keep=(wf[:,0]!=wf[:,1])&(wf[:,1]!=wf[:,2])&(wf[:,2]!=wf[:,0]);record['collapsedAtStorage']=int((~keep).sum());faces=faces[keep];parent=parent[keep];source_faces=faces.copy()
q=local[faces].astype(float);fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);effective=words.astype(float)*sign[:,None];bad=np.argwhere(np.einsum('tki,ti->tk',effective[faces],fn)<=0)
if np.any(np.all(~changed[faces[bad[:,0]]],axis=1)&np.all(faces[bad[:,0]]<source_count,axis=1)):raise ValueError('An original fixed junction face would need a different normal')
source_coordinates=v.copy()
if len(bad):
    duplicate=faces[bad[:,0],bad[:,1]].copy();flat=fn[bad[:,0]];flat/=np.linalg.norm(flat,axis=1)[:,None]
    local=np.r_[local,local[duplicate]];words=np.r_[words,np.rint(flat*sign[duplicate,None]*32767).astype('<i2')];sign=np.r_[sign,sign[duplicate]]
    faces[bad[:,0],bad[:,1]]=np.arange(count,count+len(bad))
    source_coordinates=np.r_[source_coordinates,v[duplicate]]
assert np.array_equal(local[faces],local[source_faces]),'Normal splitting moved geometric triangle corners'
assert np.array_equal(local[:count][~changed],oldlocal[~changed])and np.array_equal(words[:count][~changed],oldwords[~changed])
q=local[faces].astype(float);fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);effective=words.astype(float)*sign[:,None]
if np.any(np.einsum('tki,ti->tk',effective[faces],fn)<=0):raise ValueError('Encoded corner normal is not outward')
data=bytearray();layout=dict(vertexCount=len(local),indexCount=faces.size,triangles=len(faces),geometryOrigin=origin.tolist())
for key,array in [('positions',local.astype('<f4')),('normals',words.astype('<i2')),('indices',faces.astype('<u4')),('normalSigns',sign.astype('i1'))]:
    data.extend(b'\0'*((-len(data))%4));layout[key]=len(data);data.extend(array.tobytes())
target=out/(role+'.bin');atomic_bytes(target,bytes(data));mesh=Mesh(local.astype(float)+origin,faces,'model:local-ta98:assembled-metres',dict(normals=words,normalSigns=sign))
from geometry_tools.common import atomic_npz
atomic_npz(out/(role+'-parents.npz'),parent=parent,source_vertices=source_coordinates)
record.update(bank=str(target),bankSha256=digest_file(target),layout=layout,passed=False,sharpNormalCorners=len(bad),normalSplitMovesCoordinates=False,protectedCoordinatesAndNormalWordsLiteral=True,elapsedSeconds=time.monotonic()-started)
atomic_json(out/(role+'-pending.json'),record);print(json.dumps(dict(role=role,triangles=len(faces),constructed=True,audited=False,seconds=record['elapsedSeconds'])),flush=True)
