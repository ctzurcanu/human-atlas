"""Actual finite-wall occupied-volume check for two selected channel banks."""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'contained-hepatic-pairs-v198';out.mkdir(exist_ok=True);began=time.monotonic()
first,second=sys.argv[1:3];sources=dict(artery='v170',portal='v181',biliary='v181',hepatic_vein='v197');records={}
def load(role):
    folder=base/('hepatic-contained-banks-'+sources[role]);p=folder/(role+'-evidence.json');p=p if p.exists()else folder/(role+'-pending.json');r=json.loads(p.read_text());layout=r['layout'];raw=Path(r['bank']).read_bytes();assert digest_file(r['bank'])==r['bankSha256'];records[role]=r['bankSha256']
    v=np.ndarray((layout['vertexCount'],3),'<f4',buffer=raw,offset=layout['positions']).astype(float)+layout['geometryOrigin'];f=np.ndarray((layout['indexCount']//3,3),'<u4',buffer=raw,offset=layout['indices']);v,inv=np.unique(v,axis=0,return_inverse=True)
    solid=md.Manifold(md.Mesh64(np.ascontiguousarray(v),np.ascontiguousarray(inv[f],dtype=np.uint64)))
    if solid.status()!=md.Error.NoError:raise ValueError(role+': finite wall rejected by independent Boolean library')
    return solid
a,b=load(first),load(second);overlap=a^b
if overlap.status()!=md.Error.NoError:raise ValueError('Finite-wall intersection Boolean failed')
volume=float(overlap.volume());result=dict(first=first,second=second,banks=records,intersectionVolumeMetresCubed=volume,intersectionTriangles=int(overlap.num_tri()),passed=overlap.is_empty()and volume==0.,seconds=time.monotonic()-began,scope='Finite wall occupied volumes, not overlapping lumen envelopes')
if len(sys.argv)>3:
    result['bounds']=list(overlap.bounding_box());print(json.dumps(result));sys.exit(0)
atomic_json(out/(first+'--'+second+'.json'),result);print(json.dumps(result),flush=True)
if not result['passed']:raise ValueError('Distinct hepatic walls overlap')
