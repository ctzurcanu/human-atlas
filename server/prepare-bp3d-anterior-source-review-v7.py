"""Review only: the licensed BP3D anterior assembly in its own native frame."""
from pathlib import Path
import sys,time,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

CAT=ROOT/'public/models/atlas.json';SHA='ed04ee235eb489c60126c56b438449507722f44ab15d458ae1ca513855ab1af1'
FRAME='source:BodyParts3D-4.0:packed-native-metres-y-up'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1'
OUT=BASE/'bp3d-anterior-source-review-v7'

def main():
    began=time.monotonic();assert digest_file(CAT)==SHA;OUT.mkdir(exist_ok=False)
    a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();loaded={};ancestry=[]
    def entry(pid,color,lower_y=None):
        if pid not in loaded:loaded[pid]=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=SHA)
        m=loaded[pid];ids=np.arange(len(m.faces))if lower_y is None else np.flatnonzero(np.all(m.triangles[:,:,1]>=lower_y,axis=1));assert len(ids)
        suffix=''if lower_y is None else'-literal-root-faces';path=OUT/(pid+suffix+'.npz')
        q=Mesh(m.vertices,m.faces[ids],FRAME,m.attributes,dict(m.source,diagnostic_source_face_ids=ids.tolist(),no_source_change=True,uncapped_literal_face_selection=lower_y is not None));q.save(path)
        ancestry.append(dict(id=pid,originalName=parts[pid]['name'],originalFmaId=parts[pid]['conceptId'],licence='CC BY 4.0',sourceCatalogueSha256=SHA,sourceBankSha256=m.source['bank_sha256'],sourceGeometrySha256=m.geometry_sha,sourceAttributeSha256=m.attributes_sha,sourceFaceIds=ids.tolist(),snapshot=str(path),snapshotSha256=digest_file(path),frame=FRAME,noTransformApplied=True))
        return dict(mesh=str(path),name=pid+' '+parts[pid]['name']+suffix,color=color)
    vessels=[entry(pid,c)for pid,c in [('FJ1654',[.98,.50,.10,1]),('FJ1654M',[.90,.12,.12,1]),('FJ1655',[.78,.08,.68,1]),('FJ1682',[.55,.30,.08,1]),('FJ1682M',[.55,.08,.08,1])]]
    chiasm=[entry(pid,[.95,.86,.15,1])for pid in ['FJ1771','FJ1818']]
    brain=[entry(pid,[.56,.57,.58,1])for pid in ['FJ1758','FJ1806']]
    sphenoid=entry('FJ3394',[.66,.74,.80,1])
    roots=[entry(pid,c,1.597)for pid,c in [('FJ1682',[.55,.30,.08,1]),('FJ1682M',[.55,.08,.08,1])]]
    comm=loaded['FJ1655'];mid=comm.vertices.mean(0)
    common=dict(frame=FRAME,resolution=1100,wall_seconds=60)
    specs=[('whole-native-context',dict(**common,meshes=vessels+chiasm+brain+[sphenoid],views=[dict(name='BP3D complete anterior source and fixed BP3D brain/chiasm/sphenoid; inferior',direction=[0,-1,.12],up=[0,0,1]),dict(name='BP3D complete anterior source; anterior context',direction=[.3,-.3,1]),dict(name='BP3D complete anterior source; lateral context',direction=[1,0,.1])])),('focused-native-junction',dict(**common,meshes=vessels[:3]+roots+chiasm,views=[dict(name='BP3D ACA-ACom-ICA roots; inferior source course',direction=[0,-1,0],up=[0,0,1]),dict(name='BP3D ACA-ACom-ICA roots; anterior source course',direction=[0,0,1]),dict(name='BP3D native communicating overlap; axial uncapped section',direction=[0,1,0],up=[0,0,1],section=dict(normal=[0,1,0],point=mid.tolist(),clear_outer=True)),dict(name='BP3D native communicating overlap; coronal uncapped section',direction=[0,0,1],section=dict(normal=[0,0,1],point=mid.tolist(),clear_outer=True))]))]
    queued=[]
    for name,spec in specs:
        atomic_json(OUT/(name+'-spec.json'),spec);q=prepare(None,spec);q['purpose']=name;queued.append(q)
    atomic_json(OUT/'source-face-provenance.json',ancestry)
    proof=dict(schema='bp3d-native-anterior-source-review/v7',sourceCatalogueSha256=SHA,frame=FRAME,originalPublisherProvenance=a['provenance'],noSourceGeometryChanged=True,noRegistrationToTa98AppliedOrAssumed=True,allContextsFromExactSameNativeCatalogue=True,sourceSolidCastsDoNotEstablishPairedWallsOrPatency=True,queue=queued,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',proof);assert digest_file(CAT)==SHA
    print(json.dumps(dict(output=str(OUT/'proof.json'),queue=[dict(id=q['request_id'],purpose=q['purpose'])for q in queued],elapsedSeconds=proof['elapsedSeconds'])))

if __name__=='__main__':main()
