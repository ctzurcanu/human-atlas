"""Selective actual component/recipient tests; no new vascular identities."""
from pathlib import Path
import sys,time,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/upper-limb-native-network-opportunities-v13';OUT=BASE/'junction-and-landmark-replay-v14';CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';FRAME='model:local-ta98:assembled-metres'

def serial(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);prior=json.loads((BASE/'proof.json').read_text());assert digest_file(CAT)==prior['sourceCatalogueSha256'];a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};source={r['partId']:r for r in prior['sourceMeshes']};reader=AtlasReader();meshes={};trees={};ancestry={};tests=[];rims=[]
    def get(pid):
        if pid not in meshes:meshes[pid]=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT)
        return meshes[pid]
    def component(pid,i):
        key=pid+'#component-'+str(i)
        if key not in meshes:
            m=get(pid);ids=np.array(source[pid]['components'][i]['sourceFaceIds']);used,inv=np.unique(m.faces[ids],return_inverse=True);meshes[key]=Mesh(m.vertices[used],inv.reshape(-1,3),FRAME,{k:v[used]for k,v in m.attributes.items()},dict(m.source,sourceFaceIds=ids.tolist(),sourceVertexIds=used.tolist(),diagnosticComponentOnly=True));ancestry[key]=dict(partId=pid,sourceFaceIds=ids,sourceVertexIds=used,geometrySha256=m.geometry_sha,bankSha256=m.source['bank_sha256'])
        return key
    def tree(pid):
        if pid not in trees:trees[pid]=BVH(meshes[pid])
        return trees[pid]
    def test(x,y,contacts=True):
        closest=tree(x).closest_mesh(tree(y));closest.update(firstTriangle_m=meshes[x].triangles[closest['first_face_id']],secondTriangle_m=meshes[y].triangles[closest['second_face_id']]);row=dict(first=x,second=y,closest=closest)
        if contacts:
            d,gate=pair_contacts(meshes[x],meshes[y],0.,first_bvh=tree(x),second_bvh=tree(y),max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget'];row.update(actualContactGate=gate,actualContactFacePairs=d['pairs'])
        tests.append(row);return row
    for side in ['l','r']:
        network='O3M:Superficial veins of upper limb.'+side;course=component(network,0)
        targets=['ZA:Cephalic vein.'+side,'ZA:Axillary vein.'+side,'O3M:Subclavian vein.'+side,'ZA:Clavicle.'+side+':Bone','ZA:Clavicular part of deltoid muscle.'+side+':Abductor','ZA:Acromial part of deltoid muscle.'+side+':Abductor','ZA:Clavicular head of pectoralis major muscle.'+side+':Adductor','ZA:Pectoralis minor muscle.'+side]
        for target in targets:
            if target not in parts:continue
            get(target);test(course,target)
        for i in range(1,len(source[network]['components'])):test(course,component(network,i),False)
        for i,q in enumerate(boundary_loops(meshes[course].faces)):
            points=meshes[course].vertices[q];rid=course+'#actual-rim-'+str(i);meshes[rid]=Mesh(points,np.c_[np.arange(len(points)),np.roll(np.arange(len(points)),-1),np.arange(len(points))],FRAME,{},dict(diagnosticLineQueriesOnly=True));rr=dict(component=course,actualSourceLoop=i,sourceComponentVertexIds=q,points_m=points,centroid_m=points.mean(0),notAPairedWallAnnulus=True,actualRecipientTests=[])
            for target in targets[:3]:rr['actualRecipientTests'].append(test(rid,target,False))
            rims.append(rr)
        acromial='O3M:Thoracoacromial artery Acromial br.'+side;c0=component(acromial,0);c1=component(acromial,1);test(c0,c1)
        for curve in [c0,c1]:
            for target in ['ZA:Suprascapular artery.'+side,'ZA:Posterior circumflex humeral artery.'+side,'ZA:Clavicle.'+side+':Bone','ZA:Scapula.'+side+':Bone']:
                get(target);test(curve,target)
    result=dict(schema='upper-limb-literal-junction-feasibility/v14',sourceCatalogueSha256=prior['sourceCatalogueSha256'],priorSourceInventorySha256=digest_file(BASE/'proof.json'),sourceComponentAncestry=ancestry,tests=tests,actualVenousComponentRims=rims,readyNewTerms=[],feasibility='No identity or flow continuity is inferred merely from a nearby named artery, coarse whole-scapula surface, or a generic network carrier. Actual named branch/recipient and anatomically identified coracoid/acromion boundaries remain required.',noGeometryOrLiveMetadataChanged=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',json.loads(json.dumps(result,default=serial)));assert digest_file(CAT)==prior['sourceCatalogueSha256']
    print(json.dumps(dict(output=str(OUT/'proof.json'),selectedTests=[dict(first=r['first'],second=r['second'],gap_mm=1000*r['closest']['distance_m'],crossingPairs=r.get('actualContactGate',{}).get('crossing_face_pairs'))for r in tests if 'actual-rim'not in r['first']and('#component'not in r['second']or('Acromial br' in r['first']and'Acromial br'in r['second']))],elapsedSeconds=result['elapsedSeconds'])))
if __name__=='__main__':main()
