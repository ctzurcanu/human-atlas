"""Independent packed selector, normal-field and local crossing checks."""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts

WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1';OUT=WORK/'vein-recipient-emitted-v30';FRAME='model:local-ta98:assembled-metres'


def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();reader=AtlasReader();rows=[];meshes={}
    for role,number in [('artery',20),('vein',21)]:
        stage=WORK/('owned-region-banks-v19'if role=='artery'else'vein-recipient-bank-v29')/role;manifest=json.loads((stage/'manifest.json').read_text());records=json.loads((stage/'records.json').read_text());z=np.load(WORK/('owned-parent-cells-v18'if role=='artery'else'vein-recipient-partition-v28')/role/'partition.npz');source=Mesh.load(WORK/('fixed-target-'+str(number)+'.npz'));parent_id=manifest['parentId'];branch_id=manifest['branchId'];bank=stage/'cells-and-regions.bin';catalogue={'chunks':[dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')))]}
        for p in records.values():p['chunk']=0
        parent=reader.read_part(catalogue,records[parent_id],frame=FRAME);branch=reader.read_part(catalogue,records[branch_id],frame=FRAME);meshes[role]=dict(parent=parent,branch=branch)
        pids=z['parent_source_face_ids'];owners=z['parent_owner_ids'];shared=set(map(int,z['shared_parent_faces']));take=np.array([i not in shared for i in range(len(pids))]);order=[]
        for owner in sorted(set(owners[take])):
            for sheet in sorted(set((pids//1000000)[take&(owners==owner)])):order.extend(np.flatnonzero(take&(owners==owner)&(pids//1000000==sheet)))
        order.extend(z['shared_parent_faces']);order=np.array(order,int);native=np.isin(pids[order]//1000000,[0,2]);expected=z['parent_native_corner_words'][order[native]];actual=parent.attributes['normals'][parent.faces[native]]*parent.attributes['normalSigns'][parent.faces[native],None];assert np.array_equal(expected,actual)
        literal_faces=0;literal_corners=0
        for emitted_index,source_index in zip(np.flatnonzero(native),order[native]):
            oldface=source.triangles[pids[source_index]%1000000];new=parent.triangles[emitted_index];distance=np.linalg.norm(new[:,None]-oldface[None],axis=2);near=distance.argmin(1);literal=distance.min(1)==0
            literal_corners+=int(literal.sum())
            if literal.all() and len(set(map(int,near)))==3:literal_faces+=1
        region_faces=0;region_checks=[]
        for pid,p in records.items():
            if pid in [parent_id,branch_id]:continue
            m=reader.read_part(catalogue,p,frame=FRAME);region_faces+=len(m.faces);assert all(p[k]==records[parent_id][k]for k in ['positions','normals','normalSigns','vertexCount','geometryOrigin']);region_checks.append(dict(id=pid,triangles=len(m.faces),sameNewAttributeBank=True))
        assert region_faces==len(parent.faces)
        rows.append(dict(role=role,newBankSha256=digest_file(bank),nativeInheritedNormalCorners=int(native.sum())*3,literalSourceFaces=literal_faces,literalSourceCorners=literal_corners,nativeInheritedNormalWordsLiteral=True,allRegionIndicesUseNewParentAttributeBank=True,regionalFacesExactlyCoverParent=region_faces==len(parent.faces),regions=region_checks))
    checks=[]
    def pair(name,a,b):
        data,result=pair_contacts(a,b,0.,max_pairs=500000,max_tests=5000000);atomic_npz(OUT/(name+'.npz'),**data);checks.append(dict(name=name,**result));print(json.dumps(dict(check=name,crossings=result['crossing_face_pairs'])),flush=True)
    pair('new-cystic-artery-vs-new-cystic-vein',meshes['artery']['branch'],meshes['vein']['branch'])
    pair('new-cystic-artery-vs-new-portal-parent',meshes['artery']['branch'],meshes['vein']['parent'])
    pair('new-cystic-vein-vs-new-arterial-parent',meshes['vein']['branch'],meshes['artery']['parent'])
    inner=Mesh.load(WORK/'first-similarity/fixed-mucosa-free-surface.npz')
    for role in ['artery','vein']:pair(role+'-vs-GB-innermost',meshes[role]['branch'],inner)
    source_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';atlas=json.loads(source_path.read_text());active=[p for p in atlas['parts']if not p.get('suppressed') and (p['id']in ['FJ3080','TA98PART:FJ3079:A05.8.01.061'] or p['id']in ['TA98HEPATIC:male:bile','TA98HEPATIC:male:biliary','TA98HEPATIC:male:hepatic_vein'])]
    installed={p['id']:p for p in atlas['parts']};assert digest_file(source_path)=='5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a';assert reader.read_part(atlas,installed['TA98CYSTICBLOOD:male:artery'],frame=FRAME,catalogue_path=source_path).geometry_sha==meshes['artery']['branch'].geometry_sha;assert reader.read_part(atlas,installed['TA98HEPATIC:male:artery'],frame=FRAME,catalogue_path=source_path).geometry_sha==meshes['artery']['parent'].geometry_sha
    for p in active:
        duct=reader.read_part(atlas,p,frame=FRAME,catalogue_path=source_path)
        for role in ['artery','vein']:pair(role+'-vs-'+p['id'].replace(':','_'),meshes[role]['branch'],duct)
    passed=all(r['complete_within_pair_budget']and r['crossing_face_pairs']==0 and r['coplanar_overlap_face_pairs']==0 for r in checks)
    atomic_json(OUT/'proof.json',dict(schema='cystic-recipient-emitted/v30',catalogueSha256=digest_file(source_path),selectorsAndNormals=rows,contactChecks=checks,allRequiredLocalCrossingsZero=passed,atlasModified=False,installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began));print(json.dumps(dict(output=str(OUT/'proof.json'),allLocalCrossingsZero=passed,contactChecks=len(checks),elapsedSeconds=time.monotonic()-began)))


if __name__=='__main__':main()
