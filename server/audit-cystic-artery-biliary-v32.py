"""Read-only, source-bound localization of inherited cystic artery/duct crossings."""
from pathlib import Path
from collections import Counter
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
OUT=WORK/'fixed-artery-biliary-audit-v32'
FRAME='model:local-ta98:assembled-metres'

def main():
    began=time.monotonic();out=OUT/'stage1-localization';assert not out.exists();out.mkdir(parents=True);catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(catalogue.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();artery=reader.read_part(a,parts['TA98CYSTICBLOOD:male:artery'],frame=FRAME,catalogue_path=catalogue);duct=reader.read_part(a,parts['TA98HEPATIC:male:biliary'],frame=FRAME,catalogue_path=catalogue);assert artery.source['bank_sha256']=='d21cc2cd53b1c85713652e791260a323f4c757904969a710628616956ec2ab8b';assert duct.source['bank_sha256']=='b6af93ac5501b21703caa70f5b999510021431957b7bd7a96367e80b414e6588';data,check=pair_contacts(artery,duct,0.,max_pairs=500000,max_tests=5000000);assert check['complete_within_pair_budget'] and check['crossing_face_pairs']==603;atomic_npz(out/'actual-crossing-witnesses.npz',**data);artery.save(out/'fixed-cystic-artery.npz');duct.save(out/'fixed-hepatic-biliary-wall.npz')
    z=np.load(WORK/'owned-parent-cells-v18/artery/partition.npz');ids=z['branch_source_face_ids'];shared=z['shared_branch_faces'];take=np.ones(len(ids),bool);take[shared]=False;order=[]
    for sheet in sorted(set((ids//1000000)[take])):order.extend(np.flatnonzero(take&(ids//1000000==sheet)))
    order.extend(shared);order=np.array(order,int);assert len(order)==len(artery.faces);lineage=ids[order];arterial_sheet=lineage//1000000;station=(lineage%1000000)//36;station[arterial_sheet==6]=-1
    owner=np.full(len(duct.faces),'',object)
    for group in parts['TA98HEPATIC:male:biliary']['tissueDrawGroups']:owner[group['start']//3:(group['start']+group['count'])//3]=group.get('replacement','unowned')
    assert np.all(owner!='');pairs=data['pairs'];points=(data['first']+data['second'])/2;rows=[]
    for key,count in sorted(Counter((int(arterial_sheet[i]),int(station[i]),str(owner[j]))for i,j in pairs).items()):
        mask=np.array([(int(arterial_sheet[i]),int(station[i]),str(owner[j]))==key for i,j in pairs]);q=points[mask];rows.append(dict(arterialSheet={1:'outer',3:'inner',6:'reciprocal-tissue-annulus'}.get(key[0],str(key[0])),sourceStationInterval=key[1],biliaryRegion=key[2],crossingFacePairs=count,witnessBounds_m=[q.min(0).tolist(),q.max(0).tolist()]))
    old=json.loads((WORK/'emitted-context-v20/proof.json').read_text());new=json.loads((WORK/'vein-recipient-emitted-v30/proof.json').read_text());oldnames=[x['name']for x in old['contactChecks']];whole_checked=any('HEPATIC_male_biliary'in x or'HEPATIC:male:biliary'in x for x in oldnames);v20code=ROOT/'server/verify-cystic-emitted-context-v20.py';assert "p['id']=='TA98HEPATIC:male:bile'"in v20code.read_text();selected=[p['id']for p in a['parts']if not p.get('suppressed')and(p['id']in ['FJ3080','TA98PART:FJ3079:A05.8.01.061']or p['id']=='TA98HEPATIC:male:bile')];assert selected==['FJ3080','TA98PART:FJ3079:A05.8.01.061']
    atomic_npz(out/'face-lineage.npz',artery_source_face_ids=lineage,artery_source_stations=station,artery_source_sheets=arterial_sheet,biliary_owner_ids=np.array(owner,str))
    proof=dict(schema='inherited-cystic-artery-biliary-localization/v32',sourceCatalogueSha256=digest_file(catalogue),sourceArteryRecord=parts['TA98CYSTICBLOOD:male:artery'],sourceBiliaryRecord=parts['TA98HEPATIC:male:biliary'],fixedArteryGeometrySha256=artery.geometry_sha,fixedArteryBankSha256=artery.source['bank_sha256'],fixedBiliaryGeometrySha256=duct.geometry_sha,fixedBiliaryBankSha256=duct.source['bank_sha256'],contact=check,exactOffendingArteryFaces=np.unique(pairs[:,0]).tolist(),exactOffendingBiliaryFaces=np.unique(pairs[:,1]).tolist(),witnessBounds_m=[points.min(0).tolist(),points.max(0).tolist()],sourceIntervals=rows,sourceAliasExcludedByDistinctBanksAndStrictNoncoplanarCrossings=True,v20WholeBiliaryActuallyChecked=whole_checked,v20CheckedDuctIds=selected,v20WholeBiliaryIdBug=dict(selected='TA98HEPATIC:male:bile',actual='TA98HEPATIC:male:biliary'),v20SourceArteryBankSha256=old['selectorsAndNormals'][0]['newBankSha256'],v30SourceArteryBankSha256=new['selectorsAndNormals'][0]['newBankSha256'],previousZeroWasOmittedWholeDuctScopeNotChangedGeometry=True,v20ProofSha256=digest_file(WORK/'emitted-context-v20/proof.json'),v30ProofSha256=digest_file(WORK/'vein-recipient-emitted-v30/proof.json'),geometryModified=False,installed=False,elapsedSeconds=time.monotonic()-began);atomic_json(out/'proof.json',proof);print(json.dumps(dict(stage1=str(out),crossingPairs=check['crossing_face_pairs'],arteryFaces=len(proof['exactOffendingArteryFaces']),ductFaces=len(proof['exactOffendingBiliaryFaces']),stationIntervals=sorted(set(r['sourceStationInterval']for r in rows)),regions=sorted(set(r['biliaryRegion']for r in rows)),earlierWholeDuctCheckOmitted=not whole_checked,elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
