"""Check exact fixed venous annuli before attempting another route."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
sp=importlib.util.spec_from_file_location('v14',ROOT/'server/construct-admissible-cystic-departure-v14.py');v14=importlib.util.module_from_spec(sp);sp.loader.exec_module(v14)
WORK=v14.WORK;OUT=WORK/'vein-fixed-artery-v22';CURRENT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='be180f3bff34d92e069160e34999f94a82a07874ff30b6f71771ad5f8cd091c6'


def main():
    began=time.monotonic();assert digest_file(CURRENT)==EXPECTED;assert not OUT.exists();OUT.mkdir();a=json.loads(CURRENT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();artery=reader.read_part(a,parts['TA98CYSTICBLOOD:male:artery'],frame=v14.FRAME,catalogue_path=CURRENT);vein=Mesh.load(WORK/'owned-region-banks-v19/vein/cystic-emitted.npz');z=np.load(WORK/'owned-parent-cells-v18/vein/partition.npz');ids=z['branch_source_face_ids'];shared=z['shared_branch_faces'];not_shared=np.ones(len(ids),bool);not_shared[shared]=False;order=[]
    for role in sorted(set((ids//1000000)[not_shared])):order.extend(np.flatnonzero(not_shared&(ids//1000000==role)))
    order.extend(shared);order=np.array(order);assert len(order)==len(vein.faces);ordered=ids[order];data=np.load(WORK/'cystic-vein-v16/sections.npz');_,_,caps=v14.v10.swept_profiles(data['outer_sections']);terminal=(ordered//1000000==1)&np.isin(ordered%1000000,caps[1]);root=(ordered//1000000==6);rows=[]
    for name,selected in [('literal-external-proximal-annulus',root),('literal-distal-annulus',terminal)]:
        rim=Mesh(vein.vertices,vein.faces[selected],vein.frame,vein.attributes,source=dict(exactAcceptedVeinAnnulus=True,sourceCatalogueSha256=EXPECTED));rim.save(OUT/(name+'.npz'));contacts,proof=pair_contacts(rim,artery,0.,max_pairs=250000,max_tests=2000000);assert proof['complete_within_pair_budget'];atomic_npz(OUT/(name+'-artery-contact.npz'),**contacts);rows.append(dict(name=name,triangles=len(rim.faces),fixedGeometrySha256=rim.geometry_sha,fixedArteryGeometrySha256=artery.geometry_sha,contact=proof,actualWitnessPoints_m=((contacts['first']+contacts['second'])/2).tolist()))
    impossible=any(r['contact']['crossing_face_pairs']or r['contact']['coplanar_overlap_face_pairs']for r in rows)
    result=dict(schema='fixed-vein-endpoint-feasibility/v22',sourceCatalogueSha256=EXPECTED,stage=1,maximumStages=3,fixedArtery='TA98CYSTICBLOOD:male:artery',fixedVeinSourceStage='cystic-vein-v16',fixedParentPortal='TA98HEPATIC:male:portal',fixedRootAndDistalAnnuli=True,annuli=rows,routeConstructionFeasibleWithLiteralAnnuli=not impossible,geometryConstructed=False,atlasModified=False,installed=False,elapsedSeconds=time.monotonic()-began)
    if impossible:result['failure']='At least one literal fixed annulus already crosses the installed fixed artery. Changing only the intermediate route cannot remove that end-face intersection.'
    atomic_json(OUT/'endpoint-feasibility.json',result);print(json.dumps(dict(output=str(OUT/'endpoint-feasibility.json'),feasibleWithLiteralAnnuli=not impossible,endpointCrossings={r['name']:r['contact']['crossing_face_pairs']for r in rows},elapsedSeconds=result['elapsedSeconds'])))


if __name__=='__main__':main()
