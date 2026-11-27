"""Read-only exact original source tooth-contact ancestry comparison."""
from pathlib import Path
import json,time,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1';C=B/'candidate-v4'
def main():
 start=time.monotonic();p=json.loads((B/'plan.json').read_text());cp=ROOT/p['source_catalogue'];cat=json.loads(cp.read_text());a={k:v for k,v in np.load(B/'source-course-arrays.npz').items()};reader=AtlasReader();rows=[]
 for sid,side in[(0,'l'),(1,'r')]:
  src=next(x for x in p['sources']if x['source']==sid);sm=Mesh.load(src['path']);pid=f'ZA:Upper second molar tooth.{side}:Teeth-roots';part=next(x for x in cat['parts']if x['id']==pid);tooth=reader.read_part(cat,part,frame=sm.frame,catalogue_path=cp,catalogue_sha=p['source_master_sha256']);rr=a[f's{sid}_c6_rings'];station={int(v):i for i,ring in enumerate(rr)for v in ring};ids=np.flatnonzero(np.all(np.isin(sm.faces,rr.ravel()),axis=1));sub=Mesh(sm.vertices,sm.faces[ids],sm.frame,sm.attributes,sm.source);pairs,summary=pair_contacts(sub,tooth,0.,max_pairs=100000,max_tests=500000);assert summary['complete_within_pair_budget'];bad=pairs['pairs'][np.isin(pairs['relation'],['crossing','coplanar-overlap'])];cross=ids[bad[:,0]];strips=[min(station[int(v)]for v in sm.faces[f])for f in cross]
  rows.append(dict(source=sid,fixed_id=pid,original_c6_strict_pairs=len(bad),original_source_faces=sorted(set(map(int,cross))),original_c6_contact_strips=sorted(set(strips)),summary=summary,not_proof_that_new_cone_contacts_are_permitted=True))
 report=dict(complete=True,rows=rows,elapsed_seconds=time.monotonic()-start,scope='Original component6 surface versus actual fixed tooth roots; no geometry or catalogue mutation. New sleeve intersections require explicit anatomical disposition regardless of original nearby crossings.');(C/'tooth-native-baseline.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
if __name__=='__main__':main()
