"""Actual current preferred native route/bone check; read-only guide review."""
from pathlib import Path
from collections import Counter
import json,sys,time
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
start=time.monotonic();B=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1';O=B/'cranial-route-fit-v1/current-native-route-check-v1';O.mkdir(exist_ok=False);cp=R/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';cat=json.loads(cp.read_bytes());vis={p['id']:p for p in json.loads((B.parent/'circulatory-completion/native-mca-controls-v3/proposal-v2/resolved-installed.json').read_bytes())['effectivePartVisibility']};reader=AtlasReader(root=R);context=B/'regional-interface-v1/fixed-context';bone=Mesh.load(context/'full.npz');a=np.load(context/'full.npz');owner=a['face_owner_part_index'];rows=json.loads((context/'context.json').read_bytes())['partOwnership'];partids=[i for i,r in enumerate(rows)if r['region']=='headneck'or'Vertebra C7'in r['id']];ids=np.flatnonzero(np.isin(owner,partids));bone=Mesh(bone.vertices,bone.faces[ids],'active-male-metres');tree=BVH(bone)
selected=[]
for p in cat['parts']:
 if vis[p['id']]['suppressed']or vis[p['id']]['defaultHidden']:continue
 group=('ICA'if p['id'].startswith('TA98PART:ZA:Internal carotid artery.')else'VA'if p['id'].startswith('TA98PART:ZA:Vertebral artery.')else'Basilar'if p['id']=='ZA:Basilar artery'else None)
 if group:selected.append((p,group))
values=[];faces=[];face_owner=[];offset=0;controls=[]
for p,group in selected:
 m=reader.read_part(cat,p,frame='active-male-metres');values.append(m.vertices);faces.append(m.faces+offset);face_owner.extend([p['id']]*len(m.faces));offset+=len(m.vertices);controls.append(dict(id=p['id'],group=group,bankPath=m.source['bank_path'],bankSha256=m.source['bank_sha256'],record=p,faces=len(m.faces)))
native=Mesh(np.concatenate(values),np.concatenate(faces),'active-male-metres');contacts,summary=pair_contacts(native,bone,0.,second_bvh=tree,max_pairs=100000,max_tests=600000);pairs=contacts['pairs'][contacts['relation']=='crossing'];by_part=Counter(face_owner[int(p[0])]for p in pairs);by_pair=Counter((face_owner[int(p[0])],rows[int(owner[ids[int(p[1])]])]['id'])for p in pairs);np.savez_compressed(O/'contacts.npz',**contacts,second_full_context_face_ids=ids);report=dict(summary=summary,controls=controls,byNativePart=dict(by_part),byNativeAndBone=[dict(nativePart=k[0],bonePart=k[1],strictTrianglePairs=n)for k,n in by_pair.most_common()],qualification='Actual current preferred native source surfaces; strict triangle contacts are not lumen/wall occupancy or anatomical acceptance. No native route was copied as a fit target merely because it exists.',geometryChanged=False,installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start);(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(summary=summary,byNativePart=dict(by_part),elapsedSeconds=report['elapsedSeconds'])))
