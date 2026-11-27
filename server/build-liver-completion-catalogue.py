"""Assemble reviewed stages into one local candidate; never install here."""
from pathlib import Path
from collections import defaultdict
import copy,sys,json,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra,connected_components
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader

started_at=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'completion-v155';out.mkdir(exist_ok=False)
current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';before=json.loads(current.read_text())
a=json.loads((base/'liver-surface-controls-v143/with-landmark.json').read_text());assert a['parts']==before['parts']and a['chunks']==before['chunks'];a.pop('liverLandmarkIntegration',None)
a['concepts']=[c for c in a['concepts']if c['id']!='TA98LIVER:male:landmark:A05.8.01.016']
inventory=json.loads((base/'male-liver-section-install-v142/inventory.json').read_text());terms={t['authority']['id']:t['authority']for t in inventory['terms']}
parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';host=parts['TA98PART:ZA:Liver:capsule'];m=reader.read_part(a,host,frame=frame)
outer=m.triangles[:host['tissueDrawGroups'][0]['count']//3];v,inv=np.unique(outer.reshape(-1,3),axis=0,return_inverse=True);f=inv.reshape(-1,3);centers=outer.mean(1)
face_for_vertex={}
for i,tri in enumerate(f):
    for k,j in enumerate(tri):face_for_vertex.setdefault(int(j),(i,k))
groups=[];new_ids=[];functional=[];micro=[];ligaments=[]
def term(number):return f'A05.8.01.{number:03d}'
def title(number):return terms[term(number)]['name_en'].capitalize()
def add_concept(number,elements,representation,**extra):
    cid='TA98LIVER:male:term:'+term(number);name=title(number)
    a['concepts'].append(dict(id=cid,name=name,elements=elements,ta98Term=term(number),ta98Kind='exact',partOf='ZA:Liver',description=extra.pop('description',None)or'Local anatomical representation; detailed regional boundaries remain provisional.'))
    group=dict(term=term(number),conceptId=cid,partIds=elements,representation=representation,**extra);groups.append(group);return group
def load_pack(folder,metadata,which=None):
    spec=json.loads((folder/'packed/tissue.json').read_text());assert spec['passed'],str(folder)
    bank=Path(spec['bank']);url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));chunk=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=bank.stat().st_size))
    for key,row in spec['parts'].items():
        if which is not None and key not in which:continue
        number,name,category=metadata[key];pid=row['id'];p={**copy.deepcopy(row),'chunk':chunk,'conceptId':pid,'sourceId':'Liver:'+key,'name':name,'system':'digestive',
            'groups':['8: Visceral systems','Digestive system','Abdomen','Trunk','Liver'],
            'ta98Term':term(number),'pairedTissue':True,'sectionSolid':True,'tissueWallAssembly':'liver:male:'+category,'tissueWallLayer':key,
            'defaultHidden':True,'license':parts['TA98PART:ZA:Liver:right']['license'],
            'material':'LIVER:completion:'+key,'coverageLimitation':'Local representation; anatomical boundaries remain provisional.'}
        colors={'I':[125,75,56],'II':[158,90,67],'III':[145,68,55],'IV':[170,104,70],'V':[132,59,43],'VI':[158,78,59],'VII':[120,73,65],'VIII':[177,91,63],
            'venosum':[216,192,161],'teres':[208,172,139],'lobule':[179,103,76],'arteries':[182,54,58],'portal-venules':[94,127,169],'central-vein':[81,117,151],'bile-ductules':[156,174,87],'perivascular-sheath':[217,197,162]}
        color=colors[key];a.setdefault('materials',{})[p['material']]=dict(color=color,roughness=.62,metalness=0)
        p['tissueAppearance']=dict(color=[c/255 for c in color],pattern='fibrous'if category=='ligaments'else'solid',scale=.0003,relief=0)
        if category=='micro':
            p['geometryOrigin']=[-.068,1.27,.05];p['bounds']=(np.array(p['bounds'])+p['geometryOrigin']).tolist();micro.append(pid)
            p['coverageLimitation']='Representative magnified human-liver unit, also stored here at illustrative microscopic scale; not patient-specific histology.'
        elif category=='functional':functional.append(pid);p['liverFunctionalTerritory']=key
        else:ligaments.append(pid)
        a['parts'].append(p);parts[pid]=p;new_ids.append(pid)
    return spec
segment_numbers={'I':44,'II':39,'III':40,'IV':42,'V':47,'VI':50,'VII':51,'VIII':48}
functional_spec=load_pack(base/'functional-regions-v152',{k:(n,title(n),'functional')for k,n in segment_numbers.items()},set(segment_numbers))
for key,n in segment_numbers.items():add_concept(n,[functional_spec['part_ids'][key]],'functional-territory')
unions={33:list(segment_numbers),37:['I','II','III','IV'],38:['II','III'],41:['IV'],43:['I'],45:['V','VI','VII','VIII'],46:['V','VIII'],49:['VI','VII']}
for n,keys in unions.items():add_concept(n,[functional_spec['part_ids'][k]for k in keys],'functional-territory-union')
lig_spec=load_pack(base/'ligament-courses-v153',{'venosum':(11,title(11),'ligaments'),'teres':(15,title(15),'ligaments')})
for key,n in [('venosum',11),('teres',15)]:add_concept(n,[lig_spec['part_ids'][key]],'fibrous-cord',description='Illustrative fibrous cord following the source-bound fissure course. Attachment extent and calibre remain provisional.')
micro_terms={'perivascular-sheath':55,'lobule':56,'arteries':57,'portal-venules':58,'central-vein':59,'bile-ductules':60}
micro_spec=load_pack(base/'microscopic-unit-v151',{key:(n,title(n),'micro')for key,n in micro_terms.items()})
for key,n in micro_terms.items():add_concept(n,[micro_spec['part_ids'][key]],'illustrative-microscopic-detail',detailModel='liver-unit',description='Illustrative three-dimensional human-liver microscopic unit: radiating anastomosing plates, sinusoidal spaces and portal tracts. Dimensions are construction parameters, not patient measurements.')

# Source-boundary references use actual facet coordinates and barycentrics.
def witness(face,weights):return dict(partId=host['id'],face=int(face),barycentric=list(map(float,weights)),point=(np.array(weights)@outer[face]).tolist())
def vertex_witness(i):
    face,k=face_for_vertex[int(i)];weights=np.eye(3)[k];return witness(face,weights)
edge_owners=defaultdict(list)
for i,tri in enumerate(f):
    for x,y in zip(tri,np.roll(tri,-1)):edge_owners[tuple(sorted((int(x),int(y))))].append((i,int(x),int(y)))
edges=np.array(list(edge_owners),int);norm=np.cross(outer[:,1]-outer[:,0],outer[:,2]-outer[:,0]);norm/=np.linalg.norm(norm,axis=1)[:,None]
side=np.argmax(norm@np.array([[0,1,0],[0,0,1],[-1,0,0],[0,0,-1],[0,-1,0]]).T,axis=1)
def longest_chain(chosen):
    selected=np.array(chosen,int).reshape(-1,2);weights=np.linalg.norm(v[selected[:,0]]-v[selected[:,1]],axis=1)
    graph=coo_matrix((np.r_[weights,weights],(np.r_[selected[:,0],selected[:,1]],np.r_[selected[:,1],selected[:,0]])),shape=(len(v),len(v))).tocsr()
    _,labels=connected_components(graph,directed=False);used=np.unique(selected);label=max(set(labels[used]),key=lambda k:sum(labels[used]==k));nodes=used[labels[used]==label]
    start=int(nodes[np.argmin(v[nodes,0])]);end=int(nodes[np.argmax(v[nodes,0])]);distance,previous=dijkstra(graph,indices=start,return_predecessors=True,directed=False)
    if start==end:end=int(nodes[np.argmax(distance[nodes])])
    path=[end]
    while path[-1]!=start:path.append(int(previous[path[-1]]));assert path[-1]>=0
    return [vertex_witness(i)for i in path[::-1]]
front_border=[key for key,uses in edge_owners.items()if len(uses)==2 and{side[u[0]]for u in uses}=={1,4}]
border=longest_chain(front_border)
add_concept(24,[host['id']],'boundary-reference',curves=[dict(closed=False,witnesses=border)],anchor=border[len(border)//2],description='Source-bound inferior margin between the anterior and visceral surface charts; its regional boundary is provisional.')
teres_surface=next(g for g in a['liverFeatureIntegration']['termGroups']if g['term']==term(14));target=np.mean(outer[teres_surface['surface']['faces']],axis=(0,1))
index=int(np.argmin([np.linalg.norm(np.array(w['point'])-target)for w in border]));notch=border[max(0,index-2):min(len(border),index+3)]
assert len(notch)>=2;add_concept(25,[host['id']],'boundary-reference',curves=[dict(closed=False,witnesses=notch)],anchor=notch[len(notch)//2])
follow={u[1]:u[2]for uses in edge_owners.values()if len(uses)==1 for u in uses};todo=set(follow);curves=[]
while todo:
    start=min(todo);path=[];i=start
    while i in todo:path.append(i);todo.remove(i);i=follow[i]
    assert i==start;curves.append(dict(closed=True,witnesses=[vertex_witness(j)for j in path]))
assert len(curves)==3
anchor=curves[0]['witnesses'][0]
add_concept(16,[host['id']],'opening-configuration',curves=curves,anchor=anchor,
    description='Actual three hilar source openings and their rims. Their combined clinical porta extent remains provisional; no solid closes an opening.')
planes=json.loads((base/'functional-planes-v149/evidence.json').read_text())['planes']
for number,key in [(34,'umbilical'),(35,'main'),(36,'right')]:
    plane=next(q for q in planes if q['name']==key);normal=np.array(plane['normal']);offset=plane['offset'];curves=[]
    for face,tri in enumerate(outer):
        d=tri@normal-offset;ws=[]
        for k in range(3):
            j=(k+1)%3
            if d[k]*d[j]<0:
                s=d[k]/(d[k]-d[j]);w=np.zeros(3);w[k]=1-s;w[j]=s;ws.append(witness(face,w))
        if len(ws)==2:curves.append(dict(closed=False,witnesses=ws))
    assert curves;add_concept(number,[host['id']],'functional-boundary-reference',curves=curves,anchor=curves[len(curves)//2]['witnesses'][0],
        description='Intersection of the source-landmark functional portal plane with the repaired liver surface; a boundary reference, not a physical groove.')

# The omental tuberosity is an existing medial left-lobe prominence. Its
# footprint is selected on actual tissue, not represented by a pin or sphere.
ven_surface=next(g for g in a['liverFeatureIntegration']['termGroups']if g['term']==term(10));ven_center=np.mean(outer[ven_surface['surface']['faces']],axis=(0,1))
target=ven_center+[.013,.010,.001];distance=np.linalg.norm(centers-target,axis=1);valid=(centers[:,0]>ven_center[0])&(norm[:,1]<.7);seed=int(np.argmin(np.where(valid,distance,np.inf)))
neighbors=[[]for _ in outer]
for uses in edge_owners.values():
    if len(uses)==2:x,y=uses[0][0],uses[1][0];neighbors[x].append(y);neighbors[y].append(x)
todo=[seed];chosen=set()
while todo:
    face=todo.pop()
    if face in chosen or not valid[face]or np.linalg.norm(centers[face]-centers[seed])>.012:continue
    chosen.add(face);todo.extend(neighbors[face])
assert len(chosen)>2;face=min(chosen,key=lambda i:np.linalg.norm(centers[i]-np.mean(centers[list(chosen)],axis=0)))
surface_group=dict(term=term(17),conceptId='TA98LIVER:male:surface:'+term(17),partIds=[host['id']],representation='surface-region',surface=dict(partId=host['id'],faces=sorted(chosen),sourceBankSha256=a['liverIntegration']['bankSha256']),anchor=witness(face,[1/3]*3),boundaryProvisional=True)
a['liverFeatureIntegration']['termGroups'].append(surface_group)
a['concepts'].append(dict(id=surface_group['conceptId'],name=title(17),ta98Term=term(17),ta98Kind='exact',elements=[host['id']],partOf='ZA:Liver',description='Source-bound footprint of the existing medial left-lobe prominence beside the venous fissure. Regional extent remains provisional.'))

for pid in a['liverIntegration']['partIds']:
    if pid.endswith(':capsule'):continue
    parts[pid]['liverAlternativeParts']=functional
a['defaultHidden']=sorted(set(a.get('defaultHidden',[]))|set(new_ids))
a['liverCompletionIntegration']=dict(schema='male-liver-completion/v1',sourceCatalogueSha256=digest_file(current),
    partIds=new_ids,functionalPartIds=functional,microscopicPartIds=micro,ligamentPartIds=ligaments,
    hostPartId=host['id'],hostBankSha256=a['liverIntegration']['bankSha256'],termGroups=groups,
    numericalEvidence=[str(base/'functional-regions-v152/evidence.json'),str(base/'microscopic-unit-v151/evidence.json'),str(base/'ligament-courses-v153/evidence.json')],
    regionalBoundariesProvisional=True,microscopicUnitIllustrative=True)
represented=set(g['term']for g in a['liverIntegration']['termGroups'])|{g['term']for g in a['hepaticConnectionIntegration']['termGroups']if g['term'].startswith('A05.8.01.')}
represented|={g['term']for g in a['liverFeatureIntegration']['termGroups']}|{g['term']for g in groups}
missing=sorted(set(terms)-represented);assert len(represented)==63,(len(represented),missing)
ledger=dict(authoritativeTerms=69,representedTerms=63,percentage=round(100*63/69,2),organScaleTerms=57,illustrativeMicroscopicTerms=6,
    missing=[dict(term=t,name=terms[t]['name_en'])for t in missing],represented=sorted(represented),anatomicallyComplete=False,
    sourceCatalogueSha256=digest_file(current),allNewPartsInitiallyHidden=True,elapsedSeconds=time.monotonic()-started_at)
atomic_json(out/'coverage.json',ledger);atomic_json(out/'atlas.json',a)

# A standard local model opens the same illustrative microstructure at 1000×.
micro_rows=[copy.deepcopy(p)for p in a['parts']if p['id']in micro];bank=Path(micro_spec['bank']);data=bytearray(bank.read_bytes());oldrow=next(iter(micro_spec['parts'].values()));n=oldrow['vertexCount'];offset=oldrow['positions'];pos=np.ndarray((n,3),'<f4',buffer=data,offset=offset);pos[:]*=1000
microbank=out/'microscopic-magnified.bin';atomic_bytes(microbank,bytes(data));micro_url='/local-models/'+str(microbank.relative_to(ROOT/'.local-models'))
for p in micro_rows:
    p['chunk']=0;p['geometryOrigin']=[0,0,0];p['bounds']=(np.array(micro_spec['parts'][p['tissueWallLayer']]['bounds'])*1000).tolist();p['defaultHidden']=False;p['groups']=['Liver microscopic unit'];p.pop('liverAlternativeParts',None)
unit=dict(version='Human liver · illustrative microscopic unit',scope='Illustrative human architecture; magnification 1000×; dimensions not patient-specific',
    parts=micro_rows,concepts=[dict(id='TA98LIVER:micro:whole',name='Human liver microscopic unit',elements=micro)],
    chunks=[dict(url=micro_url,bytes=microbank.stat().st_size)],materials=a['materials'],defaultHidden=[],triangles=sum(p['triangles']for p in micro_rows))
for g in groups:
    if g['representation']=='illustrative-microscopic-detail':unit['concepts'].append(dict(id=g['conceptId'],name=terms[g['term']]['name_en'].capitalize(),elements=g['partIds'],ta98Term=g['term'],ta98Kind='exact',partOf='TA98LIVER:micro:whole',description='Illustrative human liver microscopic unit, magnified 1000×.'))
atomic_json(out/'microscopic-atlas.json',unit)
print(json.dumps(ledger),flush=True)
