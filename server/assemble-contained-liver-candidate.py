"""Assemble verified channel banks and retain literal liver tissue geometry.

Actual owned branch faces are rebuilt from their parent-face allocation. Closed
section exclusions follow the same course field; none are displayed as solids.
"""
from pathlib import Path
import copy,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_bytes,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.hepatic_course_field import HepaticCourseField,HepaticVeinShear
from geometry_tools.algorithms.adaptive_course_transport import refine_course
from geometry_tools.algorithms.transport_slivers import regularize_transport_slivers

base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/'completion-v186';out.mkdir(exist_ok=False)
a=json.loads((base/'completion-v155/atlas-ready.json').read_text());before=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();records={}
assert digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')==a['liverCompletionIntegration']['sourceCatalogueSha256']
containment=json.loads((base/'hepatic-contained-banks-v199/containment.json').read_text());assert containment['passed']
for path in (base/'contained-hepatic-pairs-v198').glob('*.json'):assert json.loads(path.read_text())['passed']
assert len(list((base/'contained-hepatic-pairs-v198').glob('*.json')))==6
def append(data,value):
    data.extend(b'\0'*((-len(data))%4));offset=len(data);data.extend(value);return offset
def chunk(path):
    url='/local-models/'+str(path.relative_to(ROOT/'.local-models'));index=len(a['chunks']);a['chunks'].append(dict(url=url,bytes=path.stat().st_size));return index,url
for role in ['artery','portal','hepatic_vein','biliary']:
    source=base/('hepatic-contained-banks-'+dict(artery='v170',portal='v181',biliary='v181',hepatic_vein='v197')[role]);proof=json.loads((source/(role+'-evidence.json')).read_text());assert proof['passed'];layout=proof['layout'];raw=Path(proof['bank']).read_bytes();assert digest_file(proof['bank'])==proof['bankSha256']
    old=copy.deepcopy(parts['TA98HEPATIC:male:'+role]);oldraw=reader.asset(a['chunks'][old['chunk']]['url']).read_bytes();oldfaces=np.ndarray((old['indexCount']//3,3),'<u4',buffer=oldraw,offset=old['indices'])
    labels=np.empty(len(oldfaces),object);labels[:]=None
    for g in old['tissueDrawGroups']:labels[g['start']//3:(g['start']+g['count'])//3]=g['replacement']
    assert all(isinstance(x,str)for x in labels)
    parent=np.load(source/(role+'-parents.npz'))['parent']if(source/(role+'-parents.npz')).exists()else np.arange(len(oldfaces));owned=labels[parent]
    faces=np.ndarray((layout['indexCount']//3,3),'<u4',buffer=raw,offset=layout['indices']);assert len(faces)==len(parent)
    positions=np.ndarray((layout['vertexCount'],3),'<f4',buffer=raw,offset=layout['positions']).astype(float)+layout['geometryOrigin'];data=bytearray(raw)
    signature=hashlib.sha256(b''.join(raw[layout[k]:layout[k]+layout['vertexCount']*size]for k,size in [('positions',12),('normals',6),('normalSigns',1)])+json.dumps(layout['geometryOrigin']).encode()).hexdigest()
    regions=[parts[pid]for pid in a['hepaticConnectionIntegration']['regionPartIds']if parts[pid]['hepaticConnectionRole']==role]
    for p in regions:
        original=reader.asset(a['chunks'][p['chunk']]['url']).read_bytes();ff=np.ndarray((p['indexCount']//3,3),'<u4',buffer=original,offset=p['indices']);assert np.array_equal(ff,oldfaces[labels==p['id']]),'Source view ownership differs'
        subset=faces[owned==p['id']];assert len(subset);p.update({k:v for k,v in layout.items()if k not in ['indices','indexCount','triangles']});p.update(indices=append(data,subset.astype('<u4').tobytes()),indexCount=subset.size,triangles=len(subset),tissueDrawGroups=[dict(start=0,count=subset.size)])
        p['bounds']=[positions[np.unique(subset)].min(0).tolist(),positions[np.unique(subset)].max(0).tolist()];p['hepaticAttributeSignature']=signature
    whole=parts[old['id']];whole.update(layout);whole['hepaticAttributeSignature']=signature;whole['bounds']=[positions[np.unique(faces)].min(0).tolist(),positions[np.unique(faces)].max(0).tolist()]
    groups=[];start=0
    for i in range(1,len(owned)+1):
        if i==len(owned)or owned[i]!=owned[start]:groups.append(dict(start=start*3,count=(i-start)*3,replacement=str(owned[start])));start=i
    whole['tissueDrawGroups']=groups;path=out/(role+'-views.bin');atomic_bytes(path,bytes(data));ci,url=chunk(path)
    for p in [whole,*regions]:p['chunk']=ci
    a['hepaticConnectionIntegration']['banks'][role]=dict(url=url,sha256=digest_file(path));records[role]=dict(sourceEvidence=str(source/(role+'-evidence.json')),bankSha256=digest_file(path),views=len(regions),triangles=len(faces))

# Retain the entire old liver bank, including every tissue position and word.
sections=a['liverSectionIntegration'];oldpath=reader.asset(sections['bankUrl']);data=bytearray(oldpath.read_bytes());oldsha=sections['bankSha256'];field=HepaticCourseField(before,base/'hepatic-contained-course-v163/field.npz',base/'separated-courses-v105/evidence.json');vein_field=HepaticVeinShear(json.loads((base/'hepatic-vein-clearance-v196/field.json').read_text()));new_exclusions=[]
for row in sections['exclusions']:
    row=copy.deepcopy(row)
    if row['id']!='retained-native-parents':
        course=vein_field if row['id']=='hepatic_vein'else field
        raw=oldpath.read_bytes();origin=np.array(row['origin']);v=np.ndarray((row['vertexCount'],3),'<f4',buffer=raw,offset=row['positions']).astype(float)+origin;faces=np.ndarray((row['indexCount']//3,3),'<u4',buffer=raw,offset=row['indices']).copy()
        v,faces,_,parent,refinement=refine_course(v,faces,np.tile([1.,0.,0.],(len(v),1)),course,deviation_m=1e-5,max_rounds=18)
        mapped=course(v);local=(mapped-origin).astype('<f4');fixed=np.all(mapped==v,axis=1);faces,slivers=regularize_transport_slivers(v,local.astype(float)+origin,faces,fixed)
        _,canonical=np.unique(local,axis=0,return_inverse=True);wf=canonical[faces];faces=faces[(wf[:,0]!=wf[:,1])&(wf[:,1]!=wf[:,2])&(wf[:,2]!=wf[:,0])]
        mesh=Mesh(local.astype(float)+origin,faces,'model:local-ta98:assembled-metres');row.update(vertexCount=len(v),indexCount=faces.size,positions=append(data,local.tobytes()),indices=append(data,faces.astype('<u4').tobytes()),geometrySha256=mesh.geometry_sha)
        mesh.save(out/('exclusion-'+row['id']+'.npz'));atomic_json(out/('exclusion-'+row['id']+'-construction.json'),dict(refinement=refinement,slivers=slivers))
    new_exclusions.append(row)
path=out/'liver-section-bank.bin';atomic_bytes(path,bytes(data));ci,url=chunk(path);sha=digest_file(path)
for pid in a['liverIntegration']['partIds']:parts[pid]['chunk']=ci
for name in ['liverIntegration','liverSectionIntegration']:a[name].update(bankUrl=url,bankSha256=sha)
sections['exclusions']=new_exclusions;a['liverFeatureIntegration']['bankSha256']=sha
for g in a['liverFeatureIntegration']['termGroups']:g['surface']['sourceBankSha256']=sha
a['liverCompletionIntegration']['hostBankSha256']=sha
functional=bytearray((base/'functional-regions-v152/packed/cells.bin').read_bytes());rows=[]
for e in new_exclusions:
    row=copy.deepcopy(e)
    for key,size in [('positions',e['vertexCount']*12),('indices',e['indexCount']*4)]:row[key]=append(functional,data[e[key]:e[key]+size])
    rows.append(row)
fp=out/'functional-section-bank.bin';atomic_bytes(fp,bytes(functional));fci,fu=chunk(fp)
for pid in a['liverCompletionIntegration']['functionalPartIds']:parts[pid]['chunk']=fci
a['liverCompletionIntegration']['functionalSections'].update(bankUrl=fu,bankSha256=digest_file(fp),sourceBankSha256=sha,exclusions=rows)
a['hepaticContainmentIntegration']=dict(schema='male-hepatic-containment/v1',sourceCatalogueSha256=a['liverCompletionIntegration']['sourceCatalogueSha256'],fields=dict(artery=str(base/'hepatic-contained-course-v163/field.npz'),portal=str(base/'hepatic-contained-course-v163/field.npz'),biliary=str(base/'hepatic-contained-course-v163/field.npz'),hepatic_vein=str(base/'hepatic-vein-clearance-v196/field.json')),channels=records,liverExteriorLiteral=True,protectedInterfacesLiteral=True)
limits='Local macroscopic liver with four gross lobes, eight functional territories, a finite capsule and selectable surface/boundary references. Regional extents remain provisional. Six liver terms still need distinct representation; the separate microscopic unit is illustrative.'
a['liverIntegration']['remainingLimits']=limits
for pid in a['liverIntegration']['partIds']:parts[pid]['description']=limits;parts[pid]['coverageLimitation']=limits
a['hepaticConnectionIntegration']['remainingLimits']='Finite communicating channel networks. Intrahepatic containment is measured against the liver; native hilar/caval entry and exit segments are retained. Named regional boundaries remain provisional.'
atomic_json(out/'atlas.json',a);atomic_json(out/'construction.json',dict(channels=records,liverBankPrefixLiteral=True,sourceLiverBankSha256=oldsha,updatedLiverBankSha256=sha,allUnrelatedPartsLiteral=True,geometryAudited=False,installed=False));print(json.dumps(dict(candidate=str(out/'atlas.json'),banks=len(a['chunks']),addedParts=16)))
