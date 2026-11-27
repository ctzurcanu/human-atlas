"""Prepare source-bound local renal catalogue revisions; never deploy."""
import copy,json,gzip
from .common import ROOT,atomic_json,digest_file
from .renal_delivery import DEST,REGIONS,ROLES

SCHEMA='renal-local-layers-v1'

def propose(sex):
    current=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json';before=json.loads(current.read_text());a=copy.deepcopy(before)
    byid={p['id']:p for p in a['parts']};concepts={c['id']:c for c in a['concepts']};oldmarker=copy.deepcopy(a['finiteRenalIntegration']);assert oldmarker['schema']=='renal-local-detail-v420'
    replacement={};newids=[];groups=[];jobs=[]
    limits='Gross renal assembly with finite adventitial, muscular and mucosal collecting coats. Relative coat depths are illustrative. Cortex is a peripheral field; complete cortical columns/zones, microscopic anatomy, vascular tissue fitting and the intramural bladder junction remain unfinished.'
    for side in ['left','right']:
        out=DEST/(sex+'-'+side);r=json.loads((out/'result.json').read_text());assert r['passed']and digest_file(r['bank'])==r['bank_sha256']
        group=copy.deepcopy(next(g for g in oldmarker['assemblies']if g['side']==side));prefix='RENALDETAIL:'+sex+':'+side+':'
        ci=len(a['chunks']);bank=__import__('pathlib').Path(r['bank']);url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'))
        a['chunks'].append(dict(url=url,gzip=url+'.gz',bytes=bank.stat().st_size,gzipBytes=__import__('pathlib').Path(str(bank)+'.gz').stat().st_size))
        for name,packed in r['packedParts'].items():
            region=name.split(':')[0];pid=r['partIds'][name];is_coat=region in REGIONS;template=byid[prefix+region]if is_coat else byid[pid];p=copy.deepcopy(template)
            for key in ['positions','normals','normalSigns','indices','uvs','colors','sourceOffset','geometryOrigin','tissueDrawGroups','uvMethod','suppressed','defaultHidden','replacementParts']:p.pop(key,None)
            p.update(copy.deepcopy(packed));p.update(chunk=ci,renalRepresentation='renal-layered-v1',sourceId='renal-layer-delivery:'+side+':'+name,coverageLimitation=limits,reviewStatus='local-collecting-coats-and-capsule-installed; complete-renal-anatomy-unfinished')
            if is_coat:
                role=name.split(':')[1];index=ROLES.index(role);term=('A08.2.01.'+f'{5+index:03d}')if region=='ureter-wall'else('A08.1.05.'+f'{9+index:03d}')
                label={'adventitia':'Adventitia','muscularis':'Muscular layer','mucosa':'Mucosa'}[role];place={'collecting-wall':'calyces and renal pelvis','upj':'pelviureteric junction','ureter-wall':'ureter'}[region]
                p.update(name=label+' of '+place+' ('+side+')',conceptId=pid,ta98Term=term,tissueWallLayer=role,material='RENAL-LAYER:'+role)
                color={'adventitia':[209,163,110],'muscularis':[156,69,77],'mucosa':[227,135,145]}[role]
                p['tissueAppearance']=dict(color=[x/255 for x in color],pattern='fibrous'if role=='adventitia'else'solid',scale=.00045,relief=0)
                a['materials'][p['material']]=dict(color=color,roughness=.75,metalness=0,surfaceMapping='object-triplanar')
                parent=concepts[template['conceptId']]['partOf'];c=dict(id=pid,name=p['name'],elements=[pid],partOf=parent,ta98Term=term,ta98Kind='exact',ta98Method='stored finite shared-interface coat; relative thickness illustrative',coverageLimitation=limits)
                a['concepts'].append(c);concepts[pid]=c;a.setdefault('ta98TermBindings',[]).append(dict(term=term,conceptId=pid,partIds=[pid],kind='exact'))
                a['parts'].append(p);byid[pid]=p;newids.append(pid)
            else:
                target=byid[pid];target.clear();target.update(p);newids.append(pid)
        for region in REGIONS:
            oldid=prefix+region;ids=[r['partIds'][region+':'+role]for role in ROLES];replacement[oldid]=ids
            byid[oldid].update(suppressed=True,defaultHidden=True,replacementParts=ids,reviewStatus='retained unpartitioned wall; replaced by three finite coats')
        group['partIds']=[r['partIds'][name]for name in r['packedParts']];group['bankSha256']=r['bank_sha256'];group['fit']=r['fit']
        for source in group['sources']:
            source['partIds']=[q for old in source['partIds']for q in replacement.get(old,[old])]
            if source.get('sourcePartId'):byid[source['sourcePartId']]['replacementParts']=source['partIds']
        groups.append(group)
        jobs.append(dict(source=str(out/'arrangement.npz'),surfaces={s['id']:s['array']for s in r['surfaces']},cells=r['cells'],part_ids=r['partIds'],frame='model:'+('local-ta98'if sex=='male'else'local-female-ta98')+':assembled-metres',origin=r['geometryOrigin']))
    for c in a['concepts']:
        c['elements']=list(dict.fromkeys(q for old in c['elements']for q in replacement.get(old,[old])))
        if 'componentElements'in c:c['componentElements']=list(dict.fromkeys(q for old in c['componentElements']for q in replacement.get(old,[old])))
        if any(pid in newids for pid in c['elements'])and c['id'].startswith(('ZA:Kidney.','LOCAL:female:','RENALDETAIL:')):c['coverageLimitation']=limits
    active=set(newids);a['defaultHidden']=sorted((set(a.get('defaultHidden',[]))|set(replacement))-active)
    a['finiteRenalIntegration']=dict(schema=SCHEMA,sex=sex,partIds=newids,assemblies=groups,previousIntegration=oldmarker,collectingWallCoatsPartitioned=True,medullaryPyramidsAndPapillaeBuilt=True,collectingToDistalUreterCommunicates=True,vascularGeometryUnchanged=True,vascularPackingPassed=False,intramuralBladderInterfaceBuilt=False,renalCortex='Peripheral field; source-vertex numerical slivers cleaned within 0.1 micrometre before recipient encoding.',donorSha256=digest_file(DEST/'donor-clean.npz'),anatomicallyAccepted=False)
    a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts'])
    changed=[p['id']for p in a['parts']if p['id']not in {q['id']for q in before['parts']}or p!=next(q for q in before['parts']if q['id']==p['id'])]
    proposal=DEST/(sex+'.json');atomic_json(proposal,a);job=DEST/(sex+'-job.json');atomic_json(job,dict(kind='tissue-pack',parameters=dict(assemblies=jobs)))
    return dict(model='local-ta98'if sex=='male'else'local-female-ta98',proposed_catalogue=str(proposal),expected_catalogue_sha256=digest_file(current),changed_parts=changed,job_spec=str(job))
