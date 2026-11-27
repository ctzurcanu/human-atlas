"""Prepare a reviewed local urinary revision with literal retained kidney tissue.

Preparation never changes the atlas. Publication uses the guarded installer,
with its source inventory, mesh audits, review and reversible transaction.
"""
from pathlib import Path
import copy,gzip,json,sqlite3
import numpy as np

from .common import ROOT,atomic_bytes,atomic_json,atomic_npz,digest_file
from .model import AtlasReader,Mesh


GEOMETRY_KEYS=('positions','normals','normalSigns','indices','uvs','colors','sourceOffset','geometryOrigin',
               'tissueDrawGroups','vertexCount','indexCount','bounds')
BLADDER_TERMS={'serosa':'008','subserosa':'009','detrusor-unstratified':'015','trigone-superficial':'012',
 'trigone-deep':'013','neck-external-longitudinal':'017','neck-circular':'018','neck-internal-longitudinal':'019',
 'submucosa':'022','mucosa':'023','mucosa-trigone':'023','mucosa-interureteric-crest':'025','mucosa-uvula':'006',
 'mucosa-ureteric-orifice-left':'026','mucosa-ureteric-orifice-right':'026','mucosa-internal-urethral-orifice':'027M'}
LABELS={'serosa':'Serosa','subserosa':'Subserosa','adventitia':'Adventitia','muscularis':'Muscular coat',
 'detrusor-unstratified':'Unstratified detrusor','trigone-superficial':'Superficial trigonal muscle',
 'trigone-deep':'Deep trigonal muscle','neck-external-longitudinal':'External longitudinal muscle',
 'neck-circular':'Circular muscle','neck-internal-longitudinal':'Internal longitudinal muscle',
 'submucosa':'Submucosa','mucosa':'Mucosa','mucosa-trigone':'Trigonal mucosa',
 'mucosa-interureteric-crest':'Interureteric crest','mucosa-uvula':'Uvula',
 'mucosa-ureteric-orifice-left':'Ureteric ostial mucosa',
 'mucosa-ureteric-orifice-right':'Ureteric ostial mucosa',
 'mucosa-internal-urethral-orifice':'Internal urethral ostial mucosa'}
LIMITS='Finite wall with shared urinary joins. Coat depths and regional boundaries remain provisional. Ligament, muscular attachments, folds and surface refinement remain unfinished.'


def bladder_term_groups(part_ids):
    """Named surface features remain part of their actual mucosal tissue."""
    ids=list(part_ids);groups={'A08.3.01.001':ids}
    region=lambda p:p.split(':')[2]
    role=lambda p:p.rsplit(':',1)[-1]
    for name,code in [('body','003'),('fundus','004'),('neck','005')]:
        groups['A08.3.01.'+code]=[p for p in ids if region(p)==name]
    for name,code in BLADDER_TERMS.items():
        if name in ('mucosa','mucosa-trigone'):continue
        members=[p for p in ids if role(p)==name]
        if members:groups.setdefault('A08.3.01.'+code,[]).extend(members)
    groups['A08.3.01.023']=[p for p in ids if role(p).startswith('mucosa')]
    muscle=[p for p in ids if role(p).startswith(('detrusor-','trigone-','neck-'))]
    groups['A08.3.01.010']=muscle
    groups['A08.3.01.011']=[p for p in muscle if role(p).startswith('trigone-')]
    groups['A08.3.01.016']=[p for p in muscle if region(p)=='neck']
    groups['A08.3.01.014']=[p for p in muscle if not role(p).startswith('trigone-')]
    groups['A08.3.01.024']=[p for p in ids if role(p).startswith('trigone-')or role(p)in
        ('mucosa-trigone','mucosa-interureteric-crest','mucosa-uvula')or 'orifice' in role(p)]
    if any(not members for members in groups.values()):raise ValueError('A declared wall group has no finite tissue')
    return groups


def prepare_job(delivery_dir,node_dir,output):
    delivery_dir,output=map(lambda p:Path(p).resolve(),(delivery_dir,output));node_dir=Path(node_dir).resolve()if node_dir else None
    if output.exists():raise FileExistsError(output)
    d=json.loads((delivery_dir/'delivery.json').read_text());nv=json.loads((node_dir/'verification.json').read_text())if node_dir else None
    forest=json.loads((node_dir/'forest-fit-check.json').read_text())if node_dir else None
    if not d.get('passed')or node_dir and(not nv.get('passed')or not forest.get('all_passed')):
        raise ValueError('Complete numerical and dependency checks are required before job preparation')
    proof=d.get('whole_geometry_verification',{})
    if digest_file(proof['path'])!=proof['sha256']or not proof.get('actual_indexed_geometry_unchanged'):
        raise ValueError('The decoded whole-wall proof changed')
    if digest_file(d['catalogue'])!=d['catalogue_sha256']or nv and nv['catalogue_sha256']!=d['catalogue_sha256']:
        raise ValueError('The source atlas changed')
    output.mkdir(parents=True);assemblies=[];names=[]
    for name,chart in d['groups'].items():
        meta=json.loads((delivery_dir/name/'delivery.json').read_text())
        if digest_file(meta['bank'])!=chart['bank_sha256']:raise ValueError('A urinary bank changed')
        cells=meta['cells'];used={u['surface']for uses in cells.values()for u in uses}
        assembly=dict(source=d['source'],surfaces={s:d['surfaces'][s]for s in sorted(used)},
            corner_words={s:d['corner_words'][s]for s in sorted(used)},cells=cells,
            part_ids={pid:pid for pid in cells},frame=d['frame'],origin=chart['origin'])
        if name=='prostate':assembly['corner_uvs']={s:d['corner_uvs'][s]for s in sorted(used)}
        assemblies.append(assembly);names.append(name)
    nodes=json.loads((node_dir/'candidates.json').read_text())if node_dir else {'rows':[]}
    for i,row in enumerate(nodes['rows']):
        pid=row['part_id'];mesh=Mesh.load(row['path']);path=output/('node-'+str(i)+'.npz')
        words=mesh.attributes['normals'][mesh.faces].astype(np.int32)
        if 'normalSigns'in mesh.attributes:words*=mesh.attributes['normalSigns'][mesh.faces,None]
        atomic_npz(path,positions=mesh.vertices,faces=mesh.faces,words=words.astype('<i2'),uvs=mesh.attributes['uvs'][mesh.faces])
        assemblies.append(dict(source=str(path),surfaces={'node':'faces'},corner_words={'node':'words'},corner_uvs={'node':'uvs'},
            cells={pid:[dict(surface='node',side=1,neighbor=None)]},part_ids={pid:pid},frame=mesh.frame,origin=row['translation_m']))
        names.append('node:'+pid)
    parameters=dict(assemblies=assemblies,source_catalogue=d['catalogue'],whole_geometry_verification=proof['path'],
        delivery_manifest=str(delivery_dir/'delivery.json'))
    if d.get('feature_plan'):parameters['feature_plan']=d['feature_plan']
    if node_dir:parameters.update(node_verification=str(node_dir/'verification.json'),node_dependencies=str(node_dir/'forest-fit-check.json'))
    spec=dict(kind='tissue-pack',parameters=parameters)
    atomic_json(output/'job.json',spec)
    atomic_json(output/'preparation.json',dict(delivery_dir=str(delivery_dir),node_dir=str(node_dir)if node_dir else None,groups=names,
        catalogue=d['catalogue'],catalogue_sha256=d['catalogue_sha256'],job_spec=str(output/'job.json'),installed=False))
    return dict(job_spec=str(output/'job.json'),assemblies=len(assemblies))


def appearance(role):
    color=([.83,.57,.55]if role=='serosa'else[.82,.67,.44]if role=='subserosa'else
           [.8,.61,.42]if role=='adventitia'else[.87,.71,.57]if role=='submucosa'else
           [.86,.52,.55]if'mucosa'in role else[.69,.32,.37])
    smooth_mucosa=role in ('mucosa-trigone','mucosa-interureteric-crest','mucosa-uvula')or 'orifice' in role
    pattern=('fibrous'if role in ('serosa','subserosa','adventitia')else'circular'if'circular'in role else
             'longitudinal'if'longitudinal'in role else'cellular'if smooth_mucosa else
             'mucosal'if role=='mucosa'else'oblique'if'muscle'in role or role in ('trigone-deep','trigone-superficial')else'solid')
    return dict(color=color,pattern=pattern,scale=.00065 if'mucosa'in role else .00045,relief=.025 if smooth_mucosa else .08)


def propose(preparation_dir,job_result,output,*,inventory=None):
    preparation_dir,job_result,output=map(lambda p:Path(p).resolve(),(preparation_dir,job_result,output))
    if output.exists():raise FileExistsError(output)
    inventory=Path(inventory).resolve()if inventory else ROOT/'.local-models/geometry-workbench/bladder-completion/source-v1/male-inventory.json'
    prep=json.loads((preparation_dir/'preparation.json').read_text());d=json.loads((Path(prep['delivery_dir'])/'delivery.json').read_text())
    job=json.loads(job_result.read_text());result=job.get('result',job)
    if not result.get('passed')or len(result['assemblies'])!=len(prep['groups']):raise ValueError('Use the complete passed packing job')
    current=Path(prep['catalogue']);before=json.loads(current.read_text())
    if digest_file(current)!=prep['catalogue_sha256']:raise ValueError('Source catalogue changed')
    a=copy.deepcopy(before);byid={p['id']:p for p in a['parts']};concepts={c['id']:c for c in a['concepts']};reader=AtlasReader()
    with sqlite3.connect(ROOT/'.local-models/ta98-authority.sqlite')as c:
        terms={r[0]:(r[1],r[2])for r in c.execute('SELECT id,name_en,parent_id FROM ta98')}
    group_rows={};new_ids=[];charts=[];protected=0;actual_meshes={}
    region_template=byid['TA98PART:ZA:Urinary bladder:A08.3.01.003'];urethra_template=byid['TA98PART:ZA:Urethra:A09.4.02.004']
    for group,r in zip(prep['groups'],result['assemblies']):
        bank=Path(r['bank'])
        if not r['passed']or digest_file(bank)!=r['bank_sha256']:raise ValueError('A packed job bank changed')
        zipped=bank.with_suffix(bank.suffix+'.gz');atomic_bytes(zipped,gzip.compress(bank.read_bytes(),mtime=0))
        url='/local-models/'+str(bank.relative_to(ROOT/'.local-models'));chunk=len(a['chunks'])
        a['chunks'].append(dict(url=url,gzip=url+'.gz',bytes=bank.stat().st_size,gzipBytes=zipped.stat().st_size,sha256=r['bank_sha256']))
        group_rows[group]=r
        if not group.startswith('node:'):charts.append(dict(name=group,partIds=list(r['parts']),bankSha256=r['bank_sha256']))
        for pid,packed in r['parts'].items():
            if group.startswith('node:'):
                expected=Mesh.load(next(q['path']for q in json.loads((Path(prep['node_dir'])/'candidates.json').read_text())['rows']if q['part_id']==pid))
            else:expected=Mesh.load(Path(prep['delivery_dir'])/group/(pid.replace(':','_')+'.npz'))
            actual=Mesh.load(bank.parent/(pid.replace(':','_')+'.npz'))
            actual_meshes[pid]=actual
            if not np.array_equal(actual.triangles,expected.triangles):raise ValueError('Job packing changed an actual oriented tissue triangle: '+pid)
            for key in ('normals','uvs'):
                if key not in expected.attributes:continue
                x=expected.attributes[key][expected.faces];y=actual.attributes[key][actual.faces]
                if key=='normals':
                    x=x.astype(np.int32)*expected.attributes.get('normalSigns',np.ones(len(expected.vertices),int))[expected.faces,None]
                    y=y.astype(np.int32)*actual.attributes.get('normalSigns',np.ones(len(actual.vertices),int))[actual.faces,None]
                if not np.array_equal(x,y):raise ValueError('Job packing changed a literal source corner attribute: '+pid)
            if pid.startswith('RENALDETAIL:')and':ureter-wall:'not in pid:protected+=1
            template=byid.get(pid,region_template if group=='bladder'else urethra_template);p=copy.deepcopy(template)
            for key in GEOMETRY_KEYS:p.pop(key,None)
            for key in ('suppressed','defaultHidden','replacementParts'):p.pop(key,None)
            p.update(copy.deepcopy(packed),chunk=chunk)
            if not group.startswith('node:'):
                for draw,use in zip(p['tissueDrawGroups'],d['cells'][pid]):
                    if use.get('neighbor'):draw['neighbor']=use['neighbor']
            if group in ('bladder','urethra'):
                region,role=pid.split(':',2)[2].split(':');p.update(name=LABELS[role]+' of '+region+' '+('urinary bladder'if group=='bladder'else'urethra'),
                    conceptId=pid,sourceId='shared-urinary:'+pid,urinaryOrgan=group,urinaryRegion=region,
                    urinarySourcePartId='ZA:Urinary bladder'if group=='bladder'else'ZA:Urethra',system='urinary',
                    pairedTissue=True,sectionSolid=True,tissueWallAssembly='urinary:male:'+group,tissueWallLayer=role,coverageLimitation=LIMITS,
                    reviewStatus='local-finite-urinary-wall; anatomical review ongoing')
                if role.startswith('mucosa-ureteric-orifice-'):p['name']+=' ('+role.rsplit('-',1)[-1]+')'
                term='A08.3.01.'+BLADDER_TERMS[role]if group=='bladder'and role in BLADDER_TERMS else None
                # The material cell is mucosa. The opening itself is the
                # separately witnessed empty boundary in its term group.
                if group=='bladder'and 'orifice'in role:term='A08.3.01.023'
                if group=='urethra':term=({'prostatic':{'muscularis':'A09.4.02.011','mucosa':'A09.4.02.015'},'intermediate':{'muscularis':'A09.4.02.018','mucosa':'A09.4.02.020'},'spongy':{'muscularis':'A09.4.02.027','mucosa':'A09.4.02.029'}}[region].get(role))
                p.pop('ta98Term',None)
                if term:p['ta98Term']=term
                ap=appearance(role);material='URINARY:male:'+group+':'+role;p['material']=material;p['tissueAppearance']=ap
                a['materials'][material]=dict(color=[round(x*255)for x in ap['color']],roughness=.78,metalness=0,pairedNormalSigns=True,tissueAppearance=ap,geometryOrigin=p['geometryOrigin'])
                c=dict(id=pid,name=p['name'],elements=[pid],partOf='ZA:Urinary bladder'if group=='bladder'else'ZA:Urethra',coverageLimitation=LIMITS)
                if term:c.update(ta98Term=term,ta98Kind='exact')
                concepts[pid]=c
            elif group=='prostate':
                material='URINARY:male:prostate';a['materials'][material]=dict(a['materials'][p['material']],pairedNormalSigns=True)
                p.update(material=material,pairedTissue=True,sectionSolid=True,tissueWallAssembly='urinary:male:prostate',tissueWallLayer='parenchyma',
                    coverageLimitation='Source exterior and UVs retained; a shared urinary bore replaces the blocked source outlet. Prostate subdivisions and duct openings are unfinished.')
                p['tissueAppearance']=dict(color=[.72,.58,.4],pattern='cellular',scale=.0005,relief=.04)
            if pid in byid:byid[pid].clear();byid[pid].update(p)
            else:a['parts'].append(p);byid[pid]=p
            new_ids.append(pid)
    if protected!=60:raise ValueError('Expected all sixty retained kidney tissue fields')
    replacements={}
    for region,code in [('body','003'),('fundus','004'),('neck','005')]:
        replacements['TA98PART:ZA:Urinary bladder:A08.3.01.'+code]=[p for p in group_rows['bladder']['parts']if p.split(':')[2]==region]
    for region,code in [('prostatic','004'),('intermediate','017'),('spongy','021')]:
        replacements['TA98PART:ZA:Urethra:A09.4.02.'+code]=[p for p in group_rows['urethra']['parts']if p.split(':')[2]==region]
    replacements['ZA:Urinary bladder']=list(group_rows['bladder']['parts']);replacements['ZA:Urethra']=list(group_rows['urethra']['parts'])
    for pid,members in replacements.items():byid[pid].update(suppressed=True,defaultHidden=True,replacementParts=members)
    for c in concepts.values():
        c['elements']=list(dict.fromkeys(p for old in c['elements']for p in replacements.get(old,[old])))
        if'componentElements'in c:c['componentElements']=list(dict.fromkeys(p for old in c['componentElements']for p in replacements.get(old,[old])))
    for source,members in replacements.items():
        if source in concepts:
            concepts[source]['elements']=members
            if 'componentElements'in concepts[source]:concepts[source]['componentElements']=members
    term_groups=[];groups=bladder_term_groups(group_rows['bladder']['parts'])
    for term,members in groups.items():
        cid='URINARY:male:bladder:'+term;concepts[cid]=dict(id=cid,name=terms[term][0].capitalize(),elements=members,componentElements=members,
            partOf='ZA:Urinary bladder',ta98Term=term,ta98Kind='exact',coverageLimitation=LIMITS)
        term_groups.append(dict(term=term,conceptId=cid,partIds=members))
    if d.get('opening_curves'):
        from .urinary_openings import bind_opening_boundary
        for term,names in [('A08.3.01.026',('left','right')),('A08.3.01.027M',('neck',))]:
            group=next(g for g in term_groups if g['term']==term);boundaries=[]
            for name in names:
                role='mucosa-internal-urethral-orifice' if name=='neck' else 'mucosa-ureteric-orifice-'+name
                support={pid:(actual_meshes[pid],byid[pid]) for pid in group['partIds']if byid[pid]['tissueWallLayer']==role}
                boundaries.append(bind_opening_boundary(name,d['opening_curves'][name],support,d['feature_plan_sha256']))
            group['openingBoundaries']=boundaries
    for region,term in [('prostatic','A09.4.02.004'),('intermediate','A09.4.02.017'),('spongy','A09.4.02.021')]:
        members=[pid for pid in group_rows['urethra']['parts']if pid.split(':')[2]==region];cid='URINARY:male:urethra:'+region
        concepts[cid]=dict(id=cid,name=terms[term][0].capitalize(),elements=members,componentElements=members,partOf='ZA:Urethra',ta98Term=term,ta98Kind='exact',coverageLimitation=LIMITS)
        term_groups.append(dict(term=term,conceptId=cid,partIds=members))
    a['concepts']=list(concepts.values())
    bindings=a.setdefault('ta98TermBindings',[])
    generated={pid for pid in new_ids if byid[pid].get('urinaryOrgan')}
    bindings[:]=[b for b in bindings if b.get('conceptId')not in generated]
    for binding in bindings:
        if binding.get('conceptId')in replacements:binding['partIds']=replacements[binding['conceptId']]
    for pid in new_ids:
        p=byid[pid];c=concepts.get(p['conceptId'])
        if p.get('urinaryOrgan')and p.get('ta98Term')and c and c.get('ta98Kind')=='exact':
            bindings.append(dict(term=p['ta98Term'],conceptId=c['id'],partIds=[pid],kind='exact'))
    for g in a['finiteRenalIntegration']['assemblies']:g['bankSha256']=group_rows[g['side']]['bank_sha256']
    a['finiteUrinaryIntegration']=dict(schema='urinary-local-coats-v3'if d.get('opening_curves')else'urinary-local-coats-v2'if 'A08.3.01.025'in groups else'urinary-local-coats-v1',sex='male',partIds=list(d['part_ids']),charts=charts,
        termGroups=term_groups,externalInterfaces=d['external_interfaces'],sourceInventory=str(inventory),
        sourceDeliverySha256=digest_file(Path(prep['delivery_dir'])/'delivery.json'),wholeGeometryProof=d['whole_geometry_verification'],
        openingPlanSha256=d.get('feature_plan_sha256')if d.get('opening_curves')else None,
        retainedKidneyTissueCells=60,renalRootTrianglesPreserved=620,ordinaryLocalAtlas=True,anatomicallyAccepted=False)
    a['defaultHidden']=sorted((set(a.get('defaultHidden',[]))|set(replacements))-set(new_ids))
    a['partsCount']=len(a['parts']);a['triangles']=sum(p['indexCount']//3 for p in a['parts'])
    previous={p['id']:p for p in before['parts']};changed=[p['id']for p in a['parts']if p!=previous.get(p['id'])]
    output.mkdir(parents=True);atomic_json(output/'atlas.json',a)
    atomic_json(output/'proposal.json',dict(model='local-ta98',proposed_catalogue=str(output/'atlas.json'),expected_catalogue_sha256=digest_file(current),
        changed_parts=changed,inventory=str(inventory),
        assembly=dict(whole_organ='Male urinary bladder with continuous ureters and outlet',
            adjoining_interfaces='Literal renal-root annuli, bladder ostia, urethral wall/prostate bore; node packing replay; skeleton fixed'),
        require_closed=True,max_pairs=1000000,job_id=job.get('job_id')))
    return dict(proposed_catalogue=str(output/'atlas.json'),changed_parts=len(changed),active_urinary_cells=len(d['part_ids']),protected_kidney_cells=protected)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare-job');s.add_argument('--delivery-dir',required=True);s.add_argument('--node-dir');s.add_argument('--output',required=True)
    s=sub.add_parser('propose');s.add_argument('--preparation-dir',required=True);s.add_argument('--job-result',required=True);s.add_argument('--output',required=True);s.add_argument('--inventory')
    args=vars(p.parse_args());command=args.pop('command');print(json.dumps(prepare_job(**args)if command=='prepare-job'else propose(**args),indent=2))
