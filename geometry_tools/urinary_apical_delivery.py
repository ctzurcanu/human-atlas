"""Prepare the source-bound apical attachment revision; never install it here."""
from pathlib import Path
import copy, gzip, json
import numpy as np
from .common import ROOT, atomic_bytes, atomic_json, digest_file, identity
from .model import AtlasReader, Mesh
from .constraints import check_constraint
from .reference_refit import verify_reference_refit


def prepare(output, *, bank_directory=None):
    output=Path(output).resolve()
    if output.exists(): raise FileExistsError(output)
    b=ROOT/'.local-models/geometry-workbench/bladder-completion'
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    sources=json.loads((b/'male-apical-sources-v1/sources-and-assembly.json').read_text())
    if digest_file(source)!=sources['source_catalogue_sha256']: raise ValueError('Apical source catalogue changed')
    proofs=[b/p for p in ['male-apical-cord-v3/whole-assembly-verification.json',
        'male-apical-cord-v3/packing-fixed-neighbors.json','male-omentum-delivery-v1/body-packing.json']]
    if not all(json.loads(p.read_text())['passed']for p in proofs): raise ValueError('Current packing evidence did not pass')
    a=json.loads(source.read_text());before=copy.deepcopy(a);reader=AtlasReader()
    by_id={p['id']:p for p in a['parts']};marker=a['finiteUrinaryIntegration'];frame=sources['frame']
    wall_id='BLADDER:male:body:adventitia';lig_id='BLADDER:male:median-umbilical-ligament';fascia_id='ZA:Transversalis fascia'
    wall=by_id[wall_id];wall_chart=next(g for g in marker['charts']if g['name']=='bladder')
    plan_path=b/'male-apical-sources-v1/plan-corridor-v3.json'
    # The construction's literal first-contact faces, not an independently
    # selected near-contact patch, determine the attached adventitial facets.
    cord=Mesh.load(b/'male-apical-cord-v3/packed/BLADDER_male_median-umbilical-ligament.npz')
    old_wall=reader.read_part(a,wall,frame=frame)
    cap={tuple(sorted(map(tuple,t)))for t in cord.triangles[:32]}
    contacts=np.array([i for i,t in enumerate(old_wall.triangles)if tuple(sorted(map(tuple,t)))in cap],int)
    if len(contacts)!=32: raise ValueError('Actual packed apical cap no longer matches its source wall')
    old_bank=reader.asset(a['chunks'][wall['chunk']]['url']);raw=bytearray(old_bank.read_bytes())
    indices=np.frombuffer(raw,dtype='<u4',count=wall['indexCount'],offset=wall['indices']).reshape(-1,3).copy()
    free=wall['tissueDrawGroups'][-1]
    if free.get('neighbor') or np.any(contacts<free['start']//3): raise ValueError('The selected attachment is not exposed adventitia')
    kept=np.array([i for i in range(len(indices))if i not in set(contacts)],int)
    reordered=indices[np.r_[kept,contacts]].astype('<u4')
    raw.extend(b'\0'*((-len(raw))%4));wall['indices']=len(raw);raw.extend(reordered.tobytes())
    free['count']-=len(contacts)*3
    wall['tissueDrawGroups'].append(dict(start=wall['indexCount']-96,count=96,neighbor=lig_id))
    output.mkdir(parents=True)
    bank_directory=Path(bank_directory).resolve()if bank_directory else output
    def bank(name,data):
        path=bank_directory/(name+'.bin')
        if path.exists():
            if path.read_bytes()!=data:raise ValueError('An existing immutable output bank differs')
        else:atomic_bytes(path,data)
        compressed=path.with_suffix('.bin.gz')
        if compressed.exists():
            if gzip.decompress(compressed.read_bytes())!=data:raise ValueError('An existing compressed bank differs')
        else:atomic_bytes(compressed,gzip.compress(data,mtime=0))
        chunk=dict(url='/local-models/'+str(path.relative_to(ROOT/'.local-models')),bytes=len(data),sha256=digest_file(path),
            gzip='/local-models/'+str(compressed.relative_to(ROOT/'.local-models')),gzipBytes=compressed.stat().st_size)
        a['chunks'].append(chunk)
        return len(a['chunks'])-1,chunk['sha256']
    chunk,wall_sha=bank('bladder',bytes(raw))
    for pid in wall_chart['partIds']: by_id[pid]['chunk']=chunk
    wall_chart['bankSha256']=wall_sha
    packed=json.loads((b/'male-apical-cord-v3/packed/tissue.json').read_text())
    lig=copy.deepcopy(packed['parts'][lig_id]);chunk,lig_sha=bank('median-umbilical-ligament',Path(packed['bank']).read_bytes())
    lig.update(chunk=chunk,name='Median umbilical ligament',conceptId=lig_id,sourceId='constructed:median-umbilical-ligament',system='urinary',
        pairedTissue=True,sectionSolid=True,tissueWallAssembly='urinary:male:bladder',tissueWallLayer='median-umbilical-ligament',
        urinaryOrgan='bladder',urinaryRegion='attachment',urinarySourcePartId='ZA:Urinary bladder',ta98Term='A08.3.01.007',ta98Kind='exact',
        groups=['8: Visceral systems','Abdomen','Pelvis','Trunk','Urinary system'],material='URINARY:male:median-umbilical-ligament',
        license='Local anatomical construction; source contact geometry retains its existing license',
        reviewStatus='local-finite-attachment; source fascia is an open reference',
        description='Fibromuscular remnant joining the bladder apex to the anterior abdominal wall beneath the umbilicus.',
        tissueAppearance=dict(color=[.75,.65,.49],pattern='fibrous',scale=.0004,relief=.08))
    lig['tissueDrawGroups'][0]['neighbor']=wall_id
    lig['tissueDrawGroups'][-1]['neighbor']=fascia_id
    a['parts'].append(lig);by_id[lig_id]=lig
    a['materials'][lig['material']]=dict(color=[.75,.65,.49],pairedNormalSigns=True,tissueAppearance=lig['tissueAppearance'])
    marker['partIds'].append(lig_id);marker['charts'].append(dict(name='attachments',partIds=[lig_id],bankSha256=lig_sha))
    a['concepts'].append(dict(id=lig_id,name=lig['name'],elements=[lig_id],ta98Term=lig['ta98Term'],ta98Kind='exact',partOf=['ZA:Urinary bladder']))
    marker['termGroups'].append(dict(term=lig['ta98Term'],conceptId=lig_id,partIds=[lig_id]))
    for c in a['concepts']:
        if c['id']=='ZA:Urinary bladder': c['elements'].append(lig_id)
    by_id['ZA:Urinary bladder']['replacementParts'].append(lig_id)
    for group in marker['termGroups']:
        if group['term']=='A08.3.01.001':
            group['partIds'].append(lig_id)
            next(c for c in a['concepts']if c['id']==group['conceptId'])['elements']=list(group['partIds'])
    # Independently replay the encoded, now contiguous, reciprocal interface.
    wm=reader.read_part(a,wall,frame=frame);lm=reader.read_part(a,lig,frame=frame)
    wall_faces=list(range(wall['indexCount']//3-32,wall['indexCount']//3));lig_faces=list(range(32))
    interface=dict(kind='shared-interface',first=wall_id,second=lig_id,
        parameters=dict(first_faces=wall_faces,second_faces=lig_faces),evidence=dict(source=str(plan_path)))
    interface_check=check_constraint(interface,wm,lm)
    if not interface_check['passed']: raise ValueError('Emitted adventitial/apical seam failed independent replay')
    def field(mesh,faces):
        words=mesh.attributes['normals'].astype(int)*mesh.attributes.get('normalSigns',np.ones(len(mesh.vertices),int))[:,None]
        return {tuple(sorted(map(tuple,t))):{tuple(q):tuple(map(int,w))for q,w in zip(t,n)}
            for t,n in zip(mesh.triangles[faces],words[mesh.faces[faces]])}
    seam=field(wm,wall_faces)
    marker['externalInterfaces'].append(dict(first=wall_id,second=lig_id,triangles=32,
        coordinate_normal_signature=identity([(k,sorted(v.items()))for k,v in sorted(seam.items())]),
        actual_decoded_positions_and_opposed_words_match=True))
    fascia=reader.read_part(a,by_id[fascia_id],frame=frame);end=lig['indexCount']//3-1
    ref_a,ref_b=field(lm,[end]),field(fascia,[1796])
    if ref_a!=ref_b: raise ValueError('Terminal contact is not literal source fascia geometry and normal field')
    marker['apicalAttachment']=dict(partId=lig_id,wallPartId=wall_id,sourcePlanSha256=digest_file(plan_path),
        referenceContact=dict(first=lig_id,second=fascia_id,firstFace=end,secondFace=1796,triangles=1,
            representation='open-reference-sheet-contact',coordinate_normal_signature=identity([(k,sorted(v.items()))for k,v in sorted(ref_a.items())]),
            actual_decoded_positions_and_same_words_match=True))
    delivery=json.loads((b/'male-omentum-delivery-v1/delivery.json').read_text());new_omentum=copy.deepcopy(delivery['part'])
    chunk,_=bank('greater-omentum',Path(delivery['bank']).read_bytes());new_omentum['chunk']=chunk
    for i,p in enumerate(a['parts']):
        if p['id']==new_omentum['id']: a['parts'][i]=new_omentum;break
    refit={key:str(b/path)for key,path in dict(refined_source='male-omentum-refined-source-v1/source.npz',
        lineage='male-omentum-refined-source-v1/source-lineage.npz',controls='male-omentum-corridor-v6/vertex-controls.npz').items()}
    refit.update(schema='open-source-z-refit/v1',source_geometry_sha256=Mesh.load(b/'male-omentum-corridor-v1/source.npz').geometry_sha,
        encoding_tolerance_m=1e-8,maximum_z_shift_m=.0026)
    for key in ['refined_source','lineage','controls']:refit[key+'_sha256']=digest_file(refit[key])
    verify_reference_refit(reader.read_part(before,by_id['ZA:Greater omentum'],frame=frame),
        reader.read_part(a,new_omentum,frame=frame),refit)
    changed=[p['id']for p in a['parts']if p!=next((q for q in before['parts']if q['id']==p['id']),None)]
    a['triangles']=sum(p.get('triangles',p['indexCount']//3)for p in a['parts'])
    catalogue=output/'atlas.json';atomic_json(catalogue,a)
    spec=dict(kind='catalogue-candidate',parameters=dict(model='local-ta98',proposed_catalogue=str(catalogue),expected_catalogue_sha256=digest_file(source),
        changed_parts=changed,inventory=str(b/'source-v4/male-inventory.json'),open_reference_refits={'ZA:Greater omentum':refit},
        assembly=dict(whole_organ='Male bladder, median umbilical ligament and source-bound local omental clearance',
            adjoining_interfaces='Literal adventitial attachment; source fascia contact; all original urinary lumens, renal roots, skeleton and bowel fixed'),
        interface_checks=[interface],evidence=dict(packing=[str(p)for p in proofs],plan=str(plan_path),
            omentum_inventory=str(b/'male-omentum-corridor-v1/inventory.json'))))
    atomic_json(output/'job.json',spec);atomic_json(output/'interface-check.json',interface_check)
    return dict(catalogue=str(catalogue),job_spec=str(output/'job.json'),changed_parts=changed,
        attachment_triangles=lig['indexCount']//3,new_raw_bytes=sum(a['chunks'][i]['bytes']for i in range(len(before['chunks']),len(a['chunks']))))
