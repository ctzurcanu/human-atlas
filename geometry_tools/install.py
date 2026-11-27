"""Local-only compare-and-swap catalogue installation with guarded rollback.

Preparation validates a complete proposed catalogue. Publication reuses the
existing atlas revision adapter so visits, TA98 and ownership forest update
together, with all overwritten files backed up. No Git or deployment calls.
"""
from pathlib import Path
import copy,gzip,json,os,importlib.util,time,fcntl
from .common import ROOT,CACHE,atomic_json,atomic_bytes,digest_file,identity
from .model import AtlasReader
from .audit import audit
from .jobs import row as job_row,artifacts,verify_bindings,implementation_signature


def replay_tissue_cell_audit(mesh,directory,*,max_pairs,max_tests=2000000,audited_layouts=(),canonical_frame=None):
    """Replay a saved full audit only for the literal emitted geometry/words.

    Metadata and names cannot authorize reuse. The numerical implementation,
    native version, frame, geometry, attributes and result checksum must match.
    A missing or superseded checkpoint falls back to the ordinary audit.
    """
    if not directory:return None
    if canonical_frame is not None:
        from .geometry_instance import instance_preimage
        from .model import Mesh
        placed=instance_preimage(mesh)
        if placed is None:return None
        local,pose=placed
        if local.attributes.keys()!=mesh.attributes.keys():raise ValueError('Instance attribute set does not replay its literal source bank')
        if not isinstance(canonical_frame,str)or not canonical_frame.strip():raise ValueError('Explicit canonical audit frame required')
        canonical=Mesh(local.vertices,local.faces,canonical_frame,local.attributes,local.source)
        replay=replay_tissue_cell_audit(canonical,directory,max_pairs=max_pairs,max_tests=max_tests,audited_layouts=audited_layouts)
        if replay is None:return None
        result=copy.deepcopy(replay);result.update(geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,frame=mesh.frame,geometry_instance=pose,canonical_audit_frame=canonical_frame)
        result['topology']['surface_area_m2']*=pose['scale']**2
        result['topology']['oriented_triangle_integral_m3']*=pose['scale']**3
        for key in ['first_candidate_patch_area_upper_bound_m2','second_candidate_patch_area_upper_bound_m2']:
            if key in result.get('self_intersections',{}):result['self_intersections'][key]*=pose['scale']**2
        return result
    import numpy as np
    import importlib.metadata
    root=Path(__file__).parent
    signature=identity('finite-tissue-validation-checkpoint/v4',np.__version__,importlib.metadata.version('manifold3d'),importlib.metadata.version('libigl'),
        [(name,digest_file(root/name)) for name in ['tissue.py','audit.py','contact.py','distance.py','spatial.py','model.py','geometry_instance.py']],
        digest_file(root/'algorithms/finite_tissue_boolean.py'),digest_file(root/'algorithms/shared_cell_disjoint.py'),max_pairs,max_tests)
    geometry=[mesh.geometry_sha,mesh.frame];key=identity(signature,'cell-audit',geometry);path=Path(directory).resolve()/(key+'.json')
    if not path.exists():return None
    record=json.loads(path.read_text());result=record.get('result',{})
    if record.get('signature')!=signature or record.get('geometry')!=geometry or record.get('kind')!='cell-audit' or record.get('result_sha256')!=identity(result):
        raise ValueError('Tissue audit checkpoint is inconsistent or edited')
    if result.get('geometry_sha256')!=mesh.geometry_sha or result.get('frame')!=mesh.frame:return None
    equivalence=None
    if result.get('attributes_sha256')!=mesh.attributes_sha:
        from .model import Mesh
        if len(np.unique(mesh.faces))!=len(mesh.vertices):return None
        # Packing retains unused controls; AtlasReader drops them. Reuse is
        # valid only if the *audited* original buffer is reproduced literally
        # and every emitted triangle and attribute corner is exactly equal.
        pid=mesh.source.get('record',{}).get('id')
        for layout in audited_layouts:
            row=layout['parts'].get(pid)
            if row is None or layout['frame']!=mesh.frame:continue
            bank=layout['_bank'];n=row['vertexCount'];reader=AtlasReader()
            vertices=reader.view(bank,row['positions'],n,'<f4',3).astype(float)+np.asarray(layout['origin'])
            faces=reader.view(bank,row['indices'],row['indexCount'],'<u4').reshape(-1,3)
            attributes={k:reader.view(bank,row[k],n,dtype,width)for k,dtype,width in
                        [('normals','<i2',3),('normalSigns','i1',1),('uvs','<f4',2),('colors','u1',3)]if k in row}
            original=Mesh(vertices,faces,mesh.frame,attributes)
            if original.attributes_sha!=result.get('attributes_sha256'):continue
            if (original.geometry_sha!=mesh.geometry_sha or not np.array_equal(original.triangles,mesh.triangles)
                or original.attributes.keys()!=mesh.attributes.keys()
                or any(not np.array_equal(v[original.faces],mesh.attributes[k][mesh.faces])for k,v in original.attributes.items())):
                raise ValueError('Audited buffer does not reproduce the literal emitted triangle/corner fields')
            equivalence=dict(audited_layout_path=layout['_path'],audited_layout_sha256=layout['_sha256'],
                audited_bank_path=layout['bank'],audited_bank_sha256=layout['bank_sha256'],
                audited_attributes_sha256=original.attributes_sha,all_emitted_triangles_and_corner_fields_literal=True,
                only_unreferenced_controls_removed=True)
            break
        if equivalence is None:return None
    if result.get('passed') and (not result.get('topology',{}).get('closed_oriented_manifold')
        or not result.get('self_intersections',{}).get('complete_within_pair_budget')
        or result['self_intersections'].get('distance_threshold_m')!=0
        or result['self_intersections'].get('defect_pairs')!=0 or result.get('failures')):
        raise ValueError('Tissue checkpoint does not contain the complete required closed-cell audit')
    return dict(result,attributes_sha256=mesh.attributes_sha,replayed_from_checkpoint=str(path),
                replayed_checkpoint_sha256=digest_file(path),buffer_equivalence=equivalence)


def catalogue_diff(before,after,allowed):
    bp={p['id']:p for p in before['parts']};ap={p['id']:p for p in after['parts']}
    if len(bp)!=len(before['parts'])or len(ap)!=len(after['parts']):raise ValueError('Duplicate catalogue IDs')
    changed={key for key in bp.keys()|ap.keys()if bp.get(key)!=ap.get(key)}
    if changed!=set(allowed):raise ValueError('Actual changed item set does not equal the explicit replacement set')
    if after['chunks'][:len(before['chunks'])]!=before['chunks']:raise ValueError('Existing banks must remain literal; append replacement banks instead')
    for pid in changed:
        if bp.get(pid,{}).get('system')=='skeletal' or ap.get(pid,{}).get('system')=='skeletal':raise ValueError('This local organ installer cannot replace the fixed skeleton')
    return bp,ap,sorted(changed)


def changed_constraint_meshes(registrations,meshes):
    """Use the proposed corner fields even when positions/faces stay literal.

    Shared-interface checks compare normal words as well as coordinates.
    Reloading the historical endpoint for an attribute-only edit would compare
    two different revisions and can both reject a valid seam and miss a bad one.
    """
    return {record['id']:meshes[pid] for pid,record in registrations.items()
            if record['geometry_sha256']!=meshes[pid].geometry_sha
            or record.get('attributes_sha256')!=meshes[pid].attributes_sha}


def transaction(changes,backup_dir,*,expected,after_write=None):
    """Generic journaled file transaction, also used by fixture tests.

    expected[path]=None means that target must be new. A failed callback
    restores literal previous bytes. Recovery rejects later user edits.
    """
    backup_dir=Path(backup_dir).resolve()
    if backup_dir.exists():raise FileExistsError(backup_dir)
    paths={Path(p).resolve():data for p,data in changes.items()};expected={str(Path(p).resolve()):sha for p,sha in expected.items()}
    for p in paths:
        if p.name=='.gitignore'or '.git'in p.parts:raise ValueError('Git configuration and metadata are outside geometry installation scope')
        sha=expected[str(p)]
        if (digest_file(p)if p.exists()else None)!=sha:raise ValueError('Compare-and-swap source mismatch: '+str(p))
    backup_dir.mkdir(parents=True);entries=[]
    for i,(p,data)in enumerate(paths.items()):
        backup=backup_dir/(str(i)+'.gz')
        if p.exists():atomic_bytes(backup,gzip.compress(p.read_bytes(),mtime=0))
        entries.append(dict(path=str(p),previous_sha256=expected[str(p)],new_sha256=identity_bytes(data),backup=str(backup)if p.exists()else None))
    manifest=backup_dir/'transaction.json';atomic_json(manifest,dict(state='prepared',files=entries))
    try:
        for p,data in paths.items():atomic_bytes(p,data,overwrite=True)
        if after_write:after_write()
        atomic_json(manifest,dict(state='installed',files=entries),overwrite=True)
    except BaseException:
        restore_transaction(manifest);raise
    return dict(manifest=str(manifest),state='installed')


def identity_bytes(data):
    import hashlib
    return hashlib.sha256(data).hexdigest()


def restore_transaction(manifest_path):
    manifest_path=Path(manifest_path);m=json.loads(manifest_path.read_text())
    for entry in m['files']:
        p=Path(entry['path']);current=digest_file(p)if p.exists()else None
        if current not in [entry['previous_sha256'],entry['new_sha256']]:raise ValueError('Rollback would overwrite a later edit: '+str(p))
    for entry in reversed(m['files']):
        p=Path(entry['path'])
        if entry['backup']:
            raw=gzip.decompress(Path(entry['backup']).read_bytes())
            if identity_bytes(raw)!=entry['previous_sha256']:raise ValueError('Corrupt rollback snapshot')
            atomic_bytes(p,raw,overwrite=True)
        elif p.exists():
            # New generated files are retained, moved into the recoverable
            # transaction folder rather than permanently deleted.
            os.replace(p,manifest_path.parent/('retained-'+identity(str(p))[:12]+p.suffix))
    atomic_json(manifest_path,dict(state='rolled-back',files=m['files']),overwrite=True)
    return dict(state='rolled-back',manifest=str(manifest_path))


def emitted_requires_closed(spec,part_id):
    """Optional open surface views cannot relax an explicitly finite wall."""
    return bool(spec.get('require_closed',True)or part_id in spec.get('closed_parts',()))


def prepare(db,spec,output):
    model=spec['model'];sex={'local-ta98':'male','local-female-ta98':'female'}.get(model)
    if not sex:raise ValueError('Only the two local TA98 catalogues may be prepared here')
    current=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json';proposed=Path(spec['proposed_catalogue']).resolve()
    if digest_file(current)!=spec['expected_catalogue_sha256']:raise ValueError('Current atlas changed; reprepare from its actual revision')
    before=json.loads(current.read_bytes());after=json.loads(proposed.read_bytes());bp,ap,changed=catalogue_diff(before,after,spec['changed_parts'])
    closed_parts=spec.get('closed_parts',[])
    if not isinstance(closed_parts,list)or len(set(closed_parts))!=len(closed_parts)or any(pid not in ap for pid in closed_parts):
        raise ValueError('Explicit finite parts must name distinct emitted catalogue parts')
    invpath=Path(spec['inventory']).resolve();inventory_bytes=invpath.read_bytes();inv=json.loads(inventory_bytes);inventory_sha=identity_bytes(inventory_bytes)
    from .inventory_validation import validate_authoritative_inventory
    inventory_validation=validate_authoritative_inventory(inv,model=model,current_catalogue=current)
    assembly=spec.get('assembly')
    if not isinstance(assembly,dict)or not assembly.get('adjoining_interfaces')or not assembly.get('whole_organ'):raise ValueError('Declare the whole assembly and every adjoining interface before installation')
    reader=AtlasReader();checks={};banks={};registrations={};mesh_by_part={};audited_layouts=[]
    # Thousands of finite cells may share one large bank and catalogue. Bind
    # each file once per unchanged stat identity; recheck all bindings before
    # returning the proposal. Geometry/corner and interface audits still run.
    proposed_sha=reader.file_sha(proposed)
    for path in spec.get('tissue_audit_layouts',[]):
        import numpy as np
        path=Path(path).resolve();layout=json.loads(path.read_text())
        verify_bindings({layout['bank']:layout['bank_sha256']})
        layout.update(_path=str(path),_sha256=digest_file(path),_bank=np.memmap(layout['bank'],dtype='u1',mode='r'))
        audited_layouts.append(layout)
    from .landmarks import validate_catalogue as validate_landmarks
    point_checks=validate_landmarks(after,reader=reader,frame='model:'+model+':assembled-metres')
    for reference in point_checks:banks[reference['bank_path']]=reference['bank_sha256']
    print(f'[install] validating {len(changed)} explicitly changed parts',flush=True)
    for ordinal,pid in enumerate(changed):
        if pid not in ap:continue
        m=reader.read_part(after,ap[pid],frame='model:'+model+':assembled-metres',catalogue_path=proposed,catalogue_sha=proposed_sha)
        if ap[pid].get('suppressed'):
            if pid not in bp:raise ValueError('New replacements must be emitted and reviewed, rather than installed hidden')
            original=reader.read_part(before,bp[pid],frame=m.frame,catalogue_path=current)
            if original.geometry_sha!=m.geometry_sha or original.attributes_sha!=m.attributes_sha:raise ValueError('A retained suppressed source must keep its literal geometry and attributes')
            checks[pid]=dict(passed=True,retained_source_geometry_literal=True,not_a_new_emitted_replacement=True)
            continue
        policy=spec.get('open_reference_refits',{}).get(pid)
        if policy:
            if pid in closed_parts:raise ValueError('An explicitly finite tissue cannot use an open-reference policy')
            if pid not in bp or any(ap[pid].get(k)or bp[pid].get(k)for k in ['pairedTissue','sectionSolid']):
                raise ValueError('Only an existing open reference, not a new finite tissue, may use a source refit policy')
            from .reference_refit import verify_reference_refit
            original=reader.read_part(before,bp[pid],frame=m.frame,catalogue_path=current)
            refit=verify_reference_refit(original,m,policy)
        replay=None
        require_closed=emitted_requires_closed(spec,pid)
        if not policy and require_closed:
            replay=replay_tissue_cell_audit(m,spec.get('tissue_audit_checkpoints'),max_pairs=spec.get('max_pairs',250000),max_tests=spec.get('max_tests',2000000),audited_layouts=audited_layouts,canonical_frame=spec.get('tissue_audit_canonical_frame'))
        checks[pid]=replay if replay is not None else audit(m,require_closed=False if policy else require_closed,max_pairs=spec.get('max_pairs',250000),max_tests=spec.get('max_tests',2000000))
        if replay is not None:banks[replay['replayed_from_checkpoint']]=replay['replayed_checkpoint_sha256']
        if replay is not None and replay.get('buffer_equivalence'):
            proof=replay['buffer_equivalence']
            banks[proof['audited_layout_path']]=proof['audited_layout_sha256'];banks[proof['audited_bank_path']]=proof['audited_bank_sha256']
        if policy:checks[pid]['reference_refit']=refit
        bank=reader.asset(after['chunks'][ap[pid]['chunk']]['url']);banks[str(bank)]=reader.file_sha(bank)
        mesh_by_part[pid]=m
        if pid in bp:
            registrations[pid]=db.mesh(pid,model)
        if (ordinal+1)%10==0:print(f'[install] emitted audits {ordinal+1}/{len(changed)}',flush=True)
    if not all(r['passed']for r in checks.values()):raise ValueError('An emitted replacement fails geometry checks; inspect its independent audit before installing')
    changed_geometry=changed_constraint_meshes(registrations,mesh_by_part)
    from .constraints import check,check_constraint,list_constraints
    overrides={};rebinding_checks=[];endpoint_overrides={};transfer_checks=[]
    requested=spec.get('shared_interface_rebindings',[])
    if requested:
        existing={r['id']:r for r in list_constraints(db,model,list(changed_geometry))}
        for binding in requested:
            old=existing.get(binding['id'])
            if not old or old['kind']!='shared-interface'or old['parameters'].get('verification_pending')or binding['id']in overrides:
                raise ValueError('Only an explicit existing verified shared interface may be rebound')
            if identity(old['parameters'])!=binding['expected_parameters_sha256']or not binding.get('evidence',{}).get('source'):
                raise ValueError('Interface rebinding must bind its exact old parameters and construction evidence')
            first=db.mesh(old['first_mesh_id']);second=db.mesh(old['second_mesh_id'])
            if [first['part_id'],second['part_id']]!=[binding['first'],binding['second']]:
                raise ValueError('A shared interface rebinding cannot change its tissue endpoints')
            source_check=check_constraint(old,reader.read_part(before,bp[first['part_id']],frame='model:'+model+':assembled-metres'),
                reader.read_part(before,bp[second['part_id']],frame='model:'+model+':assembled-metres'))
            if not source_check['passed']:raise ValueError('An unverified old seam cannot authorize a replacement seam')
            for role,peer in [('first','second'),('second','first')]:
                record=ap[binding[role]]
                expected_faces=[i for g in record.get('tissueDrawGroups',[])if g.get('neighbor')==binding[peer]
                    for i in range(g['start']//3,(g['start']+g['count'])//3)]
                if not expected_faces or sorted(binding['parameters'][role+'_faces'])!=sorted(expected_faces):
                    raise ValueError('A replacement seam must cover the complete declared reciprocal draw groups')
            overrides[binding['id']]=binding['parameters']
            rebinding_checks.append(dict(id=binding['id'],source_check=source_check,expected_parameters_sha256=binding['expected_parameters_sha256'],
                parameters=binding['parameters'],evidence=binding['evidence'],database_original_unchanged=True))
    transfers=spec.get('shared_interface_transfers',[])
    if transfers:
        from .interface_transfer import verify_transfer
        existing={r['id']:r for r in list_constraints(db,model,list(changed_geometry))}
        for binding in transfers:
            old=existing.get(binding['id'])
            if not old or old['kind']!='shared-interface' or binding['id']in overrides:
                raise ValueError('Transfer requires a distinct affected historical shared interface')
            if identity(old['parameters'])!=binding['expected_parameters_sha256'] or not binding.get('evidence',{}).get('source'):
                raise ValueError('Transfer must bind historical parameters and explicit construction evidence')
            old_ids=[db.mesh(old[k])['part_id']for k in ('first_mesh_id','second_mesh_id')]
            if old_ids!=[binding['source_first'],binding['source_second']] or not any(ap[p].get('suppressed')for p in old_ids):
                raise ValueError('Transfer must name the literal old owners and a retained suppressed source')
            candidate=[];complete_parameters={}
            for role,peer in [('first','second'),('second','first')]:
                pid=binding[role];record=ap[pid]
                if record.get('suppressed')or not record.get('pairedTissue')or not record.get('sectionSolid'):
                    raise ValueError('An interface must transfer to actual active finite tissue')
                faces=[i for g in record.get('tissueDrawGroups',[])if g.get('neighbor')==binding[peer]
                       for i in range(g['start']//3,(g['start']+g['count'])//3)]
                if not faces or not binding['parameters'][role+'_faces'] or not set(binding['parameters'][role+'_faces'])<=set(faces):
                    raise ValueError('Transferred source facets must belong to the declared reciprocal draw group')
                complete_parameters[role+'_faces']=faces
                candidate.append(mesh_by_part[pid]if pid in mesh_by_part else reader.read_part(after,record,frame='model:'+model+':assembled-metres'))
            original=[reader.read_part(before,bp[pid],frame='model:'+model+':assembled-metres')for pid in old_ids]
            complete_check=check_constraint(dict(kind='shared-interface',parameters=complete_parameters),*candidate)
            if not complete_check['passed']:raise ValueError('The complete transferred candidate interface does not pass')
            proof=verify_transfer(old,*original,*candidate,binding['parameters'])
            proof['complete_candidate_interface_check']=complete_check
            overrides[binding['id']]=binding['parameters'];endpoint_overrides[binding['id']]=candidate
            transfer_checks.append(dict(id=binding['id'],old_owners=old_ids,new_owners=[binding['first'],binding['second']],evidence=binding['evidence'],**proof))
    print(f'[install] source interface proofs: {len(rebinding_checks)} rebound, {len(transfer_checks)} transferred',flush=True)
    constraints=check(db,candidates=changed_geometry,parameter_overrides=overrides,endpoint_overrides=endpoint_overrides)if changed_geometry else dict(passed=True,constraints=[])
    constraints['source_bound_interface_rebindings']=rebinding_checks
    constraints['source_bound_interface_transfers']=transfer_checks
    print(f'[install] independently checking {len(spec.get("interface_checks",[]))} declared candidate interfaces',flush=True)
    for ordinal,interface in enumerate(spec.get('interface_checks',[])):
        def get(part):
            if part in mesh_by_part:return mesh_by_part[part]
            return reader.read_part(after,ap[part],frame='model:'+model+':assembled-metres')
        result=check_constraint(interface,get(interface['first']),get(interface['second'])if interface.get('second')else None)
        constraints['constraints'].append(result);constraints['passed']=constraints['passed']and result['passed']
        if (ordinal+1)%25==0:print(f'[install] candidate interfaces {ordinal+1}/{len(spec["interface_checks"])}',flush=True)
    if not constraints['passed']:
        failure=CACHE/'install-failures'/(identity(model,digest_file(proposed),spec['job_id'])+'.json')
        failure.parent.mkdir(parents=True,exist_ok=True)
        atomic_json(failure,dict(passed=False,model=model,proposed_catalogue=str(proposed),
                    proposed_sha256=digest_file(proposed),constraint_checks=constraints,job_id=spec['job_id'],
                    runtime_report=os.environ.get('ATLAS_RUNTIME_REPORT')))
        raise ValueError('Changed geometry violates or lacks a resolved adjoining-interface constraint; details: '+str(failure))
    verify_bindings({str(current):spec['expected_catalogue_sha256'],str(proposed):proposed_sha,str(invpath):inventory_sha,**inventory_validation['input_bindings'],**banks})
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True);target=output/'proposal.json';atomic_bytes(target,proposed.read_bytes())
    proposal=dict(schema='local-atlas-install/v1',model=model,sex=sex,current=str(current),expected_catalogue_sha256=digest_file(current),proposed_catalogue=str(target),
                   proposed_sha256=digest_file(target),changed_parts=changed,inventory=str(invpath),inventory_sha256=inventory_sha,assembly=assembly,
                  bank_bindings=banks,geometry_checks=checks,point_reference_checks=point_checks,constraint_checks=constraints,job_id=spec['job_id'],geometry_check_implementation=digest_file(Path(__file__).with_name('audit.py')),
                  validation_implementation=implementation_signature(),
                   inventory_validation=inventory_validation,inventory_bindings=inventory_validation['input_bindings'],
                   passed=True,installed=False,external_deployment=False)
    atomic_json(output/'install.json',proposal);return proposal


def apply(db,path,name):
    proposal=json.loads(Path(path).read_text());model=proposal['model'];current=Path(proposal['current']);job=job_row(db,proposal['job_id'])
    if job['state']!='reviewed':raise ValueError('A reviewed geometry job is required before local installation')
    bound_artifacts=artifacts(db,job['id']);bound_hashes={a['sha256']for a in bound_artifacts}
    if not set(proposal['bank_bindings'].values())<=bound_hashes:raise ValueError('Review job does not bind every actual encoded replacement bank')
    verify_bindings({str(current):proposal['expected_catalogue_sha256'],proposal['proposed_catalogue']:proposal['proposed_sha256'],proposal['inventory']:proposal['inventory_sha256'],**proposal['bank_bindings']})
    if proposal.get('inventory_validation',{}).get('schema')!='authoritative-inventory-integrity/v1':
        raise ValueError('Installation needs the authoritative inventory integrity gate; prepare the proposal again')
    if proposal.get('inventory_bindings')!=proposal['inventory_validation'].get('input_bindings'):
        raise ValueError('Installation inventory input bindings are missing or changed; reprepare the proposal')
    verify_bindings(proposal.get('inventory_bindings',{}))
    from .inventory_validation import validate_authoritative_inventory
    inventory_validation=validate_authoritative_inventory(json.loads(Path(proposal['inventory']).read_text()),model=model,current_catalogue=current)
    if identity(inventory_validation)!=identity(proposal['inventory_validation']):
        raise ValueError('Authoritative inventory integrity changed after preparation; reprepare the proposal')
    if digest_file(Path(__file__).with_name('audit.py'))!=proposal['geometry_check_implementation']:raise ValueError('Geometry validation code changed; prepare the proposal again')
    if proposal.get('validation_implementation'):
        if implementation_signature()!=proposal['validation_implementation']:raise ValueError('Installation validation changed; prepare the proposal again')
    if not(ROOT/'.local-models/ta98-authority.sqlite').exists():raise ValueError('Local TA98 authority cache is required for an offline derived rebuild')
    after=json.loads(Path(proposal['proposed_catalogue']).read_bytes());catalogue_diff(json.loads(current.read_bytes()),after,proposal['changed_parts'])
    lock=CACHE/'install.lock';lock.parent.mkdir(parents=True,exist_ok=True)
    with lock.open('a')as f:
        fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        verify_bindings({str(current):proposal['expected_catalogue_sha256']})
        script=ROOT/'scripts/local_catalogue_revision.py';module_spec=importlib.util.spec_from_file_location('atlas_local_revision',script);module=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
        def mutate(sex,a,before):
            if sex!=proposal['sex']:return dict(unchanged=True)
            if identity(before)!=identity(json.loads(current.read_bytes())):raise ValueError('Catalogue changed during installation')
            a.clear();a.update(copy.deepcopy(after));return dict(proposal_sha256=digest_file(path),changed_parts=proposal['changed_parts'])
        proof=module.install_revision(name,mutate)
        manifest=ROOT/'.local-models/ta98-review'/name/'manifest.json';installed=digest_file(current);iid=identity('installation',model,installed,job['id'])
        from psycopg.types.json import Jsonb
        records=json.loads(manifest.read_text());installed_files={str(ROOT/r['path']):digest_file(ROOT/r['path'])for r in records}
        # Journal the exact published state before touching database status.
        atomic_json(manifest.parent/'published.json',dict(installation_id=iid,files=installed_files,job_id=job['id']))
        try:
            with db.connect()as c:
                c.execute('INSERT INTO atlas_geom.installations(id,job_id,model_id,previous_catalogue_sha256,catalogue_sha256,manifest_path,evidence)VALUES(%s,%s,%s,%s,%s,%s,%s)',(iid,job['id'],model,proposal['expected_catalogue_sha256'],installed,str(manifest),Jsonb(dict(proposal=str(Path(path).resolve()),proof=proof,installed_files=installed_files))))
                c.execute("UPDATE atlas_geom.jobs SET state='installed',updated_at=now()WHERE id=%s",(job['id'],))
        except BaseException:
            verify_bindings(installed_files)
            for r in records:
                atomic_bytes(ROOT/r['path'],gzip.decompress((ROOT/r['backup']).read_bytes()),overwrite=True)
            raise
    return dict(id=iid,installed=True,model=model,catalogue_sha256=installed,rollback_manifest=str(manifest),external_deployment=False)


def rollback(db,installation_id):
    with db.connect()as c:r=c.execute('SELECT * FROM atlas_geom.installations WHERE id=%s',(installation_id,)).fetchone()
    if not r:raise ValueError('Unknown local installation')
    verify_bindings(r['evidence']['installed_files'])
    records=json.loads(Path(r['manifest_path']).read_text());changes={};expected={}
    for item in records:
        target=ROOT/item['path'];backup=ROOT/item['backup'];verify_bindings({str(backup):item['sha256']});changes[target]=gzip.decompress(backup.read_bytes());expected[str(target)]=digest_file(target)
    # Compare against recorded installed revisions first; never overwrite an
    # atlas edited since this install. Derived files are journaled on rollback.
    sex={'local-ta98':'male','local-female-ta98':'female'}[r['model_id']];current=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json'
    if digest_file(current)!=r['catalogue_sha256']:raise ValueError('Atlas changed since installation; merge rather than roll back')
    return transaction(changes,CACHE/'rollbacks'/(installation_id+'-'+str(time.time_ns())),expected=expected)


def add_parser(sub):
    p=sub.add_parser('install');s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('prepare');q.add_argument('spec');q.add_argument('--output',required=True)
    q=s.add_parser('apply');q.add_argument('proposal');q.add_argument('--name',required=True)
    q=s.add_parser('rollback');q.add_argument('id')
    p.set_defaults(handler=run)


def run(args,db):
    if args.action=='prepare':return prepare(db,json.loads(Path(args.spec).read_text()),args.output)
    if args.action=='apply':return apply(db,args.proposal,args.name)
    return rollback(db,args.id)
