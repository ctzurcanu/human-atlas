"""Incremental, content-addressed catalogue ingestion and local shape retrieval."""
from pathlib import Path
import json,re,time
import numpy as np
from psycopg.types.json import Jsonb
from .common import identity,digest_file
from .model import AtlasReader,Mesh

DESCRIPTOR_VERSION='surface-area-d2/v1-64'


def file_backed_source(source):
    """Retain the authoritative bank/record binding, derive index maps on load.

    AtlasReader.reload reconstructs both maps from the checked index bank.
    Repeating their potentially million-element lists in every registry row
    defeats file-backed storage and is unnecessary for exact source replay.
    Non-atlas artifacts retain their complete, format-specific provenance.
    """
    result=dict(source)
    if result.get('kind')=='atlas-part':
        deferred={key:len(result.pop(key))for key in ['source_vertex_ids','source_face_ids']if key in result}
        if deferred:result['derived_index_maps']=dict(format='checked-bank-reader/v1',counts=deferred)
    return result


def shape_descriptor(mesh, samples=512):
    """Deterministic surface descriptor; retrieval only, never anatomical proof.

    Area-weighted radial/D2/normal-angle distributions are translation,
    rotation, reflection and uniform-scale invariant up to sampling arithmetic.
    Face order remains part of the mesh identity, not a shape similarity claim.
    """
    t=mesh.triangles
    cross=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);areas=np.linalg.norm(cross,axis=1)/2
    total=float(areas.sum())
    if not total>0:return None
    rng=np.random.default_rng(1729);ids=rng.choice(len(t),samples,p=areas/total)
    u=rng.random((samples,2));u[:,0]=np.sqrt(u[:,0]);weights=np.c_[1-u[:,0],u[:,0]*(1-u[:,1]),u[:,0]*u[:,1]]
    points=np.einsum('ni,nij->nj',weights,t[ids]);center=np.einsum('i,ij->j',areas,t.mean(1))/total
    q=points-center;rms=float(np.sqrt(np.mean(np.sum(q*q,axis=1))))
    if rms==0:return None
    q/=rms
    radial=np.histogram(np.linalg.norm(q,axis=1),bins=np.linspace(0,4,17))[0].astype(float)/samples
    shift=np.roll(q,samples//2,axis=0);d2=np.histogram(np.linalg.norm(q-shift,axis=1),bins=np.linspace(0,6,33))[0].astype(float)/samples
    n=cross[ids];n/=np.linalg.norm(n,axis=1)[:,None]
    angle=np.histogram(np.einsum('ij,ij->i',n,np.roll(n,samples//2,axis=0)),bins=np.linspace(-1.0000001,1.0000001,9))[0].astype(float)/samples
    eig=np.linalg.eigvalsh(np.cov(q.T));eig/=eig.sum()
    feature=np.r_[eig,total/(rms*rms*100),float(np.quantile(np.linalg.norm(q,axis=1),.1))/4,
                  float(np.quantile(np.linalg.norm(q,axis=1),.9))/4,
                  float(np.linalg.norm(q,axis=1).max())/4, float(np.mean(np.linalg.norm(q-shift,axis=1)))/6]
    result=np.r_[radial,d2,angle,feature];assert result.shape==(64,)
    result/=np.linalg.norm(result)
    return result.tolist()


def geometry_record(mesh):
    t=mesh.triangles;area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    origin=t.reshape(-1,3).mean(0)if len(t)else np.zeros(3);p=t-origin
    volume=float(np.einsum('ij,ij->i',p[:,0],np.cross(p[:,1],p[:,2])).sum()/6)if len(t)else 0.
    return dict(sha256=mesh.geometry_sha,surface_sha256=mesh.surface_sha,vertex_count=len(mesh.vertices),face_count=len(mesh.faces),
                metrics=dict(bounds_m=mesh.bounds.tolist()if len(t)else None,surface_area_m2=float(area.sum()),oriented_triangle_integral_m3=volume,
                             volume_valid_for_closed_oriented_mesh_only=True,zero_area_faces=int(np.sum(area==0)),topology_checked=False,anatomically_accepted=False),
                descriptor=shape_descriptor(mesh)if len(t)else None,descriptor_version=DESCRIPTOR_VERSION)


def side_from_record(part):
    # Laterality is semantic metadata. Unknown remains unknown; coordinates do
    # not silently decide anatomical left/right across reflected donors.
    value=part.get('side')or part.get('laterality')
    if value in ['left','right']:return value
    found=re.search(r'\((left|right)\)\s*$',part.get('name',''),re.I)
    return found.group(1).lower()if found else None


def catalogue_revision_identity(json_sha, bank_hashes, frame_id, model_id):
    """Keep identical file revisions distinct when their owning models differ."""
    return identity('catalogue-revision/v2', model_id, json_sha, bank_hashes, frame_id)


def import_catalogue(db,path,*,model_id,frame_id=None,root=None,progress=None,max_parts=None):
    if max_parts is not None and (not isinstance(max_parts,int)or max_parts<1):
        raise ValueError('A staged import needs a positive maximum part count')
    path=Path(path).resolve();json_sha=digest_file(path);a=json.loads(path.read_bytes());frame_id=frame_id or'model:'+model_id+':assembled-metres'
    reader=AtlasReader(root)if root else AtlasReader();hidden=set(a.get('defaultHidden',[]));byconcept={c['id']:c for c in a.get('concepts',[])}
    if len({p['id']for p in a['parts']})!=len(a['parts']):raise ValueError('Duplicate source part IDs must be repaired explicitly')
    # A bank can change without its catalogue JSON changing. Bind every actual
    # referenced asset, so an unchanged import cannot silently publish stale
    # geometry or attribute caches under the same immutable revision.
    bank_hashes={str(reader.asset(a['chunks'][ci]['url'])):reader.file_sha(reader.asset(a['chunks'][ci]['url']))for ci in sorted({p['chunk']for p in a['parts']})}
    sha=catalogue_revision_identity(json_sha,bank_hashes,frame_id,model_id)
    began=time.monotonic();done=0;errors=[];reused=0;geometries={}
    with db.connect() as c:
        c.execute('SELECT pg_advisory_xact_lock(hashtext(%s))',('atlas-geometry-import:'+model_id,))
        imported=c.execute('SELECT part_id,error FROM atlas_geom.meshes WHERE catalogue_sha256=%s',(sha,)).fetchall()
        existing={r['part_id']for r in imported};expected={p['id']for p in a['parts']}
        if not existing<=expected:raise ValueError('Staged revision contains unexpected source identities')
        done=len(existing);errors=[dict(part_id=r['part_id'],error=r['error'])for r in imported if r['error']]
        if existing==expected:
            c.execute('INSERT INTO atlas_geom.current_catalogues VALUES(%s,%s)ON CONFLICT(model_id)DO UPDATE SET catalogue_sha256=EXCLUDED.catalogue_sha256',(model_id,sha))
            return dict(model_id=model_id,revision_sha256=sha,json_sha256=json_sha,parts=done,unchanged=True,complete=True,remaining_parts=0,parts_imported_this_run=0,geometry_modified=False,elapsed_seconds=time.monotonic()-began)
        c.execute("INSERT INTO atlas_geom.frames VALUES(%s,'m',%s)ON CONFLICT DO NOTHING",(frame_id,Jsonb(dict(model_id=model_id,meaning='Current assembled render coordinates; source provenance is separate and anatomical registration remains unverified'))))
        c.execute('INSERT INTO atlas_geom.catalogues(sha256,model_id,path,sex,frame_id,metadata)VALUES(%s,%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',
                  (sha,model_id,str(path),a.get('sex'),frame_id,Jsonb(dict(part_count=len(a['parts']),name=a.get('name'),source_units='m',json_sha256=json_sha,bank_sha256=bank_hashes,registration_verified=False))))
        staged=0
        for part in a['parts']:
            if part['id']in existing:continue
            if max_parts is not None and staged>=max_parts:break
            source={};error=None;geometry_sha=None;attrs_sha=None;bounds=None
            try:
                mesh=reader.read_part(a,part,frame=frame_id,catalogue_path=path,catalogue_sha=json_sha)
                source=file_backed_source(mesh.source)
                if len(mesh.faces):
                    geometry_sha=mesh.geometry_sha;attrs_sha=mesh.attributes_sha;bounds=mesh.bounds.tolist()
                    if geometry_sha not in geometries:
                        exists=c.execute('SELECT 1 FROM atlas_geom.geometry WHERE sha256=%s',(geometry_sha,)).fetchone()
                        if not exists:db.put_geometry(c,geometry_record(mesh))
                        else:reused+=1
                        geometries[geometry_sha]=True
                else:error='No indexed triangles; retained as an explicitly non-mesh source record'
            except (OSError,ValueError,KeyError,IndexError)as exc:
                error=type(exc).__name__+': '+str(exc);errors.append(dict(part_id=part['id'],error=error))
            term=part.get('ta98Term')or byconcept.get(part.get('conceptId'),{}).get('ta98Term')
            metadata=dict(record=part,concept=byconcept.get(part.get('conceptId')),registration_verified=False)
            mid=identity('catalogue-instance/v1',sha,part['id'])
            c.execute('''INSERT INTO atlas_geom.meshes(id,catalogue_sha256,part_id,concept_id,name,ta98_term,sex,side,system,source_id,license,hidden,suppressed,
                         geometry_sha256,attributes_sha256,frame_id,bbox,source,metadata,error)
                         VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,
                         CASE WHEN %s::float8[]IS NULL THEN NULL ELSE cube(%s::float8[],%s::float8[])END,%s,%s,%s)
                         ON CONFLICT(id)DO UPDATE SET source=EXCLUDED.source,error=EXCLUDED.error''',
                      (mid,sha,part['id'],part.get('conceptId'),part.get('name',part['id']),term,a.get('sex'),side_from_record(part),part.get('system'),part.get('sourceId'),part.get('license'),
                       bool(part.get('defaultHidden')or part['id']in hidden),bool(part.get('suppressed')),geometry_sha,attrs_sha,frame_id,
                       None if bounds is None else bounds[0],None if bounds is None else bounds[0],None if bounds is None else bounds[1],Jsonb(source),Jsonb(metadata),error))
            done+=1;staged+=1
            if progress and(done%100==0 or done==len(a['parts'])):progress(dict(imported=done,total=len(a['parts']),errors=len(errors),elapsed_seconds=round(time.monotonic()-began,2)))
        # Partial batches are durable but never become the current registry.
        # A later bounded call resumes only this exact JSON/bank/frame revision.
        if digest_file(path)!=json_sha or any(reader.file_sha(Path(p))!=h for p,h in bank_hashes.items()):
            raise ValueError('Catalogue or source bank changed during registry import')
        complete=done==len(a['parts'])
        if complete:c.execute('INSERT INTO atlas_geom.current_catalogues VALUES(%s,%s)ON CONFLICT(model_id)DO UPDATE SET catalogue_sha256=EXCLUDED.catalogue_sha256',(model_id,sha))
    return dict(model_id=model_id,revision_sha256=sha,json_sha256=json_sha,parts=done,distinct_geometry_records=len(geometries),reused_geometry_records=reused,errors=errors,
                complete=complete,remaining_parts=len(a['parts'])-done,parts_imported_this_run=staged,current_catalogue_updated=complete,
                geometry_modified=False,source_arrays_copied_to_database=False,elapsed_seconds=time.monotonic()-began)


def load_registered(db,reference,model=None,reader=None):
    row=db.mesh(reference,model)
    if row['error']or not row['geometry_sha256']:raise ValueError('This source has no usable triangle mesh: '+str(row['error']))
    if row['source'].get('kind')=='npz':
        source=row['source']
        if digest_file(source['path'])!=source['sha256']:raise ValueError('Registered artifact changed')
        mesh=Mesh.load(source['path'],frame=row['frame_id'])
    else:
        source_root=Path(row['source']['project_root']).resolve()
        reader=reader if reader is not None and reader.root==source_root else AtlasReader(source_root)
        mesh=reader.reload(row['source'],row['frame_id'])
    if mesh.geometry_sha!=row['geometry_sha256']or mesh.attributes_sha!=row['attributes_sha256']:
        raise ValueError('Decoded geometry/attributes differ from the registered source identity')
    return row,mesh


def search(db,*,query=None,similar=None,model=None,sex=None,side=None,term=None,system=None,license_name=None,include_hidden=False,limit=20):
    if not 1<=limit<=500:raise ValueError('Search limit must be between 1 and 500')
    db.assert_current_catalogues(model)
    params=[];where=['m.geometry_sha256 IS NOT NULL','m.error IS NULL']
    if not include_hidden:where+=['NOT m.hidden','NOT m.suppressed']
    for col,value in [('a.model_id',model),('m.sex',sex),('m.side',side),('m.ta98_term',term),('m.system',system),('m.license',license_name)]:
        if value is not None:where.append(col+'=%s');params.append(value)
    if query:
        where.append("(to_tsvector('simple',m.name)@@plainto_tsquery('simple',%s)OR m.part_id ILIKE %s OR m.name ILIKE %s)");params.extend([query,'%'+query+'%','%'+query+'%'])
    with db.connect()as c:
        select='NULL::float8 AS shape_distance';order='m.name,m.id'
        if similar:
            src=db.mesh(similar,model);vector=c.execute('SELECT descriptor::text AS v FROM atlas_geom.geometry WHERE sha256=%s',(src['geometry_sha256'],)).fetchone()
            if not vector or vector['v']is None:raise ValueError('Source mesh has no geometric descriptor')
            # Exact filtered ranking: filter before vector ordering, avoiding
            # ANN post-filter recall loss for small anatomical cohorts.
            select='g.descriptor <=> %s::vector AS shape_distance';params=[vector['v'],*params];order='shape_distance,m.id'
            where.append('g.descriptor IS NOT NULL')
        sql='''SELECT m.id,m.part_id,m.name,m.sex,m.side,m.ta98_term,m.license,m.frame_id,m.hidden,m.suppressed,m.geometry_sha256,'''+select+'''
               FROM atlas_geom.meshes m JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=m.catalogue_sha256
               JOIN atlas_geom.geometry g ON g.sha256=m.geometry_sha256 WHERE '''+' AND '.join(where)+' ORDER BY '+order+' LIMIT %s'
        rows=c.execute(sql,[*params,limit]).fetchall()
    return dict(results=rows,descriptor_version=DESCRIPTOR_VERSION,shape_similarity_is_not_registration_or_anatomical_acceptance=True)


def inventory(db,root_term,model,*,authority=None,classifications=None,sex=None):
    """Read every authoritative descendant; unclassified/missing stay explicit."""
    import sqlite3
    from .common import ROOT
    authority=Path(authority or ROOT/'.local-models/ta98-authority.sqlite').resolve()
    if not authority.is_file():raise ValueError('Local authoritative terminology cache is required; this command never downloads or invents terms')
    with sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True)as conn:
        conn.row_factory=sqlite3.Row;terms={r['id']:dict(r)for r in conn.execute('SELECT * FROM ta98')}
    from .ta98_authority import corrected_terms,CONTEXT_REFERENCES,term_sex
    terms=corrected_terms(terms)
    sex=sex or {'local-ta98':'male','local-female-ta98':'female'}.get(model)
    if sex not in [None,'male','female']:raise ValueError('Model sex must be explicit male or female')
    if root_term not in terms:raise ValueError('Unknown authoritative term')
    descendants={root_term}
    while True:
        more={code for code,t in terms.items()if t['parent_id']in descendants}
        if more<=descendants:break
        descendants|=more
    context=CONTEXT_REFERENCES.get(root_term,{})
    if any(code not in terms for code in context):raise ValueError('An explicit TA98 context reference is missing from the authority')
    descendants.update(context)
    with db.connect()as c:
        rows=c.execute('''SELECT m.id,m.part_id,m.name,m.ta98_term,m.metadata,m.geometry_sha256,m.hidden,m.suppressed FROM atlas_geom.meshes m
          JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=m.catalogue_sha256 WHERE a.model_id=%s AND m.ta98_term=ANY(%s)''',(model,sorted(descendants))).fetchall()
    byterm={}
    for row in rows:byterm.setdefault(row['ta98_term'],[]).append({k:v for k,v in row.items()if k!='metadata'}|{'mapping_kind':(row['metadata'].get('concept')or{}).get('ta98Kind')})
    classes=classifications or{};allowed={'tissue-region','surface-region','wall-layer','opening','boundary-curve','point-landmark','microscopic-structure','whole-organ','luminal-content','unclassified'}
    result=[]
    for code in sorted(descendants):
        classification=classes.get(code,'unclassified')
        if classification not in allowed:raise ValueError('Invalid anatomical representation class')
        qualified=term_sex(terms,code);applicable=not sex or not qualified or sex==qualified
        result.append(dict(authority=terms[code],representation_class=classification,source_meshes=byterm.get(code,[]),
                           applicable_to_model=applicable,authority_sex=qualified,
                           **(dict(contextual_cross_reference=dict(owner=root_term,source=context[code],canonical_parent_preserved=True))if code in context else{}),
                           distinct_representation_verified=False,missing_distinct_representation=applicable))
    return dict(schema='ta98-geometry-inventory/v1',model=model,root_term=root_term,authority_path=str(authority),authority_sha256=digest_file(authority),
                authority_corrections_path=str(Path(__file__).with_name('ta98_authority.py')),authority_corrections_sha256=digest_file(Path(__file__).with_name('ta98_authority.py')),
                terms=result,term_count=len(result),model_sex=sex,applicable_term_count=sum(row['applicable_to_model']for row in result),
                automatic_mapping_is_not_representation_acceptance=True,anatomical_acceptance=False)
