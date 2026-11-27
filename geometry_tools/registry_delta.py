"""Register a local installed revision using hash-verified unchanged buffers.

This is an alternative saved stage after a full re-index exceeds its deadline.
It reuses a previous geometry only when all indexed buffer fields, its frame,
bank path and complete bank hash match. Changed buffers are decoded normally.
The current revision is published only after the complete transaction succeeds.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from psycopg.types.json import Jsonb
from .common import identity,digest_file,atomic_json
from .db import Database
from .model import AtlasReader
from .registry import catalogue_revision_identity,geometry_record,side_from_record

LAYOUT=('vertexCount','indexCount','positions','indices','normals','normalSigns','uvs','colors','geometryOrigin','sourceOffset','geometryScale','geometryInstance')
def identical_layout(first,second):
    return all(first.get(k)==second.get(k)for k in LAYOUT)

def reusable_geometry(previous,record,frame,bank_path,bank_sha256):
    return bool(previous and not previous['error']and previous['frame_id']==frame
                and previous['bank_path']==bank_path and previous['bank_sha256']==bank_sha256
                and identical_layout(previous['record'],record))

def register_delta(db,path,model,previous_json_sha,*,frame_id=None):
    path=Path(path).resolve();data=json.loads(path.read_bytes());json_sha=digest_file(path);frame=frame_id or'model:'+model+':assembled-metres';reader=AtlasReader();started=time.monotonic()
    banks={str(reader.asset(data['chunks'][i]['url'])):reader.file_sha(reader.asset(data['chunks'][i]['url']))for i in sorted({p['chunk']for p in data['parts']})}
    revision=catalogue_revision_identity(json_sha,banks,frame,model);hidden=set(data.get('defaultHidden',[]));concepts={c['id']:c for c in data['concepts']};reused=0;decoded=0
    with db.connect()as c:
        c.execute('SELECT pg_advisory_xact_lock(hashtext(%s))',('atlas-geometry-import:'+model,))
        previous=c.execute('SELECT q.* FROM atlas_geom.catalogues q JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=q.sha256 WHERE a.model_id=%s',(model,)).fetchone()
        if not previous or previous['metadata']['json_sha256']!=previous_json_sha or previous['frame_id']!=frame:raise ValueError('The previous current revision does not match the explicit delta source')
        rows=c.execute('SELECT id,part_id,frame_id,source->\'record\'AS record,source->>\'bank_path\'AS bank_path,source->>\'bank_sha256\'AS bank_sha256,error FROM atlas_geom.meshes WHERE catalogue_sha256=%s',(previous['sha256'],)).fetchall();byid={r['part_id']:r for r in rows}
        if len({p['id']for p in data['parts']})!=len(data['parts']):raise ValueError('Duplicate part IDs')
        c.execute('INSERT INTO atlas_geom.catalogues(sha256,model_id,path,sex,frame_id,metadata)VALUES(%s,%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',(revision,model,str(path),data.get('sex'),frame,Jsonb(dict(part_count=len(data['parts']),name=data.get('name'),source_units='m',json_sha256=json_sha,bank_sha256=banks,registration_verified=False))))
        for p in data['parts']:
            old=byid.get(p['id']);bank=str(reader.asset(data['chunks'][p['chunk']]['url']));mid=identity('catalogue-instance/v1',revision,p['id']);term=p.get('ta98Term')or concepts.get(p.get('conceptId'),{}).get('ta98Term');metadata=dict(record=p,concept=concepts.get(p.get('conceptId')),registration_verified=False)
            if reusable_geometry(old,p,frame,bank,banks[bank]):
                c.execute('''INSERT INTO atlas_geom.meshes(id,catalogue_sha256,part_id,concept_id,name,ta98_term,sex,side,system,source_id,license,hidden,suppressed,geometry_sha256,attributes_sha256,frame_id,bbox,source,metadata,error)
                    SELECT %s,%s,part_id,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,geometry_sha256,attributes_sha256,frame_id,bbox,
                    source||%s::jsonb,%s,NULL FROM atlas_geom.meshes WHERE id=%s ON CONFLICT(id)DO NOTHING''',
                    (mid,revision,p.get('conceptId'),p['name'],term,data.get('sex'),side_from_record(p),p.get('system'),p.get('sourceId'),p.get('license'),bool(p.get('defaultHidden')or p['id']in hidden),bool(p.get('suppressed')),Jsonb(dict(record=p,catalogue_path=str(path),catalogue_sha256=json_sha)),Jsonb(metadata),old['id']))
                reused+=1
            else:
                m=reader.read_part(data,p,frame=frame,catalogue_path=path,catalogue_sha=json_sha)
                if not len(m.faces):raise ValueError('A changed non-mesh needs the full registry path')
                if not c.execute('SELECT 1 FROM atlas_geom.geometry WHERE sha256=%s',(m.geometry_sha,)).fetchone():db.put_geometry(c,geometry_record(m))
                bounds=m.bounds
                c.execute('''INSERT INTO atlas_geom.meshes(id,catalogue_sha256,part_id,concept_id,name,ta98_term,sex,side,system,source_id,license,hidden,suppressed,geometry_sha256,attributes_sha256,frame_id,bbox,source,metadata,error)
                    VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,cube(%s::float8[],%s::float8[]),%s,%s,NULL)ON CONFLICT(id)DO NOTHING''',
                    (mid,revision,p['id'],p.get('conceptId'),p['name'],term,data.get('sex'),side_from_record(p),p.get('system'),p.get('sourceId'),p.get('license'),bool(p.get('defaultHidden')or p['id']in hidden),bool(p.get('suppressed')),m.geometry_sha,m.attributes_sha,frame,bounds[0].tolist(),bounds[1].tolist(),Jsonb(m.source),Jsonb(metadata)))
                decoded+=1
        count=c.execute('SELECT count(*)AS n FROM atlas_geom.meshes WHERE catalogue_sha256=%s',(revision,)).fetchone()['n'];assert count==len(data['parts'])
        if digest_file(path)!=json_sha or any(reader.file_sha(Path(p))!=h for p,h in banks.items()):
            raise ValueError('Catalogue or source bank changed during delta registration')
        c.execute('INSERT INTO atlas_geom.current_catalogues VALUES(%s,%s)ON CONFLICT(model_id)DO UPDATE SET catalogue_sha256=EXCLUDED.catalogue_sha256',(model,revision))
    return dict(model_id=model,revision_sha256=revision,json_sha256=json_sha,parts=count,reused_literal_buffer_records=reused,decoded_changed_parts=decoded,all_bank_hashes_verified=True,geometry_modified=False,errors=[],elapsed_seconds=time.monotonic()-started)

def main():
    p=argparse.ArgumentParser();p.add_argument('catalogue',type=Path);p.add_argument('--model',required=True);p.add_argument('--frame',help='Explicit existing assembled frame; must match the previous registry revision');p.add_argument('--previous-json-sha',required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args();r=register_delta(Database(),a.catalogue,a.model,a.previous_json_sha,frame_id=a.frame);atomic_json(a.report,r);print(json.dumps(r))

if __name__=='__main__':main()
