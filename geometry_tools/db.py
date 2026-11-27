"""Local PostgreSQL catalogue, spatial/shape indexes and evidence cache."""
from __future__ import annotations
from pathlib import Path
from contextlib import contextmanager
import os, json
import psycopg
from psycopg.rows import dict_row
from psycopg.types.json import Jsonb
from psycopg.conninfo import conninfo_to_dict
from .common import CACHE, VERSION, digest_file, identity


class Database:
    def __init__(self, dsn=None):
        self.dsn=dsn or os.environ.get('ATLAS_GEOMETRY_DSN','dbname=human_atlas_geometry host=/tmp connect_timeout=5')
        info=conninfo_to_dict(self.dsn);host=info.get('host',os.environ.get('PGHOST',''))
        if info.get('service')or os.environ.get('PGSERVICE'):raise ValueError('Use an explicit local connection rather than a service that could redirect remotely')
        if info.get('hostaddr',os.environ.get('PGHOSTADDR',''))not in ['', '127.0.0.1','::1']:raise ValueError('Only local PostgreSQL host addresses are allowed')
        if host not in ['', 'localhost','127.0.0.1','::1'] and not host.startswith('/'):
            raise ValueError('The geometry workbench connects only to local PostgreSQL')

    @contextmanager
    def connect(self):
        with psycopg.connect(self.dsn,row_factory=dict_row) as conn:
            yield conn

    def initialize(self):
        with self.connect() as c:
            c.execute("SELECT pg_advisory_xact_lock(hashtext('human-atlas-geometry-schema-v1'))")
            c.execute(Path(__file__).with_name('schema.sql').read_text())
        return self.status()

    def status(self):
        with self.connect() as c:
            row=c.execute("SELECT current_database() AS database, current_setting('server_version') AS version").fetchone()
            row['extensions']=c.execute("SELECT extname,extversion FROM pg_extension WHERE extname IN ('vector','cube') ORDER BY extname").fetchall()
            row['meshes']=c.execute('SELECT count(*) AS total,count(*) FILTER(WHERE error IS NOT NULL)AS errors FROM atlas_geom.meshes').fetchone()
            row['current_models']=c.execute('SELECT * FROM atlas_geom.current_catalogues ORDER BY model_id').fetchall()
            return row

    def mesh(self, reference, model=None):
        with self.connect() as c:
            exact=c.execute('SELECT * FROM atlas_geom.meshes WHERE id=%s',(reference,)).fetchone()
            if exact:return exact
            rows=c.execute('''SELECT m.* FROM atlas_geom.meshes m JOIN atlas_geom.current_catalogues a
                              ON a.catalogue_sha256=m.catalogue_sha256 WHERE m.part_id=%s AND (%s::text IS NULL OR a.model_id=%s)''',(reference,model,model)).fetchall()
            if len(rows)!=1:raise ValueError(f'Mesh reference resolves to {len(rows)} current instances; provide a model or instance ID')
            self.assert_current_catalogues(model)
            return rows[0]

    def assert_current_catalogues(self,model=None):
        with self.connect()as c:
            cats=c.execute('SELECT c.path,c.metadata FROM atlas_geom.catalogues c JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=c.sha256 WHERE(%s::text IS NULL OR a.model_id=%s)',(model,model)).fetchall()
        for cat in cats:
            if not Path(cat['path']).is_file()or digest_file(cat['path'])!=cat['metadata']['json_sha256']:
                raise ValueError('Current catalogue changed; registry-import the actual revision before resolving current items')

    def cache_get(self, key):
        with self.connect() as c:row=c.execute('SELECT * FROM atlas_geom.cache WHERE key=%s',(key,)).fetchone()
        if row and Path(row['path']).is_file() and digest_file(row['path'])==row['sha256']:
            return row
        return None

    def cache_put(self, *, key, geometry_sha, kind, parameters, path, summary):
        path=Path(path).resolve()
        with self.connect() as c:
            c.execute('''INSERT INTO atlas_geom.cache(key,geometry_sha256,kind,tool_version,parameters,path,sha256,summary)
                         VALUES(%s,%s,%s,%s,%s,%s,%s,%s) ON CONFLICT(key) DO UPDATE SET
                         path=EXCLUDED.path,sha256=EXCLUDED.sha256,summary=EXCLUDED.summary''',
                      (key,geometry_sha,kind,VERSION,Jsonb(parameters),str(path),digest_file(path),Jsonb(summary)))

    def ensure_geometry(self,mesh):
        with self.connect()as c:
            found=c.execute('SELECT 1 FROM atlas_geom.geometry WHERE sha256=%s',(mesh.geometry_sha,)).fetchone()
            if not found:
                from .registry import geometry_record
                self.put_geometry(c,geometry_record(mesh))

    def register_artifact_mesh(self, mesh, path, *, part_id, name, metadata=None):
        from .registry import geometry_record
        path=Path(path).resolve();geometry=geometry_record(mesh);gid=mesh.geometry_sha
        mid=identity('artifact-instance/v1',digest_file(path),part_id,mesh.frame)
        source={'kind':'npz','path':str(path),'sha256':digest_file(path)}
        with self.connect() as c:
            c.execute("INSERT INTO atlas_geom.frames(id,units,metadata)VALUES(%s,'m',%s)ON CONFLICT DO NOTHING",(mesh.frame,Jsonb({'role':'explicit artifact coordinates'})))
            self.put_geometry(c,geometry)
            bounds=mesh.bounds.tolist()
            c.execute('''INSERT INTO atlas_geom.meshes(id,part_id,name,geometry_sha256,attributes_sha256,frame_id,bbox,source,metadata)
                         VALUES(%s,%s,%s,%s,%s,%s,cube(%s::float8[],%s::float8[]),%s,%s) ON CONFLICT(id)DO NOTHING''',
                      (mid,part_id,name,gid,mesh.attributes_sha,mesh.frame,*bounds,Jsonb(source),Jsonb(metadata or {})))
        return mid

    @staticmethod
    def put_geometry(c, row):
        c.execute('''INSERT INTO atlas_geom.geometry(sha256,surface_sha256,vertex_count,face_count,metrics,descriptor,descriptor_version)
                     VALUES(%s,%s,%s,%s,%s,%s::vector,%s) ON CONFLICT(sha256)DO NOTHING''',
                  (row['sha256'],row['surface_sha256'],row['vertex_count'],row['face_count'],Jsonb(row['metrics']),
                   None if row['descriptor']is None else json.dumps(row['descriptor']),row['descriptor_version']))
