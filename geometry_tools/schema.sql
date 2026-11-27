CREATE EXTENSION IF NOT EXISTS vector;
CREATE EXTENSION IF NOT EXISTS cube;
CREATE SCHEMA IF NOT EXISTS atlas_geom;
CREATE TABLE IF NOT EXISTS atlas_geom.schema_versions (
 version integer PRIMARY KEY, installed_at timestamptz NOT NULL DEFAULT now()
);
INSERT INTO atlas_geom.schema_versions(version) VALUES(1) ON CONFLICT DO NOTHING;
CREATE TABLE IF NOT EXISTS atlas_geom.frames (
 id text PRIMARY KEY, units text NOT NULL CHECK(units='m'), metadata jsonb NOT NULL
);
CREATE TABLE IF NOT EXISTS atlas_geom.catalogues (
 sha256 text PRIMARY KEY CHECK(length(sha256)=64), model_id text NOT NULL, path text NOT NULL,
 sex text, frame_id text NOT NULL REFERENCES atlas_geom.frames(id), metadata jsonb NOT NULL,
 imported_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE IF NOT EXISTS atlas_geom.geometry (
 sha256 text PRIMARY KEY CHECK(length(sha256)=64), surface_sha256 text NOT NULL,
 vertex_count integer NOT NULL, face_count integer NOT NULL, metrics jsonb NOT NULL,
 descriptor vector(64), descriptor_version text NOT NULL, created_at timestamptz NOT NULL DEFAULT now()
);
CREATE INDEX IF NOT EXISTS geometry_surface_idx ON atlas_geom.geometry(surface_sha256);
CREATE INDEX IF NOT EXISTS geometry_shape_hnsw ON atlas_geom.geometry USING hnsw(descriptor vector_cosine_ops);
CREATE TABLE IF NOT EXISTS atlas_geom.meshes (
 id text PRIMARY KEY CHECK(length(id)=64), catalogue_sha256 text REFERENCES atlas_geom.catalogues(sha256),
 part_id text NOT NULL, concept_id text, name text NOT NULL, ta98_term text, sex text, side text,
 system text, source_id text, license text, hidden boolean NOT NULL DEFAULT false,
 suppressed boolean NOT NULL DEFAULT false, geometry_sha256 text REFERENCES atlas_geom.geometry(sha256),
 attributes_sha256 text, frame_id text NOT NULL REFERENCES atlas_geom.frames(id), bbox cube,
 source jsonb NOT NULL, metadata jsonb NOT NULL, error text,
 UNIQUE(catalogue_sha256,part_id)
);
CREATE INDEX IF NOT EXISTS mesh_spatial_gist ON atlas_geom.meshes USING gist(bbox);
CREATE INDEX IF NOT EXISTS mesh_filter_idx ON atlas_geom.meshes(catalogue_sha256,system,ta98_term,sex,side);
CREATE INDEX IF NOT EXISTS mesh_name_fts ON atlas_geom.meshes USING gin(to_tsvector('simple',name));
CREATE TABLE IF NOT EXISTS atlas_geom.current_catalogues (
 model_id text PRIMARY KEY, catalogue_sha256 text NOT NULL REFERENCES atlas_geom.catalogues(sha256)
);
CREATE TABLE IF NOT EXISTS atlas_geom.cache (
 key text PRIMARY KEY, geometry_sha256 text REFERENCES atlas_geom.geometry(sha256), kind text NOT NULL,
 tool_version text NOT NULL, parameters jsonb NOT NULL, path text NOT NULL, sha256 text NOT NULL,
 summary jsonb NOT NULL, created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE IF NOT EXISTS atlas_geom.assemblies (
 id text PRIMARY KEY, catalogue_sha256 text NOT NULL REFERENCES atlas_geom.catalogues(sha256),
 name text NOT NULL, basis text NOT NULL, metadata jsonb NOT NULL
);
CREATE TABLE IF NOT EXISTS atlas_geom.assembly_members (
 assembly_id text NOT NULL REFERENCES atlas_geom.assemblies(id), mesh_id text NOT NULL REFERENCES atlas_geom.meshes(id),
 active boolean NOT NULL, PRIMARY KEY(assembly_id,mesh_id)
);
CREATE TABLE IF NOT EXISTS atlas_geom.constraints (
 id text PRIMARY KEY, kind text NOT NULL CHECK(kind IN ('fixed','shared-interface','disjoint','nested','adjacent-unjoined','sliding','open-conduit','clearance')),
 first_mesh_id text NOT NULL REFERENCES atlas_geom.meshes(id), second_mesh_id text REFERENCES atlas_geom.meshes(id),
 evidence jsonb NOT NULL, parameters jsonb NOT NULL, created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE IF NOT EXISTS atlas_geom.relation_declarations (
 id text PRIMARY KEY, catalogue_sha256 text NOT NULL REFERENCES atlas_geom.catalogues(sha256),
 first_mesh_id text NOT NULL REFERENCES atlas_geom.meshes(id), second_mesh_id text NOT NULL REFERENCES atlas_geom.meshes(id),
 kind text NOT NULL, co_move boolean NOT NULL, evidence jsonb NOT NULL,
 geometric_constraint_inferred boolean NOT NULL DEFAULT false CHECK(NOT geometric_constraint_inferred)
);
CREATE INDEX IF NOT EXISTS relation_declarations_filter ON atlas_geom.relation_declarations(catalogue_sha256,kind,first_mesh_id,second_mesh_id);
INSERT INTO atlas_geom.schema_versions(version) VALUES(2) ON CONFLICT DO NOTHING;
CREATE TABLE IF NOT EXISTS atlas_geom.jobs (
 id text PRIMARY KEY, kind text NOT NULL, input_signature text NOT NULL, tool_version text NOT NULL,
 parameters jsonb NOT NULL, state text NOT NULL CHECK(state IN ('queued','running','passed','failed','reviewed','installed')),
 result jsonb, error text, created_at timestamptz NOT NULL DEFAULT now(), updated_at timestamptz NOT NULL DEFAULT now(),
 UNIQUE(kind,input_signature,tool_version)
);
CREATE TABLE IF NOT EXISTS atlas_geom.job_artifacts (
 job_id text NOT NULL REFERENCES atlas_geom.jobs(id), role text NOT NULL, path text NOT NULL, sha256 text NOT NULL,
 PRIMARY KEY(job_id,role,path)
);
CREATE TABLE IF NOT EXISTS atlas_geom.installations (
 id text PRIMARY KEY, job_id text REFERENCES atlas_geom.jobs(id), model_id text NOT NULL,
 previous_catalogue_sha256 text NOT NULL, catalogue_sha256 text NOT NULL,
 manifest_path text NOT NULL, evidence jsonb NOT NULL, installed_at timestamptz NOT NULL DEFAULT now()
);
