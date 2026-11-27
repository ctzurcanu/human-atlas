"""Read-only native valve/root source topology; never creates an annulus."""
from pathlib import Path
import sys, json, time, struct, zipfile
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, atomic_npz, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from collections import defaultdict

started = time.monotonic()
out = ROOT / '.local-models/geometry-workbench/heart/male-aortic-alternative-native-sources-v1'
out.mkdir(exist_ok=False)
reader = AtlasReader()
records = []
for sex in ['male', 'female']:
    path = ROOT / f'.local-models/{sex}.json'
    atlas = json.loads(path.read_text())
    for part in atlas['parts']:
        if part.get('system') != 'cardiac' or part.get('sourceId') not in ['aortic_semilunar_valve', 'aortic_arch', 'ventricles']:
            continue
        mesh = reader.read_part(atlas, part, frame=f'model:{sex}-detail:assembled-metres')
        used = np.unique(mesh.faces)
        vertices, inverse = np.unique(mesh.vertices[used], axis=0, return_inverse=True)
        mapping = np.full(len(mesh.vertices), -1, dtype=int)
        mapping[used] = inverse
        faces = mapping[mesh.faces]
        edges = defaultdict(list)
        for fi, face in enumerate(faces):
            for aa, bb in zip(face, np.roll(face, -1)):
                edges[tuple(sorted([int(aa), int(bb)]))].append(fi)
        boundary = np.array([edge for edge, owners in edges.items() if len(owners) == 1], int).reshape(-1, 2)
        graph = defaultdict(set)
        for aa, bb in boundary:
            graph[int(aa)].add(int(bb)); graph[int(bb)].add(int(aa))
        pending = set(graph); components = []
        while pending:
            seed = min(pending); stack = [seed]; group = set(); pending.remove(seed)
            while stack:
                aa = stack.pop(); group.add(aa)
                for bb in graph[aa]:
                    if bb in pending:
                        pending.remove(bb); stack.append(bb)
            points = vertices[sorted(group)]
            components.append(dict(controls=len(group), degrees=sorted(set(len(graph[k]) for k in group)),
                                   bounds_mm=(points.min(0)*1000).tolist()+(points.max(0)*1000).tolist(),
                                   center_m=points.mean(0).tolist()))
        name = sex + '-' + part['sourceId']
        mesh.save(out / (name + '.npz'))
        atomic_npz(out/(name+'-boundary.npz'), vertices=vertices, faces=faces, boundary_edges=boundary,
                   frame=np.array(mesh.frame))
        records.append(dict(id=part['id'], source_catalogue=str(path), source_catalogue_sha256=digest_file(path),
                            geometry_sha256=mesh.geometry_sha, topology=topology(mesh),
                            boundary_components=components, declared_hinge=False,
                            declared_reciprocal_valve_root_interface=False))
glbs = []
for path in [ROOT/'.local-models/models/Male_2024.glb', ROOT/'.local-models/models/Female_2024.glb']:
    with path.open('rb') as stream:
        header=stream.read(12); length, kind=struct.unpack('<II',stream.read(8)); doc=json.loads(stream.read(length))
    glbs.append(dict(path=str(path), size=path.stat().st_size, json_chunk_size=length,
                     node_count=len(doc.get('nodes',[])), mesh_count=len(doc.get('meshes',[])),
                     semantic_aortic_nodes=[dict(index=i, **node) for i,node in enumerate(doc.get('nodes',[]))
                                            if any(key in node.get('name','').lower() for key in ['aort','valv','heart','ventric'])],
                     source_note='Node names are opaque hashes; no semantic aortic attachment annotation was found in the JSON node names.'))
archives=[]
for path in (ROOT/'.local-models/source').glob('*.zip'):
    with zipfile.ZipFile(path) as source:
        names=source.namelist()
    archives.append(dict(path=str(path), entries=len(names), first_names=names[:4],
                         semantic_matches=[name for name in names if any(key in name.lower() for key in ['aort','valv','heart','ventric'])]))
report=dict(schema='native-aortic-alternative-source-diagnostic/v1', read_only=True,
            anatomical_acceptance=False, annulus_constructed=False, records=records,
            glb_header_inventory=glbs, archive_inventory=archives, elapsed_seconds=time.monotonic()-started)
atomic_json(out/'evidence.json', report)
print(json.dumps(dict(elapsed_seconds=report['elapsed_seconds'], records=[dict(id=r['id'],
      topology=r['topology'], boundary_components=r['boundary_components'])for r in records])))
