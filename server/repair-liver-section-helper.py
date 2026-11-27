"""Repair only a render classifier; leave all anatomical tissue literal."""
from pathlib import Path
import argparse, copy, json, sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_bytes, atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import audit
from geometry_tools.hepatic_course_field import HepaticCourseField
from geometry_tools.algorithms.adaptive_course_transport import refine_course
from geometry_tools.algorithms.transport_slivers import regularize_transport_slivers

ap = argparse.ArgumentParser()
ap.add_argument('--stage', choices=['construct', 'audit'], required=True)
args = ap.parse_args()
base = ROOT / '.local-models/geometry-workbench/liver-completion'
out = base / 'completion-v200'
reader = AtlasReader()
if args.stage == 'construct':
    out.mkdir(exist_ok=False)
    original = json.loads((base / 'completion-v155/atlas-ready.json').read_text())
    a = json.loads((base / 'completion-v186/atlas.json').read_text())
    source = original['liverSectionIntegration']
    row = next(r for r in source['exclusions'] if r['id'] == 'biliary')
    raw = reader.asset(source['bankUrl']).read_bytes()
    origin = np.array(row['origin'])
    v = np.ndarray((row['vertexCount'], 3), '<f4', buffer=raw, offset=row['positions']).astype(float) + origin
    f = np.ndarray((row['indexCount'] // 3, 3), '<u4', buffer=raw, offset=row['indices']).copy()
    field = HepaticCourseField(original, base / 'hepatic-contained-course-v163/field.npz', base / 'separated-courses-v105/evidence.json')
    normals = np.tile([1., 0., 0.], (len(v), 1))
    _, _, _, parent, _ = refine_course(v, f, normals, field, deviation_m=1e-5, max_rounds=18)
    bad = json.loads((base / 'completion-v186/exclusion-biliary-audit.json').read_text())
    fine = np.unique(parent[np.unique(np.asarray(bad['self_intersections']['examples'], dtype=int))]).tolist()
    v, f, _, _, refinement = refine_course(v, f, normals, field, deviation_m=1e-5, max_rounds=18, fine_parents=fine)
    mapped = field(v)
    local = (mapped - origin).astype('<f4')
    f, slivers = regularize_transport_slivers(v, local.astype(float) + origin, f, np.all(mapped == v, axis=1))
    mesh = Mesh(local.astype(float) + origin, f, 'model:local-ta98:assembled-metres')
    mesh.save(out / 'exclusion-biliary.npz')
    atomic_json(out / 'exclusion-biliary-construction.json', dict(refinement=refinement, slivers=slivers, anatomyLiteral=True))
    print(json.dumps(dict(triangles=len(f), fineParents=fine)))
    sys.exit(0)

mesh = Mesh.load(out / 'exclusion-biliary.npz')
check = audit(mesh, require_closed=True, max_pairs=3000000, max_tests=24000000)
atomic_json(out / 'exclusion-biliary-audit.json', check)
assert check['passed'], check['failures']
a = json.loads((base / 'completion-v186/atlas.json').read_text())
for role in ['artery', 'portal', 'hepatic_vein']:
    assert json.loads((base / f'completion-v186/exclusion-{role}-audit.json').read_text())['passed']
origin = np.array([0., 1.25, 0.])
# The original classifier origin is part of its literal bank contract.
oldrow = next(e for e in a['liverSectionIntegration']['exclusions'] if e['id'] == 'biliary')
origin = np.array(oldrow['origin'])
local = (mesh.vertices - origin).astype('<f4')
def append(data, value):
    data.extend(b'\0' * ((-len(data)) % 4))
    offset = len(data)
    data.extend(value)
    return offset
def bank(name, marker):
    data = bytearray(reader.asset(marker['bankUrl']).read_bytes())
    row = copy.deepcopy(oldrow)
    row.update(vertexCount=len(local), indexCount=mesh.faces.size,
               positions=append(data, local.tobytes()),
               indices=append(data, mesh.faces.astype('<u4').tobytes()),
               geometrySha256=mesh.geometry_sha)
    path = out / name
    atomic_bytes(path, bytes(data))
    url = '/local-models/' + str(path.relative_to(ROOT / '.local-models'))
    ci = len(a['chunks'])
    a['chunks'].append(dict(url=url, bytes=len(data)))
    marker.update(bankUrl=url, bankSha256=digest_file(path),
                  exclusions=[row if e['id'] == 'biliary' else e for e in marker['exclusions']])
    return ci, url, digest_file(path)
ci, url, sha = bank('liver-section-bank.bin', a['liverSectionIntegration'])
parts = {p['id']: p for p in a['parts']}
for pid in a['liverIntegration']['partIds']:
    parts[pid]['chunk'] = ci
a['liverIntegration'].update(bankUrl=url, bankSha256=sha)
a['liverFeatureIntegration']['bankSha256'] = sha
for group in a['liverFeatureIntegration']['termGroups']:
    group['surface']['sourceBankSha256'] = sha
a['liverCompletionIntegration']['hostBankSha256'] = sha
functional = a['liverCompletionIntegration']['functionalSections']
fci, _, _ = bank('functional-section-bank.bin', functional)
functional['sourceBankSha256'] = sha
for pid in a['liverCompletionIntegration']['functionalPartIds']:
    parts[pid]['chunk'] = fci
atomic_json(out / 'atlas.json', a)
atomic_json(out / 'construction.json', dict(sourceCandidate=str(base / 'completion-v186/atlas.json'),
    anatomyLiteral=True, sectionClassifierOnly=True, passed=True,
    biliaryClassifierAudit=str(out / 'exclusion-biliary-audit.json')))
print(json.dumps(dict(candidate=str(out / 'atlas.json'), passed=True, triangles=len(mesh.faces))))
