"""Prepare uncapped source diagnostics; never modify the installed anatomy."""
from pathlib import Path
import json, sys, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json, digest_file

t = time.monotonic()
base = ROOT / '.local-models/geometry-workbench/circulatory-completion'
work = base / 'carotid-bifurcation-completion-v1'
out = work / 'source-render-review-v1'
out.mkdir(exist_ok=False)
old = base / 'carotid-bifurcation-sinus-v1c'
diag = json.loads((old / 'source-diagnostic.json').read_text())
fork = np.array(next(s for s in diag['sources'] if s['side']=='l' and s['role']=='ICA')['calibre_profiles'][0]['centre_m'])
meshes = []
for role, color in [('CCA',[.8,.36,.12,1]),('ICA',[.16,.45,.72,.8]),('ECA',[.86,.16,.13,.8])]:
    m = Mesh.load(old / ('l-'+role+'.npz'))
    tri = m.triangles
    selected = (tri[:,:,1].max(1) >= fork[1]-.012) & (tri[:,:,1].min(1) <= fork[1]+.025)
    selected &= np.all(tri.max(1) >= fork-np.array([.018,.012,.018]),axis=1)
    selected &= np.all(tri.min(1) <= fork+np.array([.018,.025,.018]),axis=1)
    m = Mesh(m.vertices,m.faces[selected],m.frame,m.attributes,
             dict(m.source,diagnostic_uncapped_crop=True,source_catalogue_sha256='668c0303671b46bb6c7b9d86e6b7461181d3f815b77cc031bb7cefd255642ac3'))
    p = out / (role+'-uncapped-diagnostic.npz')
    m.save(p)
    meshes.append(dict(mesh=str(p),name=role+' current source; uncapped diagnostic crop',color=color))
body = work / 'body-placement-evidence/anterior-body-v1/carotid-body-l-candidate.npz'
meshes.append(dict(mesh=str(body),name='Uninstalled body candidate; crotch not established',color=[.55,.29,.14,1]))
spec = dict(frame='male-native-world-metres',meshes=meshes,resolution=1100,wall_seconds=40,
            views=[dict(name='current-fork-oblique',direction=[1,.2,-1]),
                   dict(name='current-nested-root-axial-cut',direction=[0,1,0],
                        section=dict(normal=[0,1,0],point=(fork+[0,.004,0]).tolist(),clear_outer=True))])
atomic_json(out/'blender-review-spec.json',spec)
atomic_json(out/'preparation.json',dict(source_diagnostic_sha256=digest_file(old/'source-diagnostic.json'),
    source_master_sha256=digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'),
    geometry_installed=False,uncapped_diagnostic=True,elapsed_seconds=time.monotonic()-t))
print(out/'blender-review-spec.json')
