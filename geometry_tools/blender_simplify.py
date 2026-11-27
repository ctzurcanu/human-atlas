"""Reduce planar chart triangulation in the already open Blender process.

Uses only new, unlinked BMeshes. Existing objects, coordinates, scene, camera
and the user's main file are untouched. Boundary controls stay literal, and
the output needs the normal external encoded geometry/review gates.
"""
from pathlib import Path
import hashlib,json,time
import bpy,bmesh,numpy as np
from collections import Counter

ROOT=Path('/Users/christiantzurcanu/Documents/dev/human-atlas')
REQUEST=Path(globals().get('renal_simplify_request',ROOT/'.local-models/geometry-workbench/renal-completion/simplify-v438/request.json'))


def run():
    started=time.monotonic();r=json.loads(REQUEST.read_text());source=Path(r['source']);output=Path(r['output'])
    if not source.resolve().is_relative_to(ROOT/'.local-models')or not output.resolve().is_relative_to(ROOT/'.local-models'):raise ValueError('Only local model files are allowed')
    if hashlib.sha256(source.read_bytes()).hexdigest()!=r['source_sha256']or output.exists():raise ValueError('Changed source or existing output')
    before={o.name:(o.as_pointer(),tuple(x for row in o.matrix_world for x in row),o.data.as_pointer()if o.data else None)for o in bpy.data.objects}
    z=np.load(source);positions=z['positions'];origin=np.array(r['origin']);arrays={'positions':positions};rows=[]
    def boundary(f):
        e=Counter(tuple(sorted((int(a),int(b))))for tri in f for a,b in zip(tri,np.roll(tri,-1)))
        return {q for q,n in e.items()if n==1}
    for key in z.files:
        if not key.startswith('surface_'):continue
        f=z[key];used,inv=np.unique(f,return_inverse=True);bm=bmesh.new();layer=bm.verts.layers.int.new('source_control')
        for i in used:
            v=bm.verts.new(positions[i]-origin);v[layer]=int(i)
        bm.verts.ensure_lookup_table()
        for tri in inv.reshape(-1,3):bm.faces.new([bm.verts[int(i)]for i in tri])
        bm.normal_update();edges=[e for e in bm.edges if len(e.link_faces)==2]
        bmesh.ops.dissolve_limit(bm,angle_limit=r['angle_rad'],use_dissolve_boundaries=False,verts=[],edges=edges,delimit=set())
        bmesh.ops.triangulate(bm,faces=list(bm.faces),quad_method='BEAUTY',ngon_method='EAR_CLIP')
        candidate=np.array([[v[layer]for v in face.verts]for face in bm.faces],int);bm.free()
        kept=boundary(candidate)==boundary(f)
        arrays[key]=candidate if kept else f
        rows.append(dict(surface=key,before_triangles=len(f),after_triangles=len(arrays[key]),boundary_controls_literal=kept))
    after={o.name:(o.as_pointer(),tuple(x for row in o.matrix_world for x in row),o.data.as_pointer()if o.data else None)for o in bpy.data.objects}
    if before!=after:raise ValueError('Existing Blender object state changed')
    output.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(output,**arrays)
    result=dict(source=str(source),source_sha256=r['source_sha256'],output=str(output),output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                rows=rows,positions_literal=True,existing_objects_literal=True,main_file_saved=False,seconds=time.monotonic()-started,requires_external_validation=True)
    output.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Shared renal charts reduced; original Blender scene retained; external validation required.',flush=True)


run()
