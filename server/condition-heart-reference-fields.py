"""Bounded C++ signed-distance reconstruction of rejected reference buffers.

Only newly constructed reference geometry is sampled. This never processes
the fixed anatomical bank or copies triangulation from a private donor.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
import open3d as o3d
import manifold3d as md
from scipy.ndimage import gaussian_filter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json
from geometry_tools.algorithms.regular_grid_isosurface import extract
from geometry_tools.algorithms.faithful_float32 import oriented_rounding
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--only',nargs='*');args=ap.parse_args();base=ROOT/'.local-models/geometry-workbench/heart';src=base/'typed-audited-v659';out=base/args.stage;out.mkdir(exist_ok=False);started=time.monotonic();meta=copy.deepcopy(json.loads((src/'construction.json').read_text()));previous=json.loads((src/'float32-audit.json').read_text())['geometries'];rows={}
for term,g in meta['geometries'].items():
    if args.only and term not in args.only:continue
    m=Mesh.load(src/g['file']);evidence=previous[term]
    if not evidence['passed']:
        center=m.vertices.mean(0);axes=np.linalg.svd(m.vertices-center,full_matrices=False)[2];v=(m.vertices-center)@axes.T;span=np.ptp(v,axis=0);pitch=max(float(span.max()/200),float(np.cbrt(np.prod(span)/1200000)),.00006);lo=v.min(0)-2.314159*pitch;shape=np.ceil((v.max(0)+2.27*pitch-lo)/pitch).astype(int)+1
        assert int(np.prod(shape))<2000000,'Explicit bounded sampling budget exceeded'
        grid=np.indices(tuple(shape)).reshape(3,-1).T*pitch+lo;scene=o3d.t.geometry.RaycastingScene();scene.add_triangles(o3d.core.Tensor(v.astype('float32')),o3d.core.Tensor(m.faces.astype('uint32')));field=scene.compute_signed_distance(o3d.core.Tensor(grid.astype('float32')),nthreads=1,nsamples=3).numpy().reshape(tuple(shape));field=gaussian_filter(field,.35)-pitch*.017
        vv,ff,contour=extract(field,lo,pitch);vv=vv@axes+center
        signed=np.einsum('ij,ij->i',vv[ff[:,0]],np.cross(vv[ff[:,1]],vv[ff[:,2]])).sum()
        if signed<0:ff=ff[:,::-1]
        s=md.Manifold(md.Mesh64(np.array(vv,order='C'),np.array(ff,dtype='uint64',order='C'))).as_original().simplify(pitch*.025)
        assert str(s.status())=='Error.NoError'and not s.is_empty()
        b=s.to_mesh64();source=np.asarray(b.vert_properties[:,:3]);faces=np.asarray(b.tri_verts,dtype=int);n=np.cross(source[faces[:,1]]-source[faces[:,0]],source[faces[:,2]]-source[faces[:,0]]);encoded,rounding=oriented_rounding(source,source.astype('float32').astype(float),faces,np.repeat(n[:,None,:],3,axis=1),maximum_states=262144,time_limit=5.);m=Mesh(encoded,faces,'heart-native-local-metres');evidence=audit(m,require_closed=True,max_pairs=1000000,max_tests=4000000);g['faithfulBufferRounding']=rounding
        g['fieldConditioning']=dict(samplePitchMetres=pitch,samples=int(np.prod(shape)),source_reference_stage=src.name,normalSmoothingSigmaVoxels=.35,levelShiftMetres=pitch*.017,originalAnatomyChanged=False,accuracy='bounded-resolution anatomical reference; exact histological boundary remains provisional')
    m.save(out/g['file']);g.update(vertices=len(m.vertices),triangles=len(m.faces));rows[term]=evidence;atomic_json(out/'float32-audit.json',dict(complete=False,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True);print(term,evidence['passed'],evidence['failures'],len(m.faces),'triangles',round(time.monotonic()-started,2),'s',flush=True)
atomic_json(out/'construction.json',meta);atomic_json(out/'float32-audit.json',dict(complete=True,geometries=rows,elapsed_seconds=time.monotonic()-started),overwrite=True);failed=[t for t,e in rows.items()if not e['passed']];print('Complete',len(rows),'failed',failed,flush=True)
if failed:raise SystemExit(2)
