"""Construction-only tetrahedra inside the complete audited myocardial wall.

Neither tetrahedra nor interior calculation faces are deployed to the viewer.
The exact original surface controls remain fixed. Volume replay rejects
accidental filling of blood chambers or vascular passages.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import meshpy.tet as tet
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';source=base/'whole-wall-v25';assert json.loads((source/'construction.json').read_text())['passed'];out=base/'wall-calculation-v26';out.mkdir(exist_ok=False);m=Mesh.load(source/'whole-wall.npz');o=np.array([0.,1.375,.03125]);info=tet.MeshInfo();info.set_points(((m.vertices-o)*1000).tolist());info.set_facets(m.faces.tolist());result=tet.build(info,options=tet.Options('pQY'));v=np.array(result.points)/1000+o;t=np.array(result.elements,dtype='int64');assert np.linalg.norm(v[:len(m.vertices)]-m.vertices,axis=1).max()<1e-12;v[:len(m.vertices)]=m.vertices
det=np.linalg.det(np.transpose(v[t[:,1:]]-v[t[:,:1]],(0,2,1)));negative=det<0;t[negative]=t[negative][:,[0,2,1,3]];assert np.all(det!=0);volume=abs(det).sum()/6;q=m.triangles-o;reference=np.einsum('ij,ij->i',q[:,0],np.cross(q[:,1],q[:,2])).sum()/6;error=abs(volume-reference);passed=error<reference*1e-8;atomic_npz(out/'tetrahedra.npz',positions=v,tetrahedra=t,source_faces=m.faces,source_vertices=np.array(len(m.vertices)),frame=np.array(m.frame));report=dict(passed=bool(passed),source_sha256=digest_file(source/'whole-wall.npz'),source_geometry_sha256=m.geometry_sha,source_surface_controls_literal=True,source_surface_faces_literal=True,positions=len(v),tetrahedra=len(t),tissue_volume_mL=float(volume*1e6),source_tissue_volume_mL=float(reference*1e6),absolute_volume_error_m3=float(error),elapsed_seconds=time.monotonic()-start,construction_only=True,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
