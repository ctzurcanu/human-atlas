"""Sample every internal calculation edge, not just tetrahedral centres."""
from pathlib import Path
import sys,time,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_npz,atomic_json,digest_file
from geometry_tools.algorithms.tetra_edge_refinement import refine_tetrahedral_edges
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'wall-calculation-v26/tetrahedra.npz';out=b/'wall-calculation-v50';out.mkdir(exist_ok=False);z=np.load(source);r=refine_tetrahedral_edges(z['positions'],z['tetrahedra']);p=r['positions'];t=r['tetrahedra'];q=p[t];volume=np.linalg.det(np.transpose(q[:,1:]-q[:,:1],(0,2,1))).sum()/6;error=abs(volume-314.07135495017184e-6);assert error<1e-12
atomic_npz(out/'tetrahedra.npz',positions=p,tetrahedra=t,source_faces=z['source_faces'],source_vertices=z['source_vertices'],original_positions=np.array(len(z['positions'])),frame=z['frame']);report=dict(passed=True,source_sha256=digest_file(source),positions=len(p),tetrahedra=len(t),shared_edge_midpoints=len(r['edges']),original_boundary_controls_literal=bool(np.array_equal(p[:len(z['positions'])],z['positions'])),volume_mL=float(volume*1e6),absolute_volume_error_m3=float(error),elapsed_seconds=time.monotonic()-start,construction_only=True,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
