"""Qualified actual radial-envelope curves, not a fitted or copied annulus.

The preceding strictly monotone peripheral source-edge graph had no path.
This exports the actual continuous-angle source witnesses instead, retaining
any source branch jumps as an explicit obstacle rather than bridging them.
"""
from pathlib import Path
import sys,json,time,ast
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.distance import barycentric

T=time.monotonic();D=ROOT/'.local-models/geometry-workbench/heart/male-aortic-peripheral-envelope-curves-v1';D.mkdir(exist_ok=False)
prior_dir=ROOT/'.local-models/geometry-workbench/heart/male-aortic-provisional-source-hinges-v1';prior=json.loads((prior_dir/'evidence.json').read_text());center=np.array(prior['root_center_m']);axis=np.array(prior['root_axis']);u=np.cross(axis,[1,0,0]);u/=np.linalg.norm(u);v=np.cross(axis,u);meshes=[Mesh.load(prior_dir/f'cusp-{i}.npz')for i in range(3)]
# Reuse the exact meridional witness and diagnostic display helpers, without
# importing a stage or running the failed discrete source-edge solve again.
helper=ROOT/'server/propose-male-aortic-peripheral-hinge-arcs.py';parsed=ast.parse(helper.read_text());functions=[node for node in parsed.body if isinstance(node,ast.FunctionDef) and node.name in ['polar','radial_envelope','display']];exec(compile(ast.Module(body=functions,type_ignores=[]),str(helper),'exec'))
records=[];entries=[];colors=[[.15,.76,.35,1],[.2,.48,.93,1],[.94,.68,.2,1]]
for i,mesh in enumerate(meshes):
    start,end=[r for r in prior['commissures']if i in r['cusps']];aa=polar(np.array(start['point_m']))[0];bb=polar(np.array(end['point_m']))[0];middle=polar(mesh.vertices.mean(0))[0];span=float((bb-aa)%(2*np.pi))
    if float((middle-aa)%(2*np.pi))>span:start,end=end,start;aa,bb=bb,aa;span=float((bb-aa)%(2*np.pi))
    points=[];faces=[];bary=[];angles=np.linspace(aa,aa+span,257);envelope_radii=[]
    for angle in angles:
        radius,point,face,weights=radial_envelope(mesh,float(angle));points.append(point);faces.append(face);bary.append(weights);envelope_radii.append(radius)
    points=np.array(points);faces=np.array(faces);bary=np.array(bary);residual=np.linalg.norm(np.einsum('ni,nij->nj',bary,mesh.triangles[faces])-points,axis=1);assert residual.max()<1e-12
    lengths=np.linalg.norm(np.diff(points,axis=0),axis=1);height=(points-center)@axis;envelope_to_contacts=[float(np.linalg.norm(points[0]-np.array(start['point_m']))),float(np.linalg.norm(points[-1]-np.array(end['point_m'])))];largest=np.argsort(lengths)[-5:][::-1]
    witnesses=D/f'cusp-{i}-envelope.npz';atomic_npz(witnesses,points=points,source_faces=faces,source_barycentric=bary,angles=angles,
        envelope_radii=np.array(envelope_radii),adjacent_point_distance=lengths,frame=np.array(mesh.frame))
    curve=D/f'cusp-{i}-diagnostic-envelope.npz';display(points,curve);mp=D/f'cusp-{i}.npz';mesh.save(mp);entries.extend([dict(mesh=str(mp),name=prior['arcs'][i]['cusp']+' unchanged source',color=colors[i]),dict(mesh=str(curve),name='Qualified actual peripheral envelope — not a certified hinge',color=[1,.12,.55,1])])
    records.append(dict(cusp=prior['arcs'][i]['cusp'],source_geometry_sha256=mesh.geometry_sha,witness_file=str(witnesses),controls=len(points),sector_degrees=float(span*180/np.pi),
        minimum_radius_mm=float(min(envelope_radii)*1000),minimum_source_envelope_fraction=1.,angular_backtracking_segments=0,
        maximum_barycentric_residual_m=float(residual.max()),height_range_mm=[float(height.min()*1000),float(height.max()*1000)],
        source_commissure_endpoint_gaps_mm=[x*1000 for x in envelope_to_contacts],maximum_adjacent_sample_distance_mm=float(lengths.max()*1000),
        largest_step_witnesses=[dict(first_control=int(k),first_source_face=int(faces[k]),second_source_face=int(faces[k+1]),distance_mm=float(lengths[k]*1000),height_jump_mm=float(abs(height[k+1]-height[k])*1000))for k in largest],
        classification='Actual sector-monotone outer envelope; anatomical hinge status and any branch/endpoint discontinuities remain explicit'))
atomic_json(D/'review-spec.json',dict(frame=meshes[0].frame,resolution=1100,wall_seconds=50,meshes=entries,views=[dict(name='Actual peripheral radial envelope — oblique',direction=[.4,.7,1]),dict(name='Outward envelope differs from pocket floor — side',direction=[-1,.1,.4]),dict(name='Peripheral sector traces — below',direction=(-axis).tolist())]))
report=dict(schema='male-aortic-peripheral-envelope-curves/v1',source_catalogue_sha256=prior['source_catalogue_sha256'],source_geometry_modified=False,
    qualified_curve_only=True,root_tissue_constructed=False,convex_hull_used=False,source_reference_copied=False,
    explicit_preceding_obstacle='No strictly sector-monotone original source-edge path at >=80% actual radial envelope; no gate was relaxed.',
    arcs=records,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-T,
    references=[dict(url='https://www.ahajournals.org/doi/pdf/10.1161/CIRCINTERVENTIONS.108.780858',scope='Parent-reviewed anatomical distinction: semilunar crown attachment versus virtual basal ring',access='403 in this subagent fetch'),dict(url='https://pubmed.ncbi.nlm.nih.gov/28160343/',scope='Parent-provided CT orientation evidence, not copied coordinates')])
atomic_json(D/'evidence.json',report);print(json.dumps(report))
