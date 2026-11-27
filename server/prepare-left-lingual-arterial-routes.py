"""Source-bound first left lingual family, recipient and calibre contract."""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
import igl
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_sections import centripetal_curve
from geometry_tools.algorithms.tube_continuity import finite_wall_cells
from geometry_tools.audit import audit

args=argparse.ArgumentParser();args.add_argument('--version',default='v3');args=args.parse_args()
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';P=B/'male-lingual-source-plan-v1';D=B/('male-left-lingual-routes-'+args.version);D.mkdir(exist_ok=False);inputs=json.loads((P/'inputs.json').read_text());records={r['role']:r for r in inputs['fixed_sources']};meshes={role:Mesh.load(ROOT/row['mesh'])for role,row in records.items()};frame=inputs['frame']
for role,row in records.items():assert digest_file(ROOT/row['mesh'])==row['mesh_sha256'],role
eca=meshes['external-carotid-left'];station=records['external-carotid-left']['source_station_components'][0];rings=np.array(station['rings']);centers=np.array(station['centres_m']);tip=meshes['greater-horn-left'].vertices[np.argmax(meshes['greater-horn-left'].vertices[:,1])];origin_station=int(np.argmin(abs(centers[:,1]-tip[1])));origin=centers[origin_station];native_radius=float(np.linalg.norm(eca.vertices[rings[origin_station]]-origin,axis=1).mean())
whole=Mesh.load(ROOT/'.local-models/geometry-workbench/tongue-completion/basal-whole-protected-v1/whole-tongue.npz')
inner=Mesh.load(ROOT/'.local-models/geometry-workbench/tongue-completion/intrinsic-territories-v2/native-inner-envelope.npz')
whole.save(D/'whole-tongue-fixed-reference.npz');inner.save(D/'inner-tongue-fixed-reference.npz')
def recipient(mesh,xyz,radius,extra=None):
    low,high=mesh.bounds;seed=np.array(xyz);x=np.linspace(max(low[0],seed[0]-.003),min(high[0],seed[0]+.003),9);y=np.linspace(max(low[1],seed[1]-.003),min(high[1],seed[1]+.003),9);z=np.linspace(max(low[2],seed[2]-.002),min(high[2],seed[2]+.002),7);points=np.stack(np.meshgrid(x,y,z,indexing='ij'),-1).reshape(-1,3);inside=abs(igl.winding_number(mesh.vertices,mesh.faces,points))>.99;distance=np.sqrt(igl.point_mesh_squared_distance(points,mesh.vertices,mesh.faces)[0]);valid=inside&(distance>radius+.0001)
    if extra is not None:
        valid&=abs(igl.winding_number(extra.vertices,extra.faces,points))>.99
        valid&=np.sqrt(igl.point_mesh_squared_distance(points,extra.vertices,extra.faces)[0])>radius+.0002
    assert valid.any(),('No source recipient clearance',xyz,radius,float(distance.max()))
    candidates=np.flatnonzero(valid);k=candidates[np.argmin(np.linalg.norm(points[candidates]-seed,axis=1))]
    return points[k],dict(point_m=points[k].tolist(),inside_closed_source=True,source_geometry_sha256=mesh.geometry_sha,minimum_control_to_recipient_surface_m=float(distance[k]),branch_outer_radius_m=radius)

suprahyoid_target,suprahyoid_witness=recipient(meshes['hyoglossus-left'],[.0105,1.6165,.039],.00060)
dorsal0,dorsal0_witness=recipient(inner,[.0085,1.635,.0335],.00065)
dorsal1,dorsal1_witness=recipient(inner,[.009,1.638,.043],.00060)
sublingual,sublingual_witness=recipient(meshes['sublingual-gland-left'],[.012,1.620,.070],.00085)
deep,deep_witness=recipient(inner,[.0085,1.638,.075],.00085)
deep_witness['recipient_scope']='Actual closed inner tongue parenchymal envelope. The pre-existing inferior-longitudinal ROI is explicitly approximate and did not provide a common cleared recipient at the first proposed terminal point.'
hg=meshes['hyoglossus-left'];medial_x=float(hg.bounds[0,0]-.0018-.00015);anterior_z=float(hg.bounds[1,2]+.0015)
waypoints=[origin,[origin[0],1.620,.016],[.025,1.621,.024],[.015,1.622,.029],[medial_x,1.6225,.036],[medial_x,1.624,.049],[medial_x,1.626,.061],[medial_x,1.628,anterior_z]]
main,main_tangent=centripetal_curve(waypoints,chord_tolerance_m=.00001,max_step_m=.0007);main_radius=np.linspace(.0018,.0012,len(main))
def parent_site(point):
    j=int(np.argmin(np.linalg.norm(main-np.array(point),axis=1)));return main[j],j
sup_start,sup_index=parent_site([.023,1.621,.019]);d0_start,d0_index=parent_site([medial_x,1.6225,.037]);d1_start,d1_index=parent_site([medial_x,1.6235,.046])
sub_start,sub_index=parent_site([medial_x,1.6255,.060])
routes=[dict(name='lingual',term='A12.2.05.015',path=main,radii=main_radius,wall=.00022,parent='external-carotid-left',parent_station=origin_station,recipient=None),
dict(name='suprahyoid',term='A12.2.05.016',waypoints=[sup_start,[.020,1.6195,.026],suprahyoid_target],radius=.00060,parent='lingual',parent_station=sup_index,recipient=suprahyoid_witness),
dict(name='dorsal-lingual-0',term='A12.2.05.017',waypoints=[d0_start,[.0045,1.629,.034],dorsal0],radius=.00065,parent='lingual',parent_station=d0_index,recipient=dorsal0_witness),
dict(name='dorsal-lingual-1',term='A12.2.05.017',waypoints=[d1_start,[.005,1.632,.043],dorsal1],radius=.00060,parent='lingual',parent_station=d1_index,recipient=dorsal1_witness),
dict(name='sublingual',term='A12.2.05.018',waypoints=[sub_start,[.007,1.6215,.063],sublingual],radius=.00085,parent='lingual',parent_station=sub_index,recipient=sublingual_witness),
dict(name='deep-lingual',term='A12.2.05.019',waypoints=[main[-1],[.007,1.631,.074],deep],radius=.00085,parent='lingual',parent_station=len(main)-1,recipient=deep_witness)]
report=[]
for r in routes:
    if 'waypoints'in r:r['path'],_=centripetal_curve(r.pop('waypoints'),chord_tolerance_m=.00001,max_step_m=.0005);r['radii']=np.full(len(r['path']),r.pop('radius'));r['wall']=.00012
    file=D/(r['name']+'-guide.npz');atomic_npz(file,points=r['path'],outer_radii=r['radii'],inner_radii=r['radii']-r['wall'],frame=np.array(frame));report.append(dict(name=r['name'],term=r['term'],guide=str(file),guide_sha256=digest_file(file),stations=len(r['path']),length_m=float(np.linalg.norm(np.diff(r['path'],axis=0),axis=1).sum()),radius_m=[float(r['radii'].min()),float(r['radii'].max())],illustrative_wall_m=r['wall'],parent=r['parent'],parent_station=r['parent_station'],recipient=r['recipient']))
plan=dict(schema='left-lingual-source-bound-route-contract/v1',source_inputs_sha256=digest_file(P/'inputs.json'),source_catalogue_sha256=inputs['source_catalogue_sha256'],
    fixed_native_carotid=records['external-carotid-left'],source_origin_station=origin_station,native_origin_outer_radius_m=native_radius,hyoid_tip_source_point_m=tip.tolist(),
    planned_assembly='One left external-carotid communicating parent and six branch tubes representing five authoritative lingual-family terms',
    planned_parts=inputs['terms'],interfaces=[dict(parent=r['parent'],child=r['name'],blood_opening_required=True,verified=False)for r in report],
    terminal_interfaces='Five open blood terminal mouths embedded in recipient tissues; unresolved microscopic perfusion beds are not fabricated',
    main_medial_hyoglossus_clearance='At the protected medial course center_x + maximum_outer_radius <= global native HG medial extent - 0.15 mm; full emitted surfaces still require checks',
    references=inputs['anatomical_references']+[dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC8899355/',role='Original 70-side cadaveric diameter study: solitary LA mean diameter 3.6 mm; this prototype uses radius 1.8 mm proximally, not a copied mesh')],
    caliber_scope='Measured native ECA profile and a primary anatomical 3.6 mm lingual diameter; taper and small-branch calibres are illustrative',
    routes=report,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-T)
atomic_json(D/'routes.json',plan);print(json.dumps(dict(elapsed_seconds=plan['elapsed_seconds'],origin_station=origin_station,native_parent_radius_mm=native_radius*1000,routes=[dict(name=r['name'],length_mm=r['length_m']*1000,stations=r['stations'],recipient=r['recipient'])for r in report])))
