"""Stage 1/2: exact terminology and fixed licensed/native source inventory."""
from pathlib import Path
import json,sqlite3,sys,time,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_continuity import boundary_loops
args=argparse.ArgumentParser();args.add_argument('--resume',action='store_true');args=args.parse_args()
T=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-retinal-vein-source-feasibility-v2';D.mkdir(exist_ok=args.resume);native_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';donor_path=ROOT/'public/models/atlas-hra-female.json';native_sha=digest_file(native_path);donor_sha=digest_file(donor_path);native=json.loads(native_path.read_text());donor=json.loads(donor_path.read_text());assert native_sha=='5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a';npieces={p['id']:p for p in native['parts']};dpieces={p['id']:p for p in donor['parts']};reader=AtlasReader()
def save_bound(mesh,path):
    if path.exists():
        assert args.resume;old=Mesh.load(path);assert old.geometry_sha==mesh.geometry_sha and old.attributes_sha==mesh.attributes_sha
    else:mesh.save(path)
authority=ROOT/'.local-models/ta98-authority.sqlite'
with sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True)as db:
    db.row_factory=sqlite3.Row;terms=[dict(r)for r in db.execute('WITH RECURSIVE d(id)AS(SELECT ? UNION ALL SELECT t.id FROM ta98 t JOIN d ON t.parent_id=d.id)SELECT t.* FROM ta98 t JOIN d USING(id)ORDER BY t.id',('A12.3.06.111',))];assert len(terms)==10 and {'A12.3.06.111','A12.3.06.112','A12.3.06.113'}<={r['id']for r in terms};parent=dict(db.execute('SELECT * FROM ta98 WHERE id=?',(terms[0]['parent_id'],)).fetchone())
report=(ROOT/'TA98-CIRCULATION.md').read_text();assert all(r['id']in report for r in terms if r['id'].startswith('A12.'))
rows=[]
for side,letter in [('left','L'),('right','R')]:
    short=side[0];roles={
      'donor-vein':'HRA:VH_F_central_retinal_vein_'+letter,
      'donor-retina':'HRA:VH_F_retina_'+letter,
      'donor-sclera':'HRA:VH_F_sclera_'+letter,
      'donor-optic-nerve':'HRA:VH_F_'+side+'_optic_nerve',
      'donor-receiving-vein':'HRA:VH_F_superior_ophthalmic_vein_'+letter,
      'donor-retinal-artery':'HRA:VH_F_central_retinal_artery_'+letter,
      'native-retina':'ZA:Retina.'+short,
      'native-sclera':'ZA:Sclera.'+short,
      'native-optic-nerve':'ZA:Optic nerve (II).'+short,
      'native-receiving-vein':'ZA:Superior ophthalmic vein.'+short,
      'native-cavernous-sinus':'ZA:Cavernous sinus.'+short,
      'native-retinal-artery':'ZA:Central retinal artery.'+short,
    }
    for role,pid in roles.items():
        is_native=role.startswith('native');a=native if is_native else donor;source=native_path if is_native else donor_path;source_sha=native_sha if is_native else donor_sha;part=(npieces if is_native else dpieces)[pid];assert not part.get('suppressed');frame='model:local-ta98:assembled-metres'if is_native else'source:hra-female:catalogue-metres';m=reader.read_part(a,part,frame=frame,catalogue_path=source,catalogue_sha=source_sha);path=D/(side+'-'+role+'.npz');save_bound(m,path);row=dict(side=side,role=role,part_id=pid,source_record=part,frame=frame,mesh=str(path),mesh_sha256=digest_file(path),geometry_sha256=m.geometry_sha,bank_path=m.source['bank_path'],bank_sha256=m.source['bank_sha256'],bounds_m=m.bounds.tolist(),triangles=len(m.faces),fixed_context=is_native,license='CC BY-SA 4.0'if is_native else'CC BY 4.0')
        if role!='native-retinal-artery':
            points,index=np.unique(m.vertices,axis=0,return_inverse=True);w=Mesh(points,index[m.faces],frame);t=topology(w);boundary_failure=None
            try:loops=boundary_loops(w.faces)
            except ValueError as error:loops=[];boundary_failure=str(error)
            wp=D/(side+'-'+role+'-exact-topology.npz');save_bound(w,wp);row.update(topology=t,exact_coordinate_welded_analysis=str(wp),boundary_polygons_verified=boundary_failure is None,boundary_extraction_failure=boundary_failure,boundary_polygons=[dict(indices=loop.tolist(),center_m=w.vertices[loop].mean(0).tolist(),vertices=len(loop))for loop in loops])
        rows.append(row)
inventory=dict(schema='central-retinal-vein-source-and-interface-inventory/v1',whole_term='A12.3.06.111',terms=[dict(authority=r,representation_type='whole-vessel finite wall and lumen'if r['id'].endswith('111')else('source-derived regional portion of the same finite vessel; boundary unresolved'if r['id'].startswith('A12.')else'named small retinal venous vessel; distinct source geometry missing'),direct_source_meshes=[x['part_id']for x in rows if x['role']=='donor-vein']if r['id'].startswith('A12.')else[],distinct_geometry_verified=False,opening_or_generic_marker=False,installed=False)for r in terms],cross_chapter_retinal_venules_not_in_current_circulation_report=[r['id']for r in terms if r['id']not in report],authority_path=str(authority),authority_sha256=digest_file(authority),checked_in_report='TA98-CIRCULATION.md',checked_in_report_sha256=digest_file(ROOT/'TA98-CIRCULATION.md'),authoritative_receiving_parent=parent,donor_catalogue=str(donor_path),donor_catalogue_sha256=donor_sha,native_catalogue=str(native_path),native_catalogue_sha256=native_sha,sources=rows,assembly_contract=dict(whole='Bilateral central retinal venous family, with independently source-witnessed intraocular/extraocular division if available',fixed=['all existing male/female anatomical records','heart','eyeballs and retina','optic nerves','skeleton','all existing vessel courses'],planned_tissues='Licensed donor outer vessel envelopes only; finite inner walls and actual paired mouths are not yet constructed',interfaces=[dict(kind='retinal/optic-nerve entry',source_witness_required=True,accepted=False),dict(kind='receiving superior ophthalmic vein',represented_native_hosts=['ZA:Superior ophthalmic vein.l','ZA:Superior ophthalmic vein.r'],actual_outflow_mouth_required=True,accepted=False)],stop_if_no_defensible_interface=True),anatomical_acceptance=False,installed=False,geometry_or_source_catalogues_modified=False,elapsed_seconds=time.monotonic()-T)
atomic_json(D/'inputs.json',inventory);assert digest_file(native_path)==native_sha and digest_file(donor_path)==donor_sha;print(json.dumps(dict(output=str(D),terms=[r['id']for r in terms],parent=parent['id'],source_meshes=len(rows),donor_veins=[dict(side=r['side'],triangles=r['triangles'],components=len(r['topology']['components']),boundary_polygons=len(r['boundary_polygons']))for r in rows if r['role']=='donor-vein'],native_hosts=[dict(side=r['side'],boundary_polygons=len(r['boundary_polygons']),open_edges=r['topology']['open_edges'])for r in rows if r['role']=='native-receiving-vein'],elapsed_seconds=inventory['elapsed_seconds'])))
