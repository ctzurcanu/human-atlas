"""One bounded literal left-MCFA source and hip-landmark preparation stage."""
import json, re, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion'
DEST=BASE/'nonfoot-lower-limb-blood-frontier-v2'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
CENSUS=ROOT/'reports/ta98-blood-circulation-inventory-installed-v13.json'
ORD=ROOT/'.local-models/ta98-male.json'
EXPECTED='997d37f02e7b2324ff756700c5da22918753a6daec9e3c399c2d4339b6c150b8'

def main():
    began=time.monotonic();assert not DEST.exists();DEST.mkdir()
    failed=BASE/'nonfoot-lower-limb-blood-frontier-v1'
    atomic_json(failed/'failure.json',dict(status='failed-before-output',stage='inventory-parser',error="ValueError: invalid literal for int() with base 10: '7M'",elapsed_seconds=.43,runtime_receipt=str(ROOT/'.local-models/geometry-workbench/runtime/4c09ab994f3c4b0b828c547361f06c66.json'),disposition='Retained. Deterministic scalar parser correction in unique v2; no geometry trial.'))
    sha=digest_file(CAT);assert sha==EXPECTED
    atlas=json.loads(CAT.read_text());ordinary=json.loads(ORD.read_text());census=json.loads(CENSUS.read_text());parts={p['id']:p for p in atlas['parts']};oparts={p['id']:p for p in ordinary['parts']}
    frame='model:local-ta98:assembled-metres';reader=AtlasReader(ROOT)
    def number(t):return int(re.search(r'\.(\d+)[MF]?$',t).group(1))
    scoped=[t for t in census['terms']if t['male_applicable']and((t['id'].startswith('A12.2.16.')and(number(t['id'])in set(range(2,45))|{55,56,60,71,77}))or(t['id'].startswith('A12.3.11.')and(number(t['id'])in set(range(3,11))|set(range(23,34)))))]
    family=[t for t in census['terms']if t['id']in['A12.2.16.%03d'%i for i in range(21,27)]]
    atomic_json(DEST/'inventory.json',dict(authority_cache=census['authority_cache'],authority_cache_sha256=census['authority_cache_sha256'],source_census=str(CENSUS),source_census_sha256=digest_file(CENSUS),scope='Explicit non-foot lower-limb arterial/venous frontier; male-applicable authority terms. Does not change global denominator.',scope_terms=scoped,scope_count=len(scoped),scope_missing_direct_count=sum(not t['active_exact_source_part_ids']for t in scoped),whole_family_terms=family,missing_named_family_terms=5,current_distinct_named_missing_branch_meshes=0,whole_source_parent_ids=['ZA:Medial circumflex femoral artery.l','ZA:Medial circumflex femoral artery.r'],ready_literal_partitions=0,classification='All six family entries are macroscopic arterial course/network terms, not openings, wall layers, surfaces or microscopic beds.'))
    parent=parts['ZA:Medial circumflex femoral artery.l'];bb=np.array(parent['bounds']);margin=.02
    pattern=r'Obturator|Quadratus femoris|Adductor|Pectineus|Femur\.l|Hip bone\.l|hip joint|femur\.l|femoral|Femoral|Iliopsoas|Iliacus|Psoas major|Sciatic nerve|Nerve to quadratus|Ischiofemoral|Pubofemoral|iliofemoral|Acetabular|acetabular|Frenula|Gluteus minimus|Gluteus medius|Piriformis muscle|gemellus|fascia lata|Fascia lata'
    selected=[]
    for p in atlas['parts']:
        if not p['id'].startswith('ZA:')or p.get('system')=='lymphatic'or not re.search(pattern,p['id']):continue
        bounds=np.array(p['bounds'])
        if np.all(bounds[1]>=bb[0]-margin)and np.all(bounds[0]<=bb[1]+margin):selected.append(p)
    meshes={};bindings=[];banks={};(DEST/'sources').mkdir()
    for i,p in enumerate(selected):
        assert time.monotonic()-began<45,'Declared source-preparation budget nearly exhausted'
        m=reader.read_part(atlas,p,frame=frame,catalogue_path=CAT,catalogue_sha=sha);meshes[p['id']]=m
        path=DEST/'sources'/f'source-{i:02d}.npz';m.save(path)
        url=atlas['chunks'][p['chunk']]['url'];bank=reader.asset(url)
        if str(bank)not in banks:banks[str(bank)]=digest_file(bank)
        old=None
        if p['id']in oparts:
            om=reader.read_part(ordinary,oparts[p['id']],frame=frame,catalogue_path=ORD,catalogue_sha=digest_file(ORD))
            old=dict(catalogue=str(ORD),catalogue_sha256=digest_file(ORD),geometry_sha256=om.geometry_sha,attributes_sha256=om.attributes_sha,whole_geometry_exact=m.geometry_sha==om.geometry_sha,all_attributes_exact=m.attributes_sha==om.attributes_sha)
        bindings.append(dict(part_id=p['id'],mesh=str(path),mesh_sha256=digest_file(path),geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,faces=len(m.faces),vertices=len(m.vertices),source_record=p,ordinary_body_match=old))
    m=meshes[parent['id']];v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];courses=recover_station_rings(v,f);assert len(courses)==1;q=courses[0]
    rings,proof=align_station_correspondence(v,f[q['face_ids']],q['rings']);rings=np.array(rings);profiles=v[rings];centers=profiles.mean(1)
    station=np.full(len(v),-1,int)
    for j,ring in enumerate(rings):station[ring]=j
    source_interval=np.min(station[f],axis=1);assert np.all(np.ptp(station[f],axis=1)==1)
    atomic_npz(DEST/'left-parent-stations.npz',physical_vertices=v,physical_faces=f,source_face_ids=np.arange(len(f)),ordered_station_control_ids=rings,profiles=profiles,centers=centers,source_face_station_interval=source_interval)
    witnesses={};summaries=[]
    keys=['ZA:Obturator externus.l','ZA:Quadratus femoris muscle.l','ZA:Articular capsule of hip joint.l','ZA:Femur.l:Bone','ZA:Femur.l:Cartilage','ZA:Acetabular labrum.l','ZA:Deep femoral artery.l']
    for j,pid in enumerate(keys):
        assert time.monotonic()-began<52,'Declared landmark stage budget nearly exhausted'
        target=meshes[pid];tree=BVH(target);rows=[tree.closest_point(p)for p in centers]
        xyz=np.array([r['point']for r in rows]);ids=np.array([r['face_id']for r in rows]);bary=np.array([r['barycentric']for r in rows]);distance=np.array([r['distance_m']for r in rows]);tri=target.triangles[ids];norm=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm/=np.maximum(np.linalg.norm(norm,axis=1)[:,None],1e-300)
        witnesses[f'target_{j}_source_face_ids']=ids;witnesses[f'target_{j}_barycentric']=bary;witnesses[f'target_{j}_surface_points']=xyz;witnesses[f'target_{j}_distance_m']=distance;witnesses[f'target_{j}_source_face_normals']=norm;witnesses[f'target_{j}_point_minus_surface']=centers-xyz
        summaries.append(dict(key=j,part_id=pid,minimum_center_to_surface_m=float(distance.min()),minimum_at_exact_source_station=int(distance.argmin()),range_m=[float(distance.min()),float(distance.max())],qualifier='Actual nearest source-triangle witness only; not attachment, patent interface, tissue-volume proof or branch identity.'))
    atomic_npz(DEST/'left-parent-landmark-witnesses.npz',**witnesses)
    atomic_json(DEST/'source-preparation.json',dict(status='completed-read-only',catalogue=str(CAT),catalogue_sha256=sha,ordinary_catalogue=str(ORD),ordinary_catalogue_sha256=digest_file(ORD),frame=frame,source_bindings=bindings,source_banks=banks,parent_course=dict(part_id=parent['id'],components=len(courses),stations=len(rings),controls_per_ring=rings.shape[1],source_faces=len(f),mouth_station_ids=[0,len(rings)-1],correspondence_proof=proof,stations=str(DEST/'left-parent-stations.npz'),stations_sha256=digest_file(DEST/'left-parent-stations.npz')),landmark_witnesses=summaries,landmark_archive=str(DEST/'left-parent-landmark-witnesses.npz'),landmark_archive_sha256=digest_file(DEST/'left-parent-landmark-witnesses.npz'),anatomical_boundary=dict(status='not-yet-established',reason='Native parent is a single serial course without a source branch node. Actual obturator-externus source combines muscle/tendon identity and does not itself define the named tendon border. Nearest surfaces cannot promote the distal course to the deep branch.',required_semantics='Deep course posterior to obturator externus tendon, anterior to quadratus femoris, near proximal quadratus border toward posterior hip capsule; root must inspect actual landmarks before selecting an exact shared source station.',no_station_fraction_or_generic_y_bin=True),primary_references=[dict(url='https://pubmed.ncbi.nlm.nih.gov/10963165/',scope='24-cadaver study; qualitative deep-branch relationships, no patient template copied'),dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11195040/',scope='Primary periacetabular CTA; qualitative deep-branch capsule/acetabulum recipient relationships, not numeric fitting')],unchanged_master_verified=digest_file(CAT)==sha,skeleton_moved=False,new_geometry=False,installed=False,elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(output=str(DEST),sources=len(bindings),ordinary_exact=sum(bool(x['ordinary_body_match'])and x['ordinary_body_match']['whole_geometry_exact']and x['ordinary_body_match']['all_attributes_exact']for x in bindings),stations=len(rings),missing_terms=5,existing_named_missing_branch_geometry=0,elapsed_seconds=time.monotonic()-began)))

if __name__=='__main__':main()
