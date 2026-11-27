"""Bounded current-source renal blood/interface inventory; constructs no geometry."""
import json, sqlite3, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.audit import topology

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
BODY=ROOT/'.local-models/organ-models/spinal-cord/male-body.json'
AUTH=ROOT/'.local-models/ta98-authority.sqlite'
CENSUS=ROOT/'reports/ta98-blood-circulation-inventory-installed-v12.json'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
FORK=['ZA:Left renal artery','ZA:Anterior branch of renal artery.l','ZA:Posterior branch of renal artery.l']
NETWORK=['ZA:Intrarenal arteries of left kidney','ZA:Intrarenal veins of left kidney']
VESSELS=FORK+['ZA:Left renal vein']+NETWORK
TISSUE=['RENALDETAIL:male:left:'+n for n in ['capsule','cortex']+['pyramid:'+c for c in 'abcdefghijk']+['papilla:'+c for c in 'abcdefghijk']+['collecting-wall:'+c for c in ['adventitia','muscularis','mucosa']]]

def area(p):
    q=p-p.mean(0);return float(np.linalg.norm(np.cross(q,np.roll(q,-1,axis=0)).sum(0))/2)

def physical(m):
    v,inv=np.unique(m.vertices,axis=0,return_inverse=True);return v,inv[m.faces]

def mouth(p):
    return dict(center_m=p.mean(0).tolist(),controls=len(p),area_m2=area(p),
                equivalent_diameter_m=float(np.sqrt(4*area(p)/np.pi)))

def match(a,b):
    if len(a)!=len(b):return None
    best=None
    for rev in [False,True]:
        for shift in range(len(b)):
            order=np.roll(np.arange(len(b))[::-1]if rev else np.arange(len(b)),shift)
            d=np.linalg.norm(a-b[order],axis=1);r=dict(max_gap_m=float(d.max()),rms_gap_m=float(np.sqrt(np.mean(d*d))),literal=bool(np.array_equal(a,b[order])),reverse=rev,phase=shift)
            if best is None or r['max_gap_m']<best['max_gap_m']:best=r
    return best

def main():
    began=time.monotonic();assert not OUT.exists();assert digest_file(CAT)==EXPECTED
    OUT.mkdir(parents=True);a=json.loads(CAT.read_text());body=json.loads(BODY.read_text())
    parts={p['id']:p for p in a['parts']};bodyparts={p['id']:p for p in body['parts']};concepts={c['id']:c for c in a['concepts']}
    census=json.loads(CENSUS.read_text());ct={r['id']:r for r in census['terms']}
    c=sqlite3.connect(AUTH);c.row_factory=sqlite3.Row
    rows=list(c.execute("select * from ta98 where id like 'A08.1.03%' or id like 'A08.1.04%' or id between 'A12.2.12.075' and 'A12.2.12.085' or id between 'A12.3.09.009' and 'A12.3.09.012M' order by id"))
    inventory=[]
    micro={'A08.1.03.005','A08.1.03.006','A08.1.03.008','A08.1.04.005'}
    for r in rows:
        r=dict(r)
        if r['female_gender'] and not r['male_gender']:continue
        pid=[p['id']for p in parts.values()if not p.get('suppressed')and(p.get('ta98Term')or concepts.get(p.get('conceptId'),{}).get('ta98Term'))==r['id']]
        cr=ct.get(r['id']);available=cr.get('active_exact_source_part_ids',[])if cr else pid
        inventory.append(dict(id=r['id'],name=r['name_en'],parent=r['parent_id'],classification='microscopic-vessel'if r['id']in micro else 'macroscopic-vessel-network'if r['id'].startswith('A08')else 'named-arterial-or-venous-branch',
                              census_included=bool(cr),census_exact_source_ids=available,raw_current_exact_source_ids=pid,
                              missing_distinct_source=not bool(available),raw_name_alone_not_identity=True))
    atomic_json(OUT/'inventory.json',dict(authority=str(AUTH),authority_sha256=digest_file(AUTH),census=str(CENSUS),census_sha256=digest_file(CENSUS),catalogue_sha256=EXPECTED,terms=inventory,
                                         a08_omission_from_current_census=not any(x['census_included']for x in inventory if x['id'].startswith('A08')),
                                         denominator_modified=False))
    reader=AtlasReader(ROOT);meshes={};bindings=[];mouths={};physical_data={}
    for i,pid in enumerate(VESSELS+TISSUE):
        assert time.monotonic()-began<48,'Source decode nearing declared deadline'
        p=parts[pid];assert not p.get('suppressed');assert pid in bodyparts and not bodyparts[pid].get('suppressed')
        m=reader.read_part(a,p,frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED);meshes[pid]=m
        bm=reader.read_part(body,bodyparts[pid],frame=m.frame,catalogue_path=BODY,catalogue_sha=digest_file(BODY));assert bm.geometry_sha==m.geometry_sha and bm.attributes_sha==m.attributes_sha
        v,f=physical(m);physical_data[pid]=(v,f)
        source=OUT/('source-'+str(i).zfill(2)+'.npz');m.save(source)
        bindings.append(dict(part_id=pid,mesh=str(source),mesh_sha256=digest_file(source),geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,faces=len(f),bank=m.source['bank_path'],bank_sha256=m.source['bank_sha256'],normal_body_geometry_and_attributes_literal=True))
        if pid in VESSELS:
            topo=topology(m);loops=boundary_loops(f);mouths[pid]=[v[l]for l in loops]
            atomic_json(OUT/('vessel-'+str(i)+'.json'),dict(part_id=pid,topology=topo,mouths=[mouth(v[l])for l in loops],finite_wall_established=pid in FORK))
    # External fork annuli are literal authored tissue draw groups, not nearest surfaces.
    ports=[]
    for i,pid in enumerate(FORK):
        p=parts[pid];m=meshes[pid];groups=p['tissueDrawGroups'];g=next(g for g in groups if not g.get('neighbor')and g['count']==72)
        fs=m.faces[g['start']//3:(g['start']+g['count'])//3];v,inv=np.unique(m.vertices,axis=0,return_inverse=True);af=inv[fs]
        rings=[v[l]for l in boundary_loops(af)];assert len(rings)==2 and all(len(r)==12 for r in rings)
        rings.sort(key=area,reverse=True);path=OUT/('fork-port-'+str(i)+'.npz');atomic_npz(path,outer=rings[0],inner=rings[1],annulus_vertices=v,annulus_faces=af)
        comparisons=[]
        for j,r in enumerate(mouths[NETWORK[0]]):
            z=match(rings[0],r)
            if z:comparisons.append(dict(network_mouth=j,**z))
        ports.append(dict(owner=pid,port_mesh=str(path),port_sha256=digest_file(path),outer=mouth(rings[0]),inner=mouth(rings[1]),
                          anatomical_role='aortic side-entry (unrepaired)'if i==0 else 'distal '+('anterior'if i==1 else'posterior')+' outlet; hilar passage unverified',
                          native_intrarenal_outer_mouth_best_match=min(comparisons,key=lambda z:z['max_gap_m']),join_inferred=False))
    # Test literal common facets of current finite cortex and pyramids.
    cp=meshes['RENALDETAIL:male:left:cortex'];cv,cf=physical_data['RENALDETAIL:male:left:cortex']
    def face_keys(v,f):return {b''.join(sorted(v[t].astype('<f8').tobytes()[k:k+24]for k in [0,24,48]))for t in f}
    ck=face_keys(cv,cf);pyramids=[]
    for letter in 'abcdefghijk':
        pid='RENALDETAIL:male:left:pyramid:'+letter;v,f=physical_data[pid];pap=meshes['RENALDETAIL:male:left:papilla:'+letter]
        shared=len(ck&face_keys(v,f));pyramids.append(dict(part_id=pid,faces=len(f),closed_oriented=topology(meshes[pid])['closed_oriented_manifold'],
                        bounds_m=[v.min(0).tolist(),v.max(0).tolist()],papilla_center_m=pap.vertices.mean(0).tolist(),
                        exact_common_cortex_triangles=shared,anatomically_defined_basal_patch=False,
                        basal_constraint='Actual whole pyramid surface and papilla are available; no authored cortex-medulla basal patch selector is present. Do not substitute an arbitrary coordinate plane.'))
    pilot=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-renal-source-pilot-v1'
    source_replay=[]
    for idx,pid in [(4,NETWORK[0]),(5,NETWORK[1])]:
        old=json.loads((pilot/('part-'+str(idx).zfill(2)+'.json')).read_text())
        source_replay.append(dict(part_id=pid,prior_report=str(pilot/('part-'+str(idx).zfill(2)+'.json')),literal_geometry_replayed=old['geometry_sha256']==meshes[pid].geometry_sha,
                                 literal_attributes_replayed=old['attributes_sha256']==meshes[pid].attributes_sha,components=len(old['topology']['components']),outer_mouths=len(old['mouths']),finite_wall=False,
                                 old_anatomical_or_packing_pass_reused=False))
    result=dict(schema='current-renal-blood-frontier/v1',catalogue=str(CAT),catalogue_sha256=EXPECTED,body_catalogue=str(BODY),body_catalogue_sha256=digest_file(BODY),
                source_bindings=bindings,external_fork_ports=ports,native_network_replay=source_replay,pyramids=pyramids,
                cortex=dict(part_id='RENALDETAIL:male:left:cortex',closed_oriented=topology(cp)['closed_oriented_manifold'],peripheral_only=True,inter_pyramid_columns_complete=False),
                collecting_protection=dict(parts=[p for p in TISSUE if ':collecting-wall:'in p],mucosal_cell_is_tissue_not_urine_volume=True,
                    lumen_boundary_selector_present=False,need='Recover actual innermost mucosal free surface and open calyceal mouths; no vessel may cross that lining or occupy collecting lumen. Outer finite tissue contact alone does not establish this.'),
                new_geometry=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',result)
    assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(output=str(OUT/'evidence.json'),sources=len(bindings),terms=len(inventory),ports=ports,native_replay=source_replay,cortex_pyramid_exact_shared_faces=sum(x['exact_common_cortex_triangles']for x in pyramids),elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
