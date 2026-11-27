"""Bounded current-source export and anatomical FOM construction contract.

Exports only this project's current atlas guides. This is not a tissue
generator, never moves source anatomy and does not infer a frenulum cord.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sqlite3
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.spatial import BVH
from geometry_tools.ta98_authority import corrected_terms
from tongue_whole_input import WORK,CATALOGUE

FRAME='model:local-ta98:assembled-metres'
FIXED=['ZA:Mandible:Bone','ZA:Mylohyoid muscle.l','ZA:Mylohyoid muscle.r',
       'ZA:Geniohyoid muscle.l','ZA:Geniohyoid muscle.r',
       'ZA:Sublingual gland.l','ZA:Sublingual gland.r',
       'ZA:Submandibular duct.l','ZA:Submandibular duct.r',
       'ZA:Lingual nerve.l','ZA:Lingual nerve.r',
       'ZA:Mylohyoid muscle.ol','ZA:Mylohyoid muscle.or',
       'ZA:Geniohyoid muscle.ol','ZA:Geniohyoid muscle.or']
ROOT_TERMS=['A05.1.01.002','A05.1.01.101','A05.1.04.012','A05.1.04.008',
            'A05.1.04.007','A05.1.04.011','A05.1.02.008','A05.1.02.011',
            'A04.2.03.006','A04.2.03.007','A02.1.15.002','A14.2.01.081']


def bound(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=digest_file(p))


def classify(code):
    if code in ['A05.1.01.002','A05.1.04.011']:return'wall-layer'
    if code in ['A05.1.04.007']:return'surface-region'
    if code in ['A05.1.01.101']:return'whole-organ'
    return'tissue-region'


def main(output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local stage')
    atlas=json.loads(CATALOGUE.read_text());sha=digest_file(CATALOGUE);parts={p['id']:p for p in atlas['parts']}
    if any(pid not in parts for pid in FIXED):raise ValueError('A declared fixed guide is missing')
    authority_path=ROOT/'.local-models/ta98-authority.sqlite'
    with sqlite3.connect('file:'+str(authority_path)+'?mode=ro',uri=True)as conn:
        conn.row_factory=sqlite3.Row;authority=corrected_terms({r['id']:dict(r)for r in conn.execute('SELECT * FROM ta98')})
    wanted=set(ROOT_TERMS)
    while True:
        expanded=wanted|{r['id']for r in authority.values()if r['parent_id']in wanted}
        if expanded==wanted:break
        wanted=expanded
    terms=[]
    for code in sorted(wanted):
        mapped=[]
        for c in atlas['concepts']:
            if c.get('ta98Term')==code:
                mapped.extend(c['elements'])
        terms.append(dict(authority=authority[code],representation_class=classify(code),
            source_mesh_ids=sorted(set(pid for pid in mapped if pid in parts)),
            source_identity_not_anatomical_acceptance=True))
    output.mkdir(parents=True);(output/'meshes').mkdir()
    reader=AtlasReader();meshes={};exports=[]
    for index,pid in enumerate(FIXED+['ZA:Tongue']):
        part=parts[pid];mesh=reader.read_part(atlas,part,frame=FRAME,catalogue_path=CATALOGUE,catalogue_sha=sha)
        path=output/'meshes'/f'current-{index:02d}.npz';mesh.save(path);meshes[pid]=mesh
        exports.append(dict(part_id=pid,name=part['name'],mesh=bound(path),geometry_sha256=mesh.geometry_sha,
            frame=mesh.frame,source_bank=bound(mesh.source['bank_path']),source_record=part,
            fixed=True,topology=topology(mesh),source_catalogue_sha256=sha,
            visibility_hidden=pid in set(atlas.get('defaultHidden',[])),suppressed=bool(part.get('suppressed'))))
    # Explicit selected construction guides supplement the current baseline;
    # their authority is not confused with the not-yet-installed catalogue.
    selected_path=WORK/'selected-coupling-audit-v2/coupled-inputs.selected.json'
    selected=json.loads(selected_path.read_text());mucosa_path=Path(selected['covering_meshes']['mucosa']['path'])
    if digest_file(mucosa_path)!=selected['covering_meshes']['mucosa']['sha256']:raise ValueError('Selected mucosa changed')
    mucosa=Mesh.load(mucosa_path)
    tissue_path=mucosa_path.parent/'tissue.json';tissue=json.loads(tissue_path.read_text())
    with np.load(tissue['arrangement'],allow_pickle=False)as raw:
        positions=raw['positions'].copy();outer=raw[tissue['surfaces']['mucosa:outer']].copy()
    edges=np.r_[outer[:,[0,1]],outer[:,[1,2]],outer[:,[2,0]]]
    sorted_edges=np.sort(edges,axis=1);_,inv,count=np.unique(sorted_edges,axis=0,return_inverse=True,return_counts=True)
    rim_edges=edges[count[inv]==1]
    adjacency={}
    for a,b in rim_edges:adjacency.setdefault(int(a),[]).append(int(b));adjacency.setdefault(int(b),[]).append(int(a))
    if any(len(n)!=2 for n in adjacency.values()):raise ValueError('Exposed mucosal carrier lacks one simple reflection loop')
    loop=[min(adjacency)];last=None
    while True:
        now=loop[-1];choices=[q for q in adjacency[now]if q!=last]
        nxt=choices[0]
        if nxt==loop[0]:break
        loop.append(nxt);last=now
        if len(loop)>len(adjacency):raise ValueError('Reflection loop invalid')
    if len(loop)!=len(adjacency):raise ValueError('Multiple reflection loops need explicit separate anatomy')
    reflection=positions[loop];atomic_npz(output/'actual-tongue-reflection.npz',positions=reflection,
        source_position_indices=np.array(loop),source_outer_faces=outer,frame=np.array(FRAME))
    # Current source topology can justify exposed duct-rim guides, not a blind
    # tip, lumen or oral opening that has never been constructed.
    duct_boundaries=[]
    for side in ['l','r']:
        pid='ZA:Submandibular duct.'+side;m=meshes[pid];tri=m.vertices[m.faces]
        ee=np.r_[m.faces[:,[0,1]],m.faces[:,[1,2]],m.faces[:,[2,0]]]
        se=np.sort(ee,axis=1);uniq,counts=np.unique(se,axis=0,return_counts=True)
        boundary=uniq[counts==1];points=m.vertices[np.unique(boundary)]if len(boundary)else np.empty((0,3))
        path=output/f'actual-submandibular-{side}-open-edges.npz'
        atomic_npz(path,positions=m.vertices,edges=boundary,frame=np.array(FRAME))
        duct_boundaries.append(dict(part_id=pid,boundary_edges=len(boundary),boundary_vertices=len(points),
            witness=bound(path),actual_orifice_identified=False,
            finding='Open edge topology is a source rim candidate only; closed tip cannot establish a salivary ostium',
            bounds_m=None if len(points)==0 else np.array([points.min(0),points.max(0)]).tolist()))
    # Actual nearest cortical witnesses to the superior gland roof identify
    # lingual mandibular vicinity. They do not invent a fascia insertion line.
    mandible=meshes['ZA:Mandible:Bone'];index=BVH(mandible);witnesses=[]
    for side in ['l','r']:
        gland=meshes['ZA:Sublingual gland.'+side];tri=gland.triangles;centers=tri.mean(1)
        normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(normal,axis=1)
        normal/=np.maximum(length,1e-30)[:,None]
        ids=np.flatnonzero((normal[:,1]>.3)&(centers[:,1]>np.quantile(centers[:,1],.70)))
        ids=ids[np.argsort(centers[ids,2])]
        choose=np.unique(np.rint(np.linspace(0,len(ids)-1,5)).astype(int))
        for i in ids[choose]:
            hit=index.closest_point(centers[i])
            witnesses.append(dict(side=side,gland_part_id='ZA:Sublingual gland.'+side,
                gland_face=int(i),gland_barycentric=[1/3]*3,gland_point_m=centers[i].tolist(),
                mandibular_face=int(hit['face_id']),mandibular_barycentric=hit['barycentric'].tolist(),
                mandibular_point_m=hit['point'].tolist(),gap_m=float(hit['distance_m']),
                attachment_interpretation='Candidate lingual mandibular neighborhood; superior fascial/gingival attachment height not established'))
    rows=[]
    for pid in ['ZA:Mylohyoid muscle.l','ZA:Mylohyoid muscle.r','ZA:Geniohyoid muscle.l',
                'ZA:Geniohyoid muscle.r','ZA:Sublingual gland.l','ZA:Sublingual gland.r']:
        m=meshes[pid];rows.append(dict(part_id=pid,bounds_m=m.bounds.tolist(),
            superior_source_point_m=m.vertices[int(np.argmax(m.vertices[:,1]))].tolist(),
            interpretation='Bounds and superior source control only; no constructed fascial contact inferred'))
    source_refs=[dict(url='https://onlinelibrary.wiley.com/doi/10.1002/ca.23343',
        evidence='Nine fresh adult dissections: FOM fascia/mucosa reflect from broad inner mandibular arc to anterior/lateral ventral tongue; GG, sublingual glands/ducts are suspended deep to fascia.'),
        dict(url='https://onlinelibrary.wiley.com/doi/full/10.1155/2020/1820978',
        evidence='Adult histology: continuous fascia with variable thickness and mucosal incorporation; no discrete AP-oriented midline cord.')]
    contract=dict(schema='tongue-floor-of-mouth-construction-contract/v1',frame=FRAME,
        source_catalogue=bound(CATALOGUE),authority=bound(authority_path),
        corrected_authority_code=bound(ROOT/'geometry_tools/ta98_authority.py'),
        inventory=terms,inventory_roots=ROOT_TERMS,source_exports=exports,
        selected_tongue_assembly=dict(selection=bound(selected_path),mucosa=bound(mucosa_path),
            geometry_sha256=mucosa.geometry_sha,current_catalogue_installation_not_inferred=True),
        observed_guides=rows,mandibular_vicinity_witnesses=witnesses,
        tongue_reflection=dict(source=bound(output/'actual-tongue-reflection.npz'),vertices=len(loop),
            existing_outer_coordinates_literal=True,existing_anatomical_boundary_provisional=True,
            formula_reference=bound(WORK/'exposed-mucosa-v1/reflection-boundary.json'),
            warning='The exact existing computational seam is usable for coupled continuity, not a donor-confirmed frenulum attachment.'),
        duct_rim_candidates=duct_boundaries,primary_references=source_refs,
        assembly=dict(whole='Anterior and bilateral lateral floor-of-mouth mucosa/fascial roof continuous with ventral tongue',
            planned_parts=[dict(name='FOM mucosal field',term='A05.1.01.002',representation='finite-wall-layer'),
                dict(name='FOM fascia',term=None,representation='additional finite connective tissue field; no invented TA98 code'),
                dict(name='Lingual frenulum',term='A05.1.04.012',representation='one continuous selected fold of those fields, not an additional cord'),
                dict(name='Fimbriated folds',term='A05.1.04.008',representation='paired lobulated ventral mucosal relief after actual frenulum extremity is established'),
                dict(name='Sublingual caruncles',term='A05.1.01.112',representation='finite mucosal elevations around actual duct mouths'),
                dict(name='Sublingual folds',term='A05.1.01.113',representation='mucosal ridge over fixed gland/duct fields')],
            fixed=FIXED,
            boundaries=[dict(name='anterior/lateral mandibular mucogingival attachment',
                source='Actual inner mandibular cortex + current gingiva/attachment guides',
                state='Cortical vicinity witnessed; continuous mucogingival line and insertion height still require identification'),
                dict(name='tongue reflection',source='Literal current exposed mucosal loop',
                     state='Coordinates available; replace its closure rim with reciprocal FOM mucosa/fascia interfaces'),
                dict(name='posterolateral termination',source='Lateral ventral tongue near HG insertion and oral mucosal reflection',
                     state='Need paired continued side margins, not an invented complete posterior floor disk')],
            openings=[dict(name='paired submandibular duct ostia',
                rule='Actual luminal mouths must pass through finite fascia/mucosa at caruncles; annular rim only, no solid plugs',
                state='Source duct rim topology evaluated; oral mouth identity and patent finite lumen remain unresolved'),
                dict(name='minor sublingual duct ostia',
                rule='Multiple genuine openings along sublingual folds; retain as missing until ducts/orifices are constructed',state='missing')],
            neighbor_interfaces=[
                dict(parts=['mandible','FOM fascia'],relation='adjacent-unjoined until witnessed broad periosteal contact'),
                dict(parts=['FOM mucosa','FOM fascia'],relation='literal reciprocal finite interface'),
                dict(parts=['ventral tongue mucosa','FOM mucosa'],relation='continuous reciprocal joined seam; no overlapping closure cap'),
                dict(parts=['sublingual glands','FOM fascia'],relation='suspension/covering; surface contact and finite volume exclusion must be measured'),
                dict(parts=['genioglossus','FOM fascia'],relation='deep suspension; do not bisect muscle with a midline fascia wall'),
                dict(parts=['lingual nerves','FOM fascia'],relation='immediately deep; preserve source geometry and test mucosal exclusion'),
                dict(parts=['mylohyoid/geniohyoid','FOM roof'],relation='deep inferior support guides; not interchangeable with mucosa')]),
        feasibility=dict(continuous_finite_carrier_anatomically_supported=True,
            source_specific_ready_to_construct=False,
            available_geometry='Native skeletal/muscle/gland/duct/nerve guides and literal provisional tongue reflection loop',
            missing_geometry=['Actual inner mandibular mucogingival attachment curve at superior fascia level',
                'Actual paired submandibular oral mouth identity/finite lumen',
                'Donor-consistent fascia/mucosal thickness field and posterolateral fold termination',
                'Explicit reciprocal tongue rim opening and attachment normal charts'],
            next_construction='First author a curved source-guided annular/horse-shoe FOM roof from a witnessed inner-jaw curve to literal ventral-tongue seam, with duct channels and finite shared fascia/mucosa. Elevate the central existing field to form frenulum only after these boundaries are fixed.'),
        no_new_tissue_constructed=True,no_existing_mesh_modified=True,fixed_skeleton=True,
        vendor_mesh_copy=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'contract.json',contract)
    colors={'skeletal':[.72,.70,.60,.18],'muscular':[.65,.31,.33,.45],
            'digestive':[.83,.63,.34,.75],'nervous':[.78,.73,.36,.8],'attachments':[.58,.67,.76,.65]}
    render=[]
    for row in exports:
        if row['part_id']=='ZA:Tongue':continue
        render.append(dict(mesh=row['mesh']['path'],name=row['name'],color=colors.get(row['source_record']['system'],[.6,.6,.6,.5])))
    render.append(dict(mesh=str(mucosa_path),name='Selected actual ventral tongue mucosa (provisional reflection)',color=[.75,.38,.42,.4]))
    atomic_json(output/'review-spec.json',dict(frame=FRAME,resolution=640,wall_seconds=50,
        meshes=render,views=[dict(name='superior-context',direction=[0,1,0]),
            dict(name='right-context',direction=[1,0,0]),
            dict(name='sagittal-context',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.622,.060],clear_outer=True))],
        qualification='Source context only. No FOM roof, frenulum, duct ostia or new tissue constructed.'))
    print(json.dumps(dict(contract=str(output/'contract.json'),authority_terms=len(terms),current_source_exports=len(exports),
        actual_reflection_controls=len(loop),duct_open_edges=[r['boundary_edges']for r in duct_boundaries],
        source_specific_ready_to_construct=False,elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();main(a.output.resolve())
