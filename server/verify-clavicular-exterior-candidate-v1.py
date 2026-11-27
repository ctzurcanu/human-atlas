"""Independent replay of the emitted clavicular course and source ownership.

Run only on a completed reviewed course/CSG package; this never installs it.
"""
import argparse, hashlib, json, sys, time
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops

PARENT='ZA:Thoraco-acromial artery.l'
CHILD='TA98PART:ZA:Thoraco-acromial artery.l:A12.2.09.008'
FRAME='model:local-ta98:assembled-metres'
H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def canonical_rows(p,n,u):
    return min(b''.join((np.roll(p,k,axis=0).astype('<f8').tobytes(),np.roll(n,k,axis=0).astype('<i2').tobytes(),np.roll(u,k,axis=0).astype('<f4').tobytes())) for k in range(3))
def corners(mesh):
    return Counter(canonical_rows(mesh.vertices[f],mesh.attributes['normals'][f],mesh.attributes['uvs'][f]) for f in mesh.faces)
def edges(f):
    e=np.vstack((f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]));a,c=np.unique(np.sort(e,axis=1),axis=0,return_counts=True)
    return {tuple(x) for x in a[c==1]}
def run(out,side='left'):
    global PARENT,CHILD
    if side=='right':PARENT='ZA:Thoraco-acromial artery.r';CHILD='TA98PART:ZA:Thoraco-acromial artery.r:A12.2.09.008'
    started=time.monotonic();proof=json.loads((out/'packing-proof.json').read_text());assert proof['passed'] and proof['complete']
    proposed=out/'proposed-atlas.json';assert H(proposed)==proof['proposed_catalogue_sha256']
    master=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert H(master)==proof['source_catalogue_sha256']
    before=json.loads(master.read_text());after=json.loads(proposed.read_text());old={p['id']:p for p in before['parts']};new={p['id']:p for p in after['parts']}
    assert new.keys()-old.keys()=={CHILD} and not old.keys()-new.keys()
    assert all(new[k]==v for k,v in old.items() if k!=PARENT)
    assert after['materials']==before['materials'] and after['concepts'][:-1]==before['concepts'] and after['ta98TermBindings'][:-1]==before['ta98TermBindings']
    assert after['concepts'][-1]['ta98Term']=='A12.2.09.008' and after['concepts'][-1]['ta98Kind']=='exact'
    assert after['concepts'][-1]['parentId']==PARENT
    reader=AtlasReader();source=reader.read_part(before,old[PARENT],frame=FRAME);parent=reader.read_part(after,new[PARENT],frame=FRAME);child=reader.read_part(after,new[CHILD],frame=FRAME)
    z=np.load(out/'compiled-shared-geometry.npz');v,f,own=z['vertices'],z['faces'],z['face_owner'];nn,uv=z['corner_normal_words'],z['corner_uvs']
    for role,mesh in [(0,parent),(1,child)]:
        expected=Counter(canonical_rows(v[face],n,u) for face,n,u in zip(f[own==role],nn[own==role],uv[own==role]))
        assert corners(mesh)==expected,'Emitted decoded corners differ from reviewed shared geometry'
    plan=json.loads((out/'assembly-plan-before-union.json').read_text());allowed=set(plan['course_spec']['allowed_parent_face_ids']);retained=corners(parent)
    outside=[i for i in range(len(source.faces)) if i not in allowed]
    for i in outside:
        tri=source.faces[i];key=canonical_rows(source.vertices[tri],source.attributes['normals'][tri],source.attributes['uvs'][tri]);assert retained[key]==1,'Source parent chart outside ostium changed'
    ep,ec=edges(f[own==0]),edges(f[own==1]);shared=ep&ec;assert shared
    degree=Counter(x for e in shared for x in e);assert set(degree.values())=={2},'Shared ostium must be a closed curve'
    unseen=set(degree);stack=[unseen.pop()]
    while stack:
        a=stack.pop()
        for edge in shared:
            if a not in edge:continue
            b=edge[0] if edge[1]==a else edge[1]
            if b in unseen:unseen.remove(b);stack.append(b)
    assert not unseen,'Parent-child seam has multiple disjoint curves'
    for vi in degree:
        words=[]
        for role in [0,1]:
            fi,ci=np.where(f==vi);ids=own[fi]==role;words.append({tuple(x) for x in nn[fi[ids],ci[ids]]})
        assert len(words[0])==len(words[1])==1 and words[0]==words[1]
    merged=Mesh(v,f,FRAME);check=audit(merged,require_closed=False,expected_components=1,max_pairs=600000,max_tests=2000000);assert check['passed'];loops=boundary_loops(f);assert len(loops)==4
    original={frozenset(map(tuple,source.vertices[x])) for x in boundary_loops(source.faces)};actual={frozenset(map(tuple,v[x])) for x in loops};assert original<=actual,'Native end mouths changed or capped'
    bound=any(x.get('term')=='A12.2.09.008'for x in before['ta98TermBindings'])or any(x.get('ta98Term')=='A12.2.09.008'for x in before['concepts'])
    report=dict(schema='independent-clavicular-compiled-replay/v1',passed=True,complete=True,source_catalogue_sha256=H(master),proposed_catalogue_sha256=H(proposed),unchanged_existing_parts=len(old)-1,new_parts=[CHILD],new_TA98_terms=[]if bound else['A12.2.09.008'],construction_side=side,all_original_concepts_and_materials_literal=True,parent_faces_outside_ostium_literal=len(outside),single_shared_open_ostium_edges=len(shared),shared_packed_normal_words_literal=True,all_compiled_decoded_faces_owned_once=True,native_parent_mouths_literal=True,outer_mouths=4,merged_audit=check,source_outer_surface_only=True,finite_wall_or_complete_perfusion=False,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started)
    dump(out/'independent-compiled-replay.json',report);print(json.dumps({k:v for k,v in report.items() if k!='merged_audit'}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);p.add_argument('--side',choices=['left','right'],default='left');a=p.parse_args();run(a.directory.resolve(),a.side)
