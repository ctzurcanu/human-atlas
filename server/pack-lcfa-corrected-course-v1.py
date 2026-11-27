"""Pack saved LCFA course corrections while preserving source span ownership."""
from pathlib import Path
import copy, hashlib, json, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_bytes,atomic_json,digest_file
from geometry_tools.model import AtlasReader

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1'
OUT=BASE/'course-proposal-v1'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='b26860231a10b4b4382802751d6a8c15a7491e8baa83d629d73b93f5abb61707'

def main():
    began=time.monotonic();assert digest_file(SOURCE)==EXPECTED
    atlas=json.loads(SOURCE.read_bytes());proposed=copy.deepcopy(atlas);parts={p['id']:p for p in proposed['parts']};reader=AtlasReader()
    proof=dict(schema='saved-lcfa-course-literal-span-pack/v1',source_catalogue=str(SOURCE),source_catalogue_sha256=EXPECTED,passed=False,installed=False,anatomical_acceptance=False,sides=[])
    blob=bytearray();changed=[]
    for side in ['l','r']:
        owner=parts['ZA:Lateral circumflex femoral artery.'+side];child=parts['TA98NATIVE:male:LCFA:ascending:'+side]
        for field in ['positions','normals','uvs','vertexCount','chunk']:assert owner[field]==child[field]
        assert child['indices']==owner['indices']+owner['indexCount']*4
        assert not any(key in owner for key in ['geometryScale','geometryOrigin','geometryInstance','sourceOffset','normalSigns'])
        bankpath=reader.asset(atlas['chunks'][owner['chunk']]['url']);bank=reader._bank(bankpath);n=owner['vertexCount']
        oldv=reader.view(bank,owner['positions'],n,'<f4',3).copy();oldn=reader.view(bank,owner['normals'],n,'<i2',3).copy();uv=reader.view(bank,owner['uvs'],n,'<f4',2).copy()
        oldf=reader.view(bank,owner['indices'],owner['indexCount'],'<u4').reshape(-1,3).copy();ascendingf=reader.view(bank,child['indices'],child['indexCount'],'<u4').reshape(-1,3).copy()
        saved=np.load(BASE/f'bilateral-course-proof-v2/{side}-course.npz');before=np.load(BASE/'course-correction-v3/candidate.npz')['source_vertices'] if side=='l' else np.load(BASE/'source-gate-v1/r-component-0.npz')['vertices']
        lookup={tuple(x):i for i,x in enumerate(before)};assert len(lookup)==n
        mapping=np.array([lookup[tuple(x)]for x in oldv.astype(float)])
        assert np.array_equal(before[mapping][oldf],np.load(BASE/f'source-gate-v1/{side}-component-0.npz')['vertices'][np.load(BASE/f'source-gate-v1/{side}-component-0.npz')['faces']])
        newv=saved['vertices'][mapping].astype('<f4');newn=saved['normal_words'][mapping].astype('<i2');assert np.array_equal(uv,saved['uvs'][mapping])
        assert np.array_equal(newv[np.unique(ascendingf)],oldv[np.unique(ascendingf)])
        assert np.array_equal(newn[np.unique(ascendingf)],oldn[np.unique(ascendingf)])
        assert np.array_equal(newv[oldf].astype(float),saved['vertices'][saved['faces']])
        assert np.array_equal(newn[oldf],saved['normal_words'][saved['faces']])
        offsets={}
        for name,array in [('positions',newv),('normals',newn),('uvs',uv),('indices',np.concatenate([oldf,ascendingf]).astype('<u4'))]:
            while len(blob)%4:blob.extend(b'\0')
            offsets[name]=len(blob);blob.extend(np.ascontiguousarray(array).tobytes())
        for p,f in [(owner,oldf),(child,ascendingf)]:
            p.update(chunk=len(atlas['chunks']),positions=offsets['positions'],normals=offsets['normals'],uvs=offsets['uvs'],indices=offsets['indices']+(owner['indexCount']*4 if p is child else 0))
            points=newv[np.unique(f)].astype(float);p['bounds']=[points.min(0).tolist(),points.max(0).tolist()];changed.append(p['id'])
        owner['sourceOwnershipQualification']='Coarse root/transverse course corrected around the fixed femur; zero detected femur crossings and at least 0.25 mm sampled control clearance. Original parent/fork contacts retained. No transverse-branch identity claim; finite wall and communicating lumen unverified.'
        proof['sides'].append(dict(side=side,source_bank=str(bankpath),source_bank_sha256=digest_file(bankpath),source_vertex_order_retained=True,all_original_index_spans_retained=True,ascending_positions_normals_uvs_literal=True,course_corners_equal_saved_arrays=True,vertex_count=n,triangles=len(oldf)+len(ascendingf),changed_controls=int(np.count_nonzero(np.any(newv!=oldv,axis=1))),maximum_displacement_m=float(np.linalg.norm(newv.astype(float)-oldv,axis=1).max())))
    sha=hashlib.sha256(blob).hexdigest();bankpath=ROOT/'.local-models/ta98-runtime'/('lcfa-course-'+sha[:16])/'geometry.bin';atomic_bytes(bankpath,bytes(blob));proposed['chunks'].append(dict(bytes=len(blob),url='/local-models/'+str(bankpath.relative_to(ROOT/'.local-models'))))
    data=(json.dumps(proposed,separators=(',',':'))+'\n').encode();atomic_bytes(OUT/'source-atlas.json',SOURCE.read_bytes());atomic_bytes(OUT/'proposed-atlas.json',data)
    loaded=json.loads(data);emitted={p['id']:p for p in loaded['parts']};original={p['id']:p for p in atlas['parts']};assert {k for k in original if original[k]!=emitted[k]}==set(changed)
    assert {k: v for k,v in atlas.items()if k not in ['parts','chunks']}=={k:v for k,v in loaded.items()if k not in ['parts','chunks']}
    for side in ['l','r']:
        saved=np.load(BASE/f'bilateral-course-proof-v2/{side}-course.npz');owner=reader.read_part(loaded,emitted['ZA:Lateral circumflex femoral artery.'+side],frame=str(saved['frame']));assert np.array_equal(owner.triangles,saved['vertices'][saved['faces']]);assert np.array_equal(owner.attributes['normals'][owner.faces],saved['normal_words'][saved['faces']])
        child_id='TA98NATIVE:male:LCFA:ascending:'+side;before=reader.read_part(atlas,original[child_id],frame=str(saved['frame']));after=reader.read_part(loaded,emitted[child_id],frame=str(saved['frame']));assert before.geometry_sha==after.geometry_sha and before.attributes_sha==after.attributes_sha
    inventory=json.loads((BASE/'source-ownership-proposal-v1/source-bound-inventory.json').read_bytes());inventory['source_catalogue']=str(OUT/'source-atlas.json');inventory['source_catalogue_sha256']=EXPECTED
    for term in inventory['terms']:
        term['missing_distinct_representation']=term['authority']['id']=='A12.2.16.030'
        if term['authority']['id']=='A12.2.16.029':term['source_meshes']=[dict(part_id=p['id'],source_catalogue=str(OUT/'source-atlas.json'))for p in atlas['parts']if 'Descending' in p.get('name','')and 'circumflex femoral' in p.get('name','').lower()]
        else:
            term['source_meshes']=[dict(part_id=k,source_catalogue=str(OUT/'source-atlas.json'))for k in (['TA98NATIVE:male:LCFA:ascending:'+s for s in ['l','r']] if term['authority']['id']=='A12.2.16.028' else ['ZA:Lateral circumflex femoral artery.'+s for s in ['l','r']])]
        term['planned_representation']['scope']='Saved coarse-course bone clearance correction; original outer surfaces and source interfaces retained. Finite wall and communicating lumen unverified.'
    assert all(t['source_meshes'] for t in inventory['terms'] if t['authority']['id']!='A12.2.16.030')
    atomic_json(OUT/'source-bound-inventory.json',inventory)
    proof.update(passed=True,bank=str(bankpath),bank_sha256=sha,bank_bytes=len(blob),proposed_catalogue=str(OUT/'proposed-atlas.json'),proposed_catalogue_sha256=digest_file(OUT/'proposed-atlas.json'),changed_parts=sorted(changed),elapsed_seconds=time.monotonic()-began,no_unrelated_catalogue_changes=True)
    atomic_json(OUT/'packing-proof.json',proof);print(json.dumps(proof))

if __name__=='__main__':main()
