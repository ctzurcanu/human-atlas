"""Emit the bounded14-control female core-only partial; no fitting or install."""
import copy,json,sys,time,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'server')]
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from tongue_basal_contract import require
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
SOURCE=W/'female-partial-core-native-emitted-v1';OUT=W/'female-core-only-emitted-v1'
FRAME='model:local-female-ta98:assembled-metres'
NATIVE=['LOCAL:female:genioglossus:239','LOCAL:female:l_hyoglossus:242','LOCAL:female:r_hyoglossus:255','LOCAL:female:l_styloglossus:243','LOCAL:female:r_styloglossus:256']
LIMIT=('Partial female tongue core:14 of35 ordinary-adult TA98 controls;21 adult terms remain missing. '
       'Male-derived mucosa/aponeurosis and approximate intrinsic carriers are installed as a partial preview. '
       'Five original female extrinsic muscle sources remain preserved but suppressed as unfinished work; '
       'their placement/attachments and their terms are not represented by this new assembly. '
       'Intrinsic regions are non-owning macroscopic envelopes, not exclusive histological tissues.')
def bind(p):return dict(path=str(Path(p).resolve()),sha256=digest_file(p))
def main():
    started=time.monotonic();prior=SOURCE/'candidate-atlas.json';require(digest_file(prior)=='0aa5e329dc0441eb2bd6c01d80ef3525a14a72dde5d83007c3eac846fb5ef872','Frozen six-part source packet changed')
    old=json.loads(prior.read_text());e=json.loads((SOURCE/'transfer-evidence.json').read_text());a=copy.deepcopy(old);m=a['finiteTongueIntegration'];new=e['new_part_ids'];require(len(new)==6 and not OUT.exists(),'Use the exact six new parts and a new stage')
    drop={'genioglossus','hyoglossus','styloglossus'};removed={g['conceptId'] for g in m['termGroups'] if g['key'] in drop}
    m['termGroups']=[g for g in m['termGroups'] if g['key'] not in drop];m['partIds']=new;m['activePhysicalMeshes']=6
    for g in m['termGroups']:g['partIds']=[p for p in g['partIds'] if p in new]
    require(len(m['termGroups'])==14 and all(g['partIds'] for g in m['termGroups']),'Core controls must retain actual carriers')
    m.update(partialMode='female-core-only',representedAdultControls=14,missingAdultTerms=21,
        unfinishedExtrinsicSourcesSuppressed=True,suppressedRetainedWorkPartIds=NATIVE,
        coverageLimitation=LIMIT,anatomicalAcceptance=False,tissueOwnershipClaim=False)
    m.pop('retainedNativePartIds',None);m.pop('partialNativeExtrinsicsRetained',None)
    m['carrierBankHashes']={p:h for p,h in m['carrierBankHashes'].items() if p in new}
    m['transfer'].update(new_physical_meshes=6,retained_native_meshes=0,suppressed_unfinished_native_meshes=5,
        qualifying_candidate=False,actual_render_review_required=True,target_context_checks_pending=True)
    a['concepts']=[c for c in a['concepts'] if c['id'] not in removed]
    groups={g['conceptId']:g for g in m['termGroups']}
    aliases={'LOCAL:female:tongue':new,'TA98FEMALE:A05.1.04.001':new,
             'HRA:VH_F_dorsal_tongue':['TA98TONGUE:female:mucosa'],'TA98FEMALE:A05.1.04.004':['TA98TONGUE:female:mucosa']}
    for c in a['concepts']:
        if c['id'] in groups:c.update(elements=groups[c['id']]['partIds'],coverageLimitation=LIMIT)
        elif c['id'] in aliases:c.update(elements=aliases[c['id']],coverageLimitation=LIMIT)
    for p in a['parts']:
        if p['id'] in NATIVE:p['suppressed']=True
        if p['id'] in new:p['coverageLimitation']=LIMIT
    a['defaultHidden']=list(dict.fromkeys(a.get('defaultHidden',[])+NATIVE));a['triangles']=sum(p['indexCount']//3 for p in a['parts'] if not p.get('suppressed'))
    OUT.mkdir();(OUT/'meshes').mkdir();reader=AtlasReader();parts={p['id']:p for p in a['parts']};records=[]
    for row in e['actual_world_meshes']:
        if row['part_id'] not in new:continue
        actual=reader.read_part(a,parts[row['part_id']],frame=FRAME);previous=Mesh.load(row['mesh']['path'])
        require(np.array_equal(actual.vertices,previous.vertices) and np.array_equal(actual.faces,previous.faces)
                and actual.attributes_sha==previous.attributes_sha,'Final scope change altered checked geometry or normal fields')
        path=OUT/'meshes'/Path(row['mesh']['path']).name;shutil.copyfile(row['mesh']['path'],path)
        require(digest_file(path)==row['mesh']['sha256'],'Passed decoded NPZ must remain byte-identical')
        records.append(dict(row,mesh=bind(path),array_and_attribute_replay_literal=True))
    out=OUT/'candidate-atlas.json';atomic_json(out,a)
    evidence=dict(schema='tongue-female-core-only-partial/v1',candidate=bind(out),source_packet=bind(prior),source=copy.deepcopy(e['source']),target=copy.deepcopy(e['target']),
        source_evidence=bind(SOURCE/'transfer-evidence.json'),new_part_ids=new,active_part_ids=new,suppressed_retained_work_part_ids=NATIVE,
        suppressed_source_parts=e['suppressed_source_parts']+NATIVE,legacy_alias_elements=aliases,actual_world_meshes=records,fixed_context_meshes=e['fixed_context_meshes'],
        source_bank_bindings=e['source_bank_bindings'],IL_source_proof=e['IL_source_proof'],controls=14,adult_terms_missing=21,physical_meshes=6,
        fixed_skeleton=True,all_six_decoded_arrays_byte_identical=True,qualifying_candidate=False,anatomical_acceptance=False,installed=False)
    atomic_json(OUT/'transfer-evidence.json',evidence)
    meshes=[]
    for row in records:
        p=parts[row['part_id']];color=a['materials'][p['material']]['color'];color=[v/255 for v in color[:3]] if max(color)>1 else color[:3]
        item=dict(mesh=row['mesh']['path'],name=p['name'],color=color+[1.],part_id=p['id'])
        if p.get('tissueDrawGroups'):item['tissue_draw_groups']=copy.deepcopy(p['tissueDrawGroups'])
        meshes.append(item)
    for pid in ['LOCAL:female:mandible:960','LOCAL:female:hyoid:1002']:
        row=next(r for r in evidence['fixed_context_meshes'] if r['part_id']==pid);meshes.append(dict(mesh=row['mesh']['path'],name=parts[pid]['name'],color=[.82,.80,.72,.38]))
    common=dict(frame=FRAME,assembled_tissues=True,meshes=meshes,resolution=1000,wall_seconds=45)
    atomic_json(OUT/'review-surface-right.json',dict(common,views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[1,0,0])]))
    muco=reader.read_part(a,parts['TA98TONGUE:female:mucosa'],frame=FRAME);point=muco.bounds.mean(0).tolist()
    atomic_json(OUT/'review-sagittal.json',dict(common,views=[dict(name='sagittal',direction=[1,0,0],section=dict(normal=[1,0,0],point=point,clear_outer=True))]))
    print(json.dumps(dict(candidate=str(out),sha256=digest_file(out),controls=14,meshes=6,adult_missing=21,
        surface_review=str(OUT/'review-surface-right.json'),section_review=str(OUT/'review-sagittal.json'),elapsed_seconds=time.monotonic()-started)))
if __name__=='__main__':main()
