"""Pack one source-bound cystic/parent pair and all affected named views."""
from pathlib import Path
from collections import defaultdict
import sys,json,time,copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_bytes,atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.tissue import normal_words,decode_packed,validate_cells
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

WORK=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-cystic-blood-donor-v1'
FRAME='model:local-ta98:assembled-metres'


def main():
    began=time.monotonic();role=sys.argv[1];assert role in ['artery','vein']
    source=WORK/'owned-parent-cells-v18'/role;out=WORK/'owned-region-banks-v19'/role;assert not out.exists();out.mkdir(parents=True)
    z=np.load(source/'partition.npz');proof=json.loads((source/'proof.json').read_text());old=json.loads((source/'source-records.json').read_text())
    positions=z['positions'];pf=z['parent_faces'];bf=z['branch_faces'];pids=z['parent_source_face_ids'];bids=z['branch_source_face_ids'];owners=z['parent_owner_ids'];parent_id=proof['parentId'];branch_id=proof['branchId'];oldparent=old['records'][parent_id];origin=np.array(oldparent['geometryOrigin'])
    surfaces={};cells={'parent':[],'cystic':[]};owner_for_surface={};native_for_surface={}
    shared=z['shared_parent_faces'];shared_branch=z['shared_branch_faces'];take_parent=np.ones(len(pf),bool);take_parent[shared]=False;take_branch=np.ones(len(bf),bool);take_branch[shared_branch]=False
    for owner in sorted(set(owners[take_parent])):
        for sheet in sorted(set((pids//1000000)[take_parent&(owners==owner)])):
            selected=np.flatnonzero(take_parent&(owners==owner)&(pids//1000000==sheet));sid='parent:'+owner+':sheet:'+str(sheet);surfaces[sid]=pf[selected];cells['parent'].append(dict(surface=sid,side=1));owner_for_surface[sid]=owner;native_for_surface[sid]=selected
    for sheet in sorted(set((bids//1000000)[take_branch])):
        selected=np.flatnonzero(take_branch&(bids//1000000==sheet));sid='cystic:sheet:'+str(sheet);surfaces[sid]=bf[selected];cells['cystic'].append(dict(surface=sid,side=1))
    sid='reciprocal-parent-cystic-annulus';surfaces[sid]=pf[shared];cells['parent'].append(dict(surface=sid,side=1,neighbor='cystic'));cells['cystic'].append(dict(surface=sid,side=-1,neighbor='parent'));owner_for_surface[sid]=proof['parentRecipientRegion']
    encoded=(positions-origin).astype('<f4').astype(float);words=normal_words(encoded,surfaces)
    literal_native_corner_count=0
    for sid,selected in native_for_surface.items():
        native=np.isin(pids[selected]//1000000,[0,2]);words[sid][native]=z['parent_native_corner_words'][selected[native]];literal_native_corner_count+=int(native.sum())*3
    packed=pack_finite_cells(positions,surfaces,words,cells,{'parent':parent_id,'cystic':branch_id},origin)
    emitted=decode_packed(packed,FRAME,origin)
    # Real occupied-volume checks supplement the exact external-face replay.
    ps,_=make_solid(emitted['parent'].vertices,emitted['parent'].faces);bs,_=make_solid(emitted['cystic'].vertices,emitted['cystic'].faces)
    overlap=float((ps^bs).volume());extent=np.ptp(np.r_[emitted['parent'].vertices,emitted['cystic'].vertices],axis=0).max();bound=2048*np.finfo(float).eps*extent**3*(len(pf)+len(bf));assert overlap<=bound
    checks=validate_cells(emitted,max_pairs=1000000,max_tests=12000000,checkpoint_directory=out/'validation-checkpoints');atomic_json(out/'checks.json',checks);assert checks['passed']
    data=bytearray(packed.data);records={};parent=copy.deepcopy(oldparent)
    for key in ['geometryScale','geometryInstance','sourceOffset','uvs','colors']:parent.pop(key,None)
    parent.update(packed.cells['parent']);parent['pairedTissue']=True;parent['sectionSolid']=True
    for draw,ref in zip(parent['tissueDrawGroups'],cells['parent']):
        if not draw.get('neighbor'):draw['replacement']=owner_for_surface[ref['surface']]
    records[parent_id]=parent
    new=copy.deepcopy(oldparent)
    for key in ['geometryScale','geometryInstance','sourceOffset','uvs','colors','hepaticRegionOf','hepaticAttributeSignature','hepaticConnectionRole']:new.pop(key,None)
    new.update(packed.cells['cystic']);term=proof['newTerm'];new.update(id=branch_id,conceptId=branch_id,name='Cystic artery'if role=='artery'else'Cystic vein',ta98Term=term,sourceId='HRA:VH_F_cystic_'+role,pairedTissue=True,sectionSolid=True,sectionAssembly=oldparent['sectionAssembly'],tissueWallAssembly=oldparent['tissueWallAssembly'],defaultHidden=False,license='HRA united-female v1.10, CC BY 4.0; native parent credits retained',provenance=dict(label='HRA source-specific cystic vessel, local male placement',url='https://humanatlas.io/',detail='Original 18-control proximal profile and separately transported source sections; actual native right-hepatic/right-portal blood-space junction. Finite wall thickness is illustrative.'),description='Macroscopic source-specific cystic blood vessel, connected to native parent blood space and gallbladder wall. Wall thickness and microscopic terminal bed are illustrative.',coverageLimitation='Local review geometry; microscopic perfusion and exact wall thickness are not certified.')
    records[branch_id]=new
    indexed_positions=packed.positions.astype(float)+origin
    for owner in proof['affectedRegionIds']:
        selected=[packed.surfaces[sid]for sid in owner_for_surface if owner_for_surface[sid]==owner]
        fs=np.concatenate(selected);data.extend(b'\0'*((-len(data))%4));offset=len(data);data.extend(fs.astype('<u4').tobytes());record=copy.deepcopy(old['records'][owner]);record.update({k:parent[k]for k in ['positions','normals','normalSigns','vertexCount','geometryOrigin']});record.update(indices=offset,indexCount=int(fs.size),triangles=len(fs),bounds=[indexed_positions[np.unique(fs)].min(0).tolist(),indexed_positions[np.unique(fs)].max(0).tolist()]);groups=[];start=0
        for sid in owner_for_surface:
            if owner_for_surface[sid]!=owner:continue
            count=int(packed.surfaces[sid].size);groups.append(dict(start=start,count=count,**({'neighbor':branch_id}if sid=='reciprocal-parent-cystic-annulus'else{})));start+=count
        record['tissueDrawGroups']=groups;records[owner]=record
    bank=out/'cells-and-regions.bin';atomic_bytes(bank,bytes(data))
    for name,mesh in emitted.items():mesh.save(out/(name+'-emitted.npz'))
    atomic_json(out/'records.json',records)
    manifest=dict(schema='cystic-owned-region-bank/v19',role=role,bank=str(bank),bankSha256=digest_file(bank),bytes=len(data),records=list(records),parentId=parent_id,branchId=branch_id,newTerm=term,sourceRecords=str(source/'source-records.json'),sourceCatalogueSha256=old['catalogue_sha256'],sourceParentBankSha256=old['parent_bank_sha256'],positionOrigin=origin.tolist(),literalNativeSourceNormalCorners=literal_native_corner_count,allAffectedNamedViewsRebuilt=True,reciprocalParentCysticAnnulusTriangles=len(shared),occupiedOverlap_m3=overlap,numericalOverlapBound_m3=float(bound),emittedGeometryChecksPassed=True,parentGeometrySha256=emitted['parent'].geometry_sha,cysticGeometrySha256=emitted['cystic'].geometry_sha,parentPartitionProof=str(source/'proof.json'),checks=str(out/'checks.json'),atlasModified=False,installed=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(out/'manifest.json',manifest);print(json.dumps({k:manifest[k]for k in ['role','bytes','records','occupiedOverlap_m3','emittedGeometryChecksPassed','elapsedSeconds']}))


if __name__=='__main__':main()
