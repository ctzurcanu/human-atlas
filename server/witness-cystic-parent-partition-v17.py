"""Partition passing cystic wall candidates with explicit native face ancestry.

Writes local construction evidence only. It neither edits an atlas nor installs.
"""
from pathlib import Path
from collections import defaultdict
import sys,json,time,importlib.util
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.audit import topology
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank,oriented_key,barycentric
sp=importlib.util.spec_from_file_location('tags',ROOT/'server/cystic-blood-handle-witness-v11.py');tags=importlib.util.module_from_spec(sp);sp.loader.exec_module(tags)
v10=tags.v10;WORK=v10.WORK;OUT=WORK/'owned-parent-cells-v17-topology-witness';B=tags.B


def extracted(solid,origin,bank):
    raw=solid.to_mesh64();p=np.asarray(raw.vert_properties)[:,:3]/1000+origin
    mapped=np.array([bank.insert(q)for q in p],int)
    return mapped[np.asarray(raw.tri_verts,int)],np.asarray(raw.face_id,np.int64)


def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir()
    catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    atlas=json.loads(catalogue.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();rows=[]
    cases=[('artery',20,'artery','A12.2.12.030','A12.2.12.031','admissible-departure-v14','artery-joined-wall.npz'),
           ('vein',21,'portal','A12.3.12.002','A12.3.12.012','cystic-vein-v16','vein-joined-wall.npz')]
    for role,number,parent_role,recipient_term,term,stage,reference in cases:
        saved=json.loads((WORK/stage/'gates.json').read_text());assert saved['allRequiredGeometryGatesPassed']
        parent_id='TA98HEPATIC:male:'+parent_role;recipient_id='TA98HEPATIC:male:region:'+recipient_term;branch_id='TA98CYSTICBLOOD:male:'+role
        part=parts[parent_id];host=reader.read_part(atlas,part,frame=v10.FRAME,catalogue_path=catalogue)
        old=Mesh.load(WORK/('fixed-target-'+str(number)+'.npz'));assert host.geometry_sha==old.geometry_sha,'Current parent is no longer the evaluated native source'
        folder=OUT/role;folder.mkdir();z=np.load(WORK/stage/'sections.npz');root=z['centers'][0];blood=Mesh.load(WORK/'native-parent-envelopes-v9'/(role+'-blood.npz'))
        lineage=json.loads((WORK/'native-parent-envelopes-v9'/(role+'-inner-lineage.json')).read_text());lookup={tuple(p):i for i,p in enumerate(blood.vertices)}
        host_faces=np.array([[lookup[tuple(p)]for p in tri]for tri in host.triangles]);inner=np.array(lineage['sourceFaceIds']);outer=np.setdiff1d(np.arange(len(host.faces)),inner);nd=lineage['temporaryRootDiskFaces'];disk=blood.faces[-nd:]
        ho=tags.tagged(blood.vertices,np.r_[host_faces[outer],disk],np.r_[outer,4*B+np.arange(nd)],root)
        hb=tags.tagged(blood.vertices,blood.faces,np.r_[2*B+inner,5*B+np.arange(nd)],root)
        ov,of,oc=v10.swept_profiles(z['outer_sections']);iv,iff,ic=v10.swept_profiles(z['fluid_sections'])
        bo=tags.tagged(ov,of,B+np.arange(len(of)),root);bb=tags.tagged(iv,iff,3*B+np.arange(len(iff)),root)
        whole=(ho+bo)-(hb+bb);original=Mesh.load(WORK/stage/reference);baseline=v10.previous.solid(original.vertices,original.faces,root)
        symmetric=float(((whole-baseline)+(baseline-whole)).volume())/1e9;assert abs(symmetric)<1e-20,'Tagged replay differs physically from passing candidate'
        # Cut at the already constructed external donor port. This cut is an
        # annular tissue interface, not a disk across the blood space.
        port=7;cv,cf,cc=v10.swept_profiles(z['outer_sections'][port:]);side_count=len(cf)-36
        cut_ids=np.r_[B+np.arange(port*36,port*36+side_count),6*B+np.arange(18),B+oc[1]]
        cutter=tags.tagged(cv,cf,cut_ids,root)
        overlap=float(((ho-hb)^cutter).volume())/1e9
        assert abs(overlap)<1e-20,'External donor partition would take native parent tissue'
        parent=whole-cutter;branch=whole^cutter
        residual=float(((parent+branch)-whole).volume()+((whole-(parent+branch)).volume()))/1e9;assert abs(residual)<1e-20
        occupied=float((parent^branch).volume())/1e9;assert abs(occupied)<1e-20
        position_bank=CommonPositionBank(np.r_[host.vertices,z['outer_sections'].reshape(-1,3),z['fluid_sections'].reshape(-1,3)],tolerance=2e-12)
        pf,pids=extracted(parent,root,position_bank);bf,bids=extracted(branch,root,position_bank);positions=position_bank.positions
        assert not np.any(np.isin(pids//B,[4,5])) and not np.any(np.isin(bids//B,[0,2,4,5]))
        parent_shared=np.flatnonzero(pids//B==6);branch_shared=np.flatnonzero(bids//B==6)
        left={oriented_key(pf[i]):int(i)for i in parent_shared};right={oriented_key(bf[i,[0,2,1]]):int(i)for i in branch_shared}
        assert left and left.keys()==right.keys(),'Reciprocal annular cut subdivisions differ'
        old_owner=np.full(len(host.faces),'',object)
        for g in part['tissueDrawGroups']:
            assert g.get('replacement');old_owner[g['start']//3:(g['start']+g['count'])//3]=g['replacement']
        assert np.all(old_owner!='');owners=np.full(len(pf),recipient_id,object);native=np.isin(pids//B,[0,2]);owners[native]=old_owner[pids[native]%B]
        parent_top=topology(Mesh(positions,pf,v10.FRAME));branch_top=topology(Mesh(positions,bf,v10.FRAME))
        atomic_npz(folder/'failed-topology.npz',positions=positions,parent_faces=pf,parent_source_face_ids=pids,branch_faces=bf,branch_source_face_ids=bids)
        atomic_json(folder/'failed-topology.json',dict(parent=parent_top,branch=branch_top,parentManifoldComponentVolumes_m3=[float(c.volume())/1e9 for c in parent.decompose()],branchManifoldComponentVolumes_m3=[float(c.volume())/1e9 for c in branch.decompose()],maximumSnap_m=position_bank.maximum_snap,symmetricDifference_m3=symmetric,sharedTriangles=len(parent_shared)))
        print(json.dumps(dict(parent=parent_top,branch=branch_top,parentComponentVolumes=[float(c.volume())/1e9 for c in parent.decompose()],branchComponentVolumes=[float(c.volume())/1e9 for c in branch.decompose()],maximumSnap=position_bank.maximum_snap)),flush=True)
        return
        # Preserve every literal native corner normal word; new intersection
        # controls inherit the actual native parent chart barycentrically.
        source_words=host.attributes['normals'].astype(np.int64)*host.attributes['normalSigns'][:,None]
        nwords=np.zeros((len(pf),3,3),'<i2');literal_faces=0;split_faces=0;changed_source=[]
        for i in np.flatnonzero(native):
            source_face=int(pids[i]%B);tri=host.triangles[source_face];word=source_words[host.faces[source_face]];q=positions[pf[i]];d=np.linalg.norm(q[:,None]-tri[None,:],axis=2);nearest=d.argmin(1)
            exact=d.min(1)<2e-12
            if exact.all() and len(set(map(int,nearest)))==3:literal_faces+=1
            else:split_faces+=1;changed_source.append(source_face)
            interpolation=barycentric(q,tri)@word;length=np.linalg.norm(interpolation,axis=1);interpolation=np.rint(interpolation/length[:,None]*32767).astype('<i2');interpolation[exact]=word[nearest[exact]].astype('<i2');nwords[i]=interpolation
        # The untouched native triangles and their normal fields are replayed
        # literally; the patch is the only allowed changed native ancestry.
        atomic_npz(folder/'partition.npz',positions=positions,parent_faces=pf,parent_source_face_ids=pids,parent_owner_ids=np.array(owners,str),branch_faces=bf,branch_source_face_ids=bids,parent_native_corner_words=nwords,shared_parent_faces=parent_shared,shared_branch_faces=branch_shared)
        Mesh(positions,pf,v10.FRAME).save(folder/'parent-wall.npz');Mesh(positions,bf,v10.FRAME).save(folder/'cystic-wall.npz')
        selected={pid:parts[pid]for pid in set(old_owner)};selected[parent_id]=part
        atomic_json(folder/'source-records.json',dict(catalogue_sha256=digest_file(catalogue),records=selected,record_sha256={pid:identity(record)for pid,record in selected.items()},parent_geometry_sha256=host.geometry_sha,parent_bank_sha256=host.source['bank_sha256']))
        proof=dict(role=role,parentId=parent_id,branchId=branch_id,newTerm=term,parentRecipientRegion=recipient_id,sourceCatalogueSha256=digest_file(catalogue),passingCandidateStage=stage,passingCandidateReplayedPhysicalDifference_m3=symmetric,nativeParentTissueTakenByBranch_m3=overlap,parentBranchOccupiedOverlap_m3=occupied,wholeReplayDifference_m3=residual,reciprocalInterfaceTriangles=len(parent_shared),reciprocalIndicesAndCoordinates=True,newInterfaceIsAnnulusNotBloodDisk=True,exposedComputationalRootDisks=0,parentTopology=parent_top,branchTopology=branch_top,nativeAncestryOutputFaces=int(native.sum()),literalNativeFaces=literal_faces,splitNativeFaces=split_faces,splitNativeSourceFaces=sorted(set(changed_source)),allAffectedRegionalViewsRebuildFromAncestry=True,affectedRegionIds=sorted(set(old_owner)),maximumIntersectionCanonicalization_m=position_bank.maximum_snap,atlasModified=False,installed=False,anatomicalAcceptance=False)
        atomic_json(folder/'proof.json',proof);rows.append({k:proof[k]for k in ['role','parentId','branchId','affectedRegionIds','reciprocalInterfaceTriangles','literalNativeFaces','splitNativeFaces']});print(json.dumps(rows[-1]),flush=True)
    atomic_json(OUT/'proof.json',dict(schema='cystic-owned-parent-cells/v17',rows=rows,sourceCatalogueSha256=digest_file(catalogue),elapsedSeconds=time.monotonic()-began,atlasModified=False,installed=False))


if __name__=='__main__':main()
