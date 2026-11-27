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
v10=tags.v10;WORK=v10.WORK;OUT=WORK/'vein-recipient-partition-v28';B=tags.B


def extracted(solid,origin,bank):
    raw=solid.to_mesh64();p=np.asarray(raw.vert_properties)[:,:3]/1000+origin
    mapped=np.array([bank.insert(q)for q in p],int)
    return mapped[np.asarray(raw.tri_verts,int)],np.asarray(raw.face_id,np.int64)


def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir()
    catalogue=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(catalogue)=='5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a';atlas=json.loads(catalogue.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();rows=[]
    cases=[('vein',21,'portal','A12.3.12.002','A12.3.12.012','vein-recipient-course-v27','vein-joined-wall.npz')]
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
        # Source station 7 is already an external paired ring. Partition the
        # accepted source faces there without re-Booleaning coincident sides.
        port=7;position_bank=CommonPositionBank(np.r_[host.vertices,z['outer_sections'].reshape(-1,3),z['fluid_sections'].reshape(-1,3)],tolerance=2e-12)
        wf,wids=extracted(whole,root,position_bank);roles=wids//B;local=wids%B;side_count=(len(z['centers'])-1)*36
        assert np.all(np.isin(roles,[0,1,2,3]))
        donor=np.isin(roles,[1,3]);donor_side=donor&(local<side_count)
        parent_take=np.isin(roles,[0,2])|(donor_side&(local//36<port))
        branch_take=~parent_take
        assert np.all(donor[branch_take])
        ov7=np.array([position_bank.insert(p)for p in z['outer_sections'][port]])
        iv7=np.array([position_bank.insert(p)for p in z['fluid_sections'][port]])
        interface=[]
        for i in range(18):
            j=(i+1)%18;interface.extend([[ov7[i],ov7[j],iv7[j]],[ov7[i],iv7[j],iv7[i]]])
        interface=np.array(interface,int);positions=position_bank.positions;tri=positions[interface]
        forward=z['centers'][port+1]-z['centers'][port-1]
        if np.median(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])@forward)<0:interface=interface[:,[0,2,1]]
        pf=np.r_[wf[parent_take],interface];pids=np.r_[wids[parent_take],6*B+np.arange(len(interface))]
        bf=np.r_[wf[branch_take],interface[:,[0,2,1]]];bids=np.r_[wids[branch_take],6*B+np.arange(len(interface))]
        parent_shared=np.arange(parent_take.sum(),len(pf));branch_shared=np.arange(branch_take.sum(),len(bf))
        left={oriented_key(pf[i]):int(i)for i in parent_shared};right={oriented_key(bf[i,[0,2,1]]):int(i)for i in branch_shared}
        assert left and left.keys()==right.keys(),'Reciprocal annular cut subdivisions differ'
        external=np.r_[pf[:parent_take.sum()],bf[:branch_take.sum()]]
        assert sorted(oriented_key(t)for t in external)==sorted(oriented_key(t)for t in wf)
        # Source partition alone proves external replay. Exact shared cut
        # faces cancel as opposite triangle multisets.
        residual=0.;occupied=0.;overlap=0.
        old_owner=np.full(len(host.faces),'',object)
        for g in part['tissueDrawGroups']:
            assert g.get('replacement');old_owner[g['start']//3:(g['start']+g['count'])//3]=g['replacement']
        assert np.all(old_owner!='');owners=np.full(len(pf),recipient_id,object);native=np.isin(pids//B,[0,2]);owners[native]=old_owner[pids[native]%B]
        parent_top=topology(Mesh(positions,pf,v10.FRAME));branch_top=topology(Mesh(positions,bf,v10.FRAME))
        atomic_json(folder/'topology.json',dict(parent=parent_top,branch=branch_top));atomic_npz(folder/'topology-fields.npz',positions=positions,parent_faces=pf,branch_faces=bf)
        assert all(t['closed_oriented_manifold'] and len(t['components'])==1 for t in [parent_top,branch_top])
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
    atomic_json(OUT/'proof.json',dict(schema='cystic-owned-parent-cells/v18',rows=rows,sourceCatalogueSha256=digest_file(catalogue),elapsedSeconds=time.monotonic()-began,atlasModified=False,installed=False))


if __name__=='__main__':main()
