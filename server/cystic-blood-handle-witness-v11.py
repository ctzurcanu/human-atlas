"""Trace actual composed blood handles and native/donor junction curves.

Replays the saved v10 geometry exactly; it never changes a course or seals a
terminal. Face IDs distinguish every native sheet and donor classifier face.
"""
from pathlib import Path
import sys,json,time,importlib.util
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
import manifold3d as m
spec=importlib.util.spec_from_file_location('v10',ROOT/'server/correct-cystic-profiles-and-host-mouths-v10.py');v10=importlib.util.module_from_spec(spec);spec.loader.exec_module(v10)
WORK=v10.WORK;OUT=WORK/'blood-handle-witness-v11';B=1000000

def tagged(v,f,ids,origin):
    used=np.unique(f);points,inv=np.unique(v[used],axis=0,return_inverse=True);mapping=np.full(len(v),-1,int);mapping[used]=inv
    solid=m.Manifold(m.Mesh64((points-origin)*1000,mapping[f].astype(np.uint64),face_id=np.array(ids,np.uint64)))
    assert solid.status()==m.Error.NoError;return solid

def save(solid,name,origin):
    packed=solid.to_mesh64();points=np.asarray(packed.vert_properties)[:,:3]/1000+origin;faces=np.asarray(packed.tri_verts,np.int64);ids=np.asarray(packed.face_id,np.uint64)
    mesh=Mesh(points,faces,v10.FRAME);mesh.save(OUT/(name+'.npz'));atomic_npz(OUT/(name+'-source-face-ids.npz'),source_face_ids=ids)
    top=topology(mesh);return mesh,ids,dict(components=len(top['components']),genera=[c['closed_orientable_genus']for c in top['components']],euler=[c['euler']for c in top['components']],volume_m3=float(solid.volume())/1e9,closed=top['closed_oriented_manifold'])

def interface_curves(mesh,face_ids,first,second):
    v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];edges=defaultdict(list)
    for fi,tri in enumerate(f):
        for a,b in zip(tri,np.roll(tri,-1)):edges[tuple(sorted((int(a),int(b))))].append(fi)
    selected=[];provenance=[]
    for edge,owners in edges.items():
        if len(owners)==2 and {int(face_ids[i])//B for i in owners}=={first,second}:selected.append(edge);provenance.append([int(face_ids[i])for i in owners])
    adj=defaultdict(set)
    for a,b in selected:adj[a].add(b);adj[b].add(a)
    groups=[];seen=set()
    for vertex in sorted(adj):
        if vertex in seen:continue
        stack=[vertex];members=[]
        while stack:
            node=stack.pop()
            if node in seen:continue
            seen.add(node);members.append(node);stack.extend(adj[node]-seen)
        members=set(members);ids=[i for i,e in enumerate(selected)if e[0]in members]
        coords=v[sorted(members)];groups.append(dict(vertices=len(members),edges=len(ids),simpleClosedCurve=all(len(adj[i])==2 for i in members),bounds_m=[coords.min(0).tolist(),coords.max(0).tolist()],center_m=coords.mean(0).tolist(),sourceFaceIds=sorted({pid for i in ids for pid in provenance[i]})))
    return dict(curves=groups,segments=np.array(selected,int).reshape(-1,2),points=v,edgeSourceFaceIds=np.array(provenance,np.uint64).reshape(-1,2))

def main():
    began=time.monotonic();assert not OUT.exists();OUT.mkdir();role='artery';data=np.load(WORK/'profile-mouth-correction-v10/artery-sections.npz');centers=data['centers'];root=centers[0]
    host=Mesh.load(WORK/'fixed-target-20.npz');blood=Mesh.load(WORK/'native-parent-envelopes-v9/artery-blood.npz');lineage=json.loads((WORK/'native-parent-envelopes-v9/artery-inner-lineage.json').read_text())
    lookup={tuple(p):i for i,p in enumerate(blood.vertices)};f=np.array([[lookup[tuple(p)]for p in tri]for tri in host.triangles]);innerids=np.array(lineage['sourceFaceIds']);outerids=np.setdiff1d(np.arange(len(f)),innerids);nd=lineage['temporaryRootDiskFaces'];disk=blood.faces[-nd:]
    ho=tagged(blood.vertices,np.r_[f[outerids],disk],np.r_[outerids,4*B+np.arange(nd)],root)
    hb=tagged(blood.vertices,blood.faces,np.r_[2*B+innerids,5*B+np.arange(nd)],root)
    ov,of,oc=v10.swept_profiles(data['outer_sections']);iv,iff,ic=v10.swept_profiles(data['fluid_sections']);bo=tagged(ov,of,B+np.arange(len(of)),root);bb=tagged(iv,iff,3*B+np.arange(len(iff)),root)
    rows={}
    for name,solid in [('host-outer',ho),('host-blood',hb),('donor-outer',bo),('donor-blood',bb),('composed-outer',ho+bo),('composed-blood',hb+bb)]:
        _,_,rows[name]=save(solid,name,root)
    overlap=hb^bb
    components=overlap.decompose();overlaps=[]
    for i,component in enumerate(components):
        mesh,ids,top=save(component,'blood-overlap-'+str(i),root);overlaps.append(dict(index=i,topology=top,bounds_m=mesh.bounds.tolist(),center_m=mesh.vertices.mean(0).tolist(),sourceFacesByRole={str(r):sorted(set(map(int,ids[ids//B==r]%B)))for r in np.unique(ids//B)}))
    joined=(ho+bo)-(hb+bb);mesh,ids,rows['composed-wall']=save(joined,'composed-wall',root)
    original=Mesh.load(WORK/'profile-mouth-correction-v10/artery-joined-wall.npz');assert mesh.surface_sha==original.surface_sha,'Tagged replay changed the saved candidate surface'
    curves={}
    for name,a,b in [('outer-native-donor',0,1),('inner-native-donor',2,3)]:
        curve=interface_curves(mesh,ids,a,b);atomic_npz(OUT/(name+'-curves.npz'),points=curve.pop('points'),segments=curve.pop('segments'),edgeSourceFaceIds=curve.pop('edgeSourceFaceIds'));curves[name]=curve
    native_root=json.loads((WORK/'native-parent-envelopes-v9/proof.json').read_text())['rows'][0]['proof']['annularRootFaceIds']
    source_caps=dict(donorOuterStart=list(map(int,oc[0])),donorOuterTerminal=list(map(int,oc[1])),donorBloodStart=list(map(int,ic[0])),donorBloodTerminal=list(map(int,ic[1])))
    exposed=dict(nativeRootFaces=sorted(set(map(int,ids[(ids//B==0)&np.isin(ids,native_root)]))),donorOuterStartFaces=sorted(set(map(int,ids[(ids//B==1)&np.isin(ids%B,oc[0])]%B))),donorOuterTerminalFaces=sorted(set(map(int,ids[(ids//B==1)&np.isin(ids%B,oc[1])]%B))),donorBloodStartFaces=sorted(set(map(int,ids[(ids//B==3)&np.isin(ids%B,ic[0])]%B))),donorBloodTerminalFaces=sorted(set(map(int,ids[(ids//B==3)&np.isin(ids%B,ic[1])]%B))),temporaryParentClassifierDiskFaces=int(np.sum(np.isin(ids//B,[4,5]))))
    result=dict(schema='cystic-arterial-blood-handle-witness/v11',savedCandidateReplayedExactly=True,betti=dict((k,dict(b0=v['components'],b1=sum(g or 0 for g in v['genera'])))for k,v in rows.items()),topology=rows,bloodOverlapComponents=overlaps,junctionCurves=curves,sourceCaps=source_caps,exposedFaces=exposed,intendedPorts=['One pre-existing native upstream blood mouth','One cystic branch terminal at gallbladder serosal bed'],legitimateDonorDistalTerminals=1,sourceDonorSectionLoopsPerStation=1,geometryChanged=False,installed=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'witness.json',result);print(json.dumps(dict(output=str(OUT/'witness.json'),betti=result['betti'],bloodOverlapComponents=len(overlaps),junctionCurveCounts={k:len(v['curves'])for k,v in curves.items()},exposedFaces=exposed,elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
