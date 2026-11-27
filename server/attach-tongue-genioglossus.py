"""One bone-aware local GG origin candidate after bounded basal untangling.

Only a small native genial-tip window extends toward the witnessed cortex.
Literal fixed mandibular/hyoid surfaces trim that extension and native basal
bone penetration, yielding opposed source-chart contacts instead of a snap
onto a single face. Neither fixed bone bank nor catalogue is modified.
"""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import subtract_occupied
from geometry_tools.algorithms.conforming_fields import replay_corner_fields
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

BASE=ROOT/'.local-models/geometry-workbench/tongue-completion'
SOURCE=BASE/'source-v1/meshes';PREVIOUS=BASE/'genioglossus-v2-stage-a';OUT=BASE/'genioglossus-v2-stage-b'
CATALOGUE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='cbed9740aca9cd20ca12bf43a051eb7b36f60f90dc1919b64818dde8848f963f'
START=time.monotonic();DEADLINE=START+50


def check():
    if time.monotonic()>DEADLINE:raise TimeoutError('GG attachment stage exceeded its 50 second native deadline')
    if hashlib.sha256(CATALOGUE.read_bytes()).hexdigest()!=EXPECTED:raise ValueError('Authoritative source catalogue changed')


def norm_words(v,f):
    t=v[f];n=np.zeros_like(v);cross=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    for i in range(3):np.add.at(n,f[:,i],cross)
    lengths=np.linalg.norm(n,axis=1);n[lengths>0]/=lengths[lengths>0,None]
    return np.rint(np.clip(n,-1,1)*32767).astype('<i2')


def local_mandible_classifier(bone):
    """Literal genial cortex with remote computational backing and side walls."""
    t=bone.triangles;c=t.mean(1);area=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);unit=area/np.linalg.norm(area,axis=1)[:,None]
    witness=np.array([0.,1.623,.077]);center=c[11088].copy();center[0]=0
    # Include both genial origins on a single actual posterior-cortex patch.
    allowed=(np.linalg.norm(c-center,axis=1)<.010)&(unit[:,2]<-.45)
    _,first,alias=np.unique(bone.vertices,axis=0,return_index=True,return_inverse=True)
    ff=alias[bone.faces];edge_faces=defaultdict(list)
    for i,face in enumerate(ff):
        for a,b in zip(face,np.roll(face,-1)):edge_faces[tuple(sorted([int(a),int(b)]))].append(i)
    adj=defaultdict(set)
    for values in edge_faces.values():
        for i in values:adj[i].update(j for j in values if j!=i)
    selected={11088};todo=[11088]
    while todo:
        i=todo.pop()
        for j in adj[i]:
            if j not in selected and allowed[j]:selected.add(j);todo.append(j)
    ids=np.array(sorted(selected),int);loops=boundary_loops(ff[ids])
    if len(loops)!=1:raise ValueError('Witnessed genial cortex is not one local disk; do not invent classifier closure')
    # Extrude only the computational classifier into bone. Its actual cortex
    # and original normals remain literal, while all virtual faces are marked
    # and must prove disjoint from the complete host before subtraction.
    inward=-area[ids].sum(0);inward/=np.linalg.norm(inward)
    v=np.vstack([bone.vertices,bone.vertices+.005*inward]);n=len(bone.vertices)
    cap=bone.faces[ids];inner=cap[:,[0,2,1]]+n;loop=first[loops[0]];side=[]
    for a,b in zip(loop,np.roll(loop,-1)):side.extend([[int(a),int(a+n),int(b+n)],[int(a),int(b+n),int(b)]])
    f=np.vstack([cap,inner,np.array(side,int)]);virtual=np.arange(len(cap),len(f))
    proof=topology(Mesh(v,f,bone.frame))
    if not proof['closed_oriented_manifold']:raise ValueError('Local computational classifier closure is not a valid cell')
    return v,f,np.vstack([bone.attributes['normals']]*2),np.vstack([bone.attributes['uvs']]*2),virtual,dict(
        actual_cortex_source_face_ids=ids.tolist(),actual_cortex_faces_literal=True,
        virtual_faces=virtual.tolist(),temporary_backing_depth_m=.005,topology=proof,
        source_bone_modified=False,virtual_faces_must_not_be_emitted=True)


def main():
    check();OUT.mkdir(exist_ok=True)
    if (OUT/'candidate.json').exists():raise FileExistsError('Keep the one completed attachment candidate immutable')
    before=json.loads((PREVIOUS/'candidate.json').read_text())
    if not all(r['numerical_self_repair_passed']for r in before['candidates']):raise ValueError('Do not extend a rejected self-repair candidate')
    bones=[Mesh.load(SOURCE/f'native-{i:02d}.npz')for i in [5,13,14,15,16,17]]
    mv,mf,mn,mu,virtual,classifier=local_mandible_classifier(bones[0])
    bv=[mv];bf=[mf];bn=[mn];bu=[mu]
    owners=['ZA:Mandible:Bone']*len(mf)
    for i in virtual:owners[i]='COMPUTE:mandible-genial-classifier'
    offset=len(mv)
    for bone in bones[1:]:
        bv.append(bone.vertices);bf.append(bone.faces+offset);bn.append(bone.attributes['normals']);bu.append(bone.attributes.get('uvs',np.zeros((len(bone.vertices),2),'<f4')))
        owners.extend([bone.source['record']['id']]*len(bone.faces));offset+=len(bone.vertices)
    bonev,bonef,bonen,boneuv=np.vstack(bv),np.vstack(bf),np.vstack(bn),np.vstack(bu)
    cutter=Mesh(bonev,bonef,bones[0].frame)
    if not topology(cutter)['closed_oriented_manifold']:raise ValueError('Fixed bone cutter is not a closed assembled compartment')
    original_contacts=json.loads((BASE/'native-attachments-v1/contact-audit.json').read_text())['tests'];rows=[]
    for side,index in [('left',1),('right',2)]:
        check();native=Mesh.load(SOURCE/f'native-{index:02d}.npz');prior=Mesh.load(PREVIOUS/f'genioglossus-{side}.npz')
        witness=next(r['closest']for r in original_contacts if r['musclePartId']==native.source['record']['id']and r['bonePartId']=='ZA:Mandible:Bone')
        center=np.array(witness['first']);direction=np.array(witness['second'])-center;direction/=np.linalg.norm(direction)
        radius=.0045;distance=np.linalg.norm(prior.vertices-center,axis=1);weight=np.maximum(1-(distance/radius)**2,0)**2
        weight[prior.vertices[:,1]>=1.630]=0
        v=prior.vertices+.0009*weight[:,None]*direction
        moved=np.flatnonzero(weight>0);attrs={k:a.copy()for k,a in prior.attributes.items()}
        touched=np.flatnonzero(np.any(np.isin(prior.faces,moved),axis=1));shade=np.unique(prior.faces[touched]);attrs['normals'][shade]=norm_words(v,prior.faces)[shade]
        extension=Mesh(v,prior.faces,prior.frame,attrs,{'record':prior.source['record'],'local_origin_extension':True})
        extension.save(OUT/f'genioglossus-{side}-extended-before-bone-trim.npz')
        prediag=diagnostics(extension,self_intersections=True)
        if prediag['self_intersection_pairs']:raise ValueError(f'{side}: local extension introduced self-crossings; reject without carving')
        surfaces={'GG:outer':prior.faces};surface_rows={'GG:outer':dict(id='GG:outer',kind='outer',array='outer')}
        cells={'GG':prior.faces};parents={'GG':[('GG:outer',i,1)for i in range(len(prior.faces))]};uses={'GG':[dict(surface='GG:outer',side=1)]}
        check()
        result=subtract_occupied(v,surfaces,surface_rows,cells,parents,uses,bonev,bonef,
            classifier_cap_faces=virtual,occupied_face_owners=owners,exact_source_replay=True,reconcile_shared_subdivisions=True,
            diagnostic_directory=OUT/f'{side}-boolean-proof')
        chart_rows={r['id']:r for r in result['surfaceRows']};fields={}
        for sid,chart in result['charts'].items():
            bone_contact=chart_rows[sid].get('parentSurface')=='occupied-outer'
            fields[sid]=dict(source_surfaces=['fixed:bone'if bone_contact else'GG:outer']*len(result['surfaces'][sid]),
                source_faces=chart['parentFaces'],normal_signs=np.full(len(result['surfaces'][sid]),-1 if bone_contact else 1),barycentrics=chart['barycentrics'])
        rr=dict(positions=result['positions'],surfaces=result['surfaces'],charts=fields)
        nn,uu=replay_corner_fields(np.vstack([v,bonev]),{'GG:outer':prior.faces,'fixed:bone':bonef+len(v)},
            {'GG:outer':attrs['normals'][prior.faces],'fixed:bone':bonen[bonef]},
            {'GG:outer':attrs['uvs'][prior.faces],'fixed:bone':boneuv[bonef]},rr)
        ordered=[];normal=[];uv=[]
        for use in result['cells']['GG']:
            if use['side']!=1:raise ValueError('Unexpected negative GG surface use requires explicit replay handling')
            sid=use['surface'];ordered.append(result['surfaces'][sid]);normal.append(nn[sid]);uv.append(uu[sid])
        faces=np.vstack(ordered);corners=result['positions'][faces];packed_faces=np.arange(corners.size//3).reshape(-1,3)
        candidate=Mesh(corners.reshape(-1,3),packed_faces,prior.frame,
            {'normals':np.vstack(normal).reshape(-1,3),'uvs':np.vstack(uv).reshape(-1,2)},
            {'record':native.source['record'],'source_catalogue_sha256':EXPECTED,
             'stage_a_geometry_sha256':prior.geometry_sha,'proposal':'Local source extension with literal bone-aware reciprocal origin trim','anatomical_acceptance':False})
        path=OUT/f'genioglossus-{side}.npz';candidate.save(path)
        tests=[]
        for bone in bones:
            check();contact,summary=pair_contacts(candidate,bone,distance_m=0.,max_pairs=25000,max_tests=200000)
            name=bone.source['record']['id'];np.savez_compressed(OUT/f'{side}-{name.replace(":","_")}-contacts.npz',**contact)
            tests.append(dict(bone_part_id=name,**summary))
        after=diagnostics(candidate,self_intersections=True)
        # Compare actual oriented corner triples, not bounding boxes, with all
        # source upper fibre faces outside this declared local origin window.
        fixed_face=np.flatnonzero(np.all(native.triangles[:,:,1]>=1.630,axis=1))
        def key(tri):
            rows=[tuple(map(tuple,np.roll(tri,i,axis=0)))for i in range(3)];return min(rows)
        emitted={key(t)for t in candidate.triangles};literal=sum(key(native.triangles[i])in emitted for i in fixed_face)
        row=dict(side=side,mesh=str(path),mesh_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            local_extension_radius_m=radius,maximum_local_extension_m=float(np.linalg.norm(v-prior.vertices,axis=1).max()),
            extended_source_vertex_ids=moved.tolist(),witness_source_muscle_face=witness['first_face_id'],witness_source_bone_face=witness['second_face_id'],
            actual_fixed_bone_geometry_sha256=[b.geometry_sha for b in bones],source_bones_modified=False,
            upper_source_triangles=len(fixed_face),upper_source_triangles_retained_literal=literal,
            boolean=result['evidence'],topology=topology(candidate),open3d=after,bone_contacts=tests,
            fixed_mandible_classifier=classifier,
            numerical_geometry_passed=bool(after['self_intersection_pairs']==0 and all(t['crossing_face_pairs']==0 for t in tests)and literal==len(fixed_face)),
            installed=False,anatomical_acceptance=False)
        rows.append(row)
    report=dict(schema='native-gg-bone-aware-attachment/v2-stage-b',elapsed_seconds=time.monotonic()-START,
        native_deadline_seconds=50,source_catalogue_sha256=EXPECTED,skeleton_modified=False,candidates=rows,installed=False,anatomical_acceptance=False)
    (OUT/'candidate.json').write_text(json.dumps(report,indent=2)+'\n')
    spec=json.loads((BASE/'native-attachments-v1/review-spec.json').read_text())
    spec['meshes']=[r for r in spec['meshes']if not any(r['mesh'].endswith(f'native-{i:02d}.npz')for i in[6,7])and'styloglossus-'not in r['mesh']]
    for row in spec['meshes']:
        for side,index in [('left',1),('right',2)]:
            if row['mesh'].endswith(f'native-{index:02d}.npz'):row['mesh']=str(OUT/f'genioglossus-{side}.npz');row['name']=f'GG bone-aware attachment candidate ({side})'
    (OUT/'review-spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
