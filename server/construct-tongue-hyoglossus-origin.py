"""One original broad lower-HG sheet proposal between literal source disks.

The native superior muscle and actual fixed hyoid footprint remain literal.
New intervening wall controls form a shallow curved muscular sheet. This is
one saved review candidate, not an installation or an acceptance assertion.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.tube_rim_bridge import bridge_boundaries

BASE=ROOT/'.local-models/geometry-workbench/tongue-completion'
SOURCE=BASE/'source-v1/meshes'
GUIDES=BASE/'native-attachment-candidate-v1'
OUT=BASE/'hyoglossus-v1'
START=time.monotonic()
DEADLINE=START+50


def check():
    if time.monotonic()>DEADLINE:raise TimeoutError('HG construction exceeded 50 second native deadline')


def load(i):return Mesh.load(SOURCE/f'native-{i:02d}.npz')


def normal_words(v,f):
    n=np.zeros_like(v);t=v[f];cross=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    for j in range(3):np.add.at(n,f[:,j],cross)
    length=np.linalg.norm(n,axis=1);n[length>0]/=length[length>0,None]
    return np.rint(np.clip(n,-1,1)*32767).astype('<i2')


def subdivide(v,band,start_t):
    # One deterministic subdivision introduces a shared midpoint surface;
    # finite edges remain shared and original boundary coordinates unchanged.
    points=v.tolist();phase=list(start_t);cache={};out=[]
    def middle(a,b):
        key=tuple(sorted((int(a),int(b))))
        if key not in cache:
            cache[key]=len(points);points.append(((np.array(points[a])+points[b])/2).tolist());phase.append((phase[a]+phase[b])/2)
        return cache[key]
    for a,b,c in band:
        ab,bc,ca=middle(a,b),middle(b,c),middle(c,a)
        out.extend([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]])
    return np.array(points),np.array(out,int),np.array(phase)


def compact(mesh):
    ids,inv=np.unique(mesh.faces,return_inverse=True)
    return Mesh(mesh.vertices[ids],inv.reshape(-1,3),mesh.frame,
        {k:v[ids]for k,v in mesh.attributes.items()},mesh.source)


def main():
    OUT.mkdir(parents=True,exist_ok=False)
    proposals=json.loads((GUIDES/'hyoglossus-origin-proposal.json').read_text())['proposals']
    bones=[load(i)for i in [5,13,14,15,16,17]]
    all_rows=[]
    source_diagnostics={}
    for i in [1,2]:source_diagnostics[f'native-{i:02d}']=diagnostics(load(i),self_intersections=True)
    (OUT/'native-gg-self-audit.json').write_text(json.dumps(source_diagnostics,indent=2)+'\n')
    for guide,mi,sgn in zip(proposals,[3,4],[1,-1]):
        check();side=guide['side'];native=load(mi)
        upper=compact(Mesh.load(GUIDES/f'hyoglossus-{side}-preserved-upper.npz'))
        origin=Mesh.load(guide['origin_surface']['mesh'])
        v=np.vstack([upper.vertices,origin.vertices]);offset=len(upper.vertices)
        cap=origin.faces[:,[0,2,1]]+offset
        upperloop=boundary_loops(upper.faces);originloop=boundary_loops(cap)
        if len(upperloop)!=1 or len(originloop)!=1:raise ValueError('Explicit preservation and origin disks must each have one boundary')
        band,loft=bridge_boundaries(v,upperloop[0],originloop[0],align_projected_phase=False)
        parameter=np.r_[np.ones(len(upper.vertices)),np.zeros(len(origin.vertices))]
        p,wall,phase=subdivide(v,band,parameter)
        new_ids=np.arange(len(v),len(p));t=phase[new_ids]
        # A common shallow lateral bow changes the course of both wall sides,
        # preserving a broad sheet, rather than creating a cylindrical cord.
        p[new_ids,0]+=sgn*.0007*np.sin(np.pi*t)
        f=np.vstack([upper.faces,cap,wall]);words=normal_words(p,f)
        words[:len(upper.vertices)]=upper.attributes['normals']
        # Native boundary edges must stay unsplit in the assembled cell. The
        # refinement boundary midpoints are not emitted on the cap constraints;
        # reconcile them by leaving these source edge triangles unsplit.
        # The deterministic step below triangulates the original coarse rim
        # against its deformed cross-edge midpoints, preserving source chains.
        boundary_edges={tuple(sorted((int(a),int(b))))for ring in [upperloop[0],originloop[0]]for a,b in zip(ring,np.roll(ring,-1))}
        # Rebuild each coarse rim triangle with its non-source edge midpoints.
        # This avoids a T-junction against either retained literal source disk.
        wall_rows=[];mid_by_point={tuple(point):i for i,point in enumerate(p[:len(v)])}
        edge_mid={}
        for a,b in {(tuple(sorted((int(a),int(b)))))for tri in band for a,b in zip(tri,np.roll(tri,-1))}:
            if (a,b) in boundary_edges:continue
            raw=(v[a]+v[b])/2;pt=(parameter[a]+parameter[b])/2
            raw[0]+=sgn*.0007*np.sin(np.pi*pt)
            k=int(np.argmin(np.linalg.norm(p-raw,axis=1)));edge_mid[(a,b)]=k
        for tri in band:
            poly=[]
            for a,b in zip(tri,np.roll(tri,-1)):
                poly.append(int(a));key=tuple(sorted((int(a),int(b))))
                if key in edge_mid:poly.append(edge_mid[key])
            # Each coarse triangle owns a curved interior centroid. Its fan
            # shares the exact subdivided cross edges and literal cap edges.
            center=v[tri].mean(0);pt=parameter[tri].mean();center[0]+=sgn*.0007*np.sin(np.pi*pt)
            center_id=len(p);p=np.vstack([p,center])
            wall_rows.extend([[center_id,a,b]for a,b in zip(poly,np.roll(poly,-1))])
        wall=np.asarray(wall_rows,int);f=np.vstack([upper.faces,cap,wall]);words=normal_words(p,f)
        corner_words=words[f].copy()
        corner_words[:len(upper.faces)]=upper.attributes['normals'][upper.faces]
        owners=guide['origin_surface']['actual_source_face_owners']
        sources={load(i).source['record']['id']:load(i)for i in [13,14,15,16,17]}
        contact_words=[]
        for owner in owners:
            host=sources[owner['part']];face=host.faces[owner['face']]
            contact_words.append(np.clip(-host.attributes['normals'][face[[0,2,1]]].astype(np.int32),-32767,32767).astype('<i2'))
        corner_words[len(upper.faces):len(upper.faces)+len(cap)]=np.asarray(contact_words)
        # Corner bank records native shading charts without corrupting exact
        # physical welding or averaging different bone-part normal words.
        flat=f.reshape(-1);candidate=Mesh(p[flat],np.arange(len(flat)).reshape(-1,3),native.frame,
            {'normals':corner_words.reshape(-1,3)},source={'record':native.source['record'],
             'parent_geometry_sha256':native.geometry_sha,'proposal':'Curved broad lower-origin sheet from actual hyoid body/greater horn to literal superior HG',
             'anatomical_acceptance':False})
        path=OUT/f'hyoglossus-{side}.npz';candidate.save(path)
        tests=[]
        for bone in bones:
            check();contacts,summary=pair_contacts(candidate,bone,distance_m=0.,max_pairs=25000,max_tests=200000)
            tests.append(dict(bone_part=bone.source['record']['id'],**summary))
            np.savez_compressed(OUT/f'{side}-{bone.source["record"]["id"].replace(":","_")}.npz',**contacts)
        row=dict(side=side,mesh=str(path),mesh_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            native_part_id=native.source['record']['id'],topology=topology(candidate),open3d=diagnostics(candidate,self_intersections=True),
            fixed_bone_checks=tests,loft=loft,origin_source_face_owners=owners,
            actual_origin_triangles=len(cap),literal_retained_superior_triangles=len(upper.faces),
            native_retained_triangles_sha256=hashlib.sha256(upper.triangles.astype('<f8').tobytes()).hexdigest(),
            native_retained_corner_normal_words_sha256=hashlib.sha256(upper.attributes['normals'][upper.faces].tobytes()).hexdigest(),
            reciprocal_native_bone_contact_words_saved=True,curved_new_wall_lateral_bow_m=.0007,
            sheet_width_basis='Actual 29 mm body/greater-horn superior footprint and actual native superior insertion boundary; no microanatomical precision claimed.',
            installed=False,anatomical_acceptance=False)
        all_rows.append(row)
    report=dict(schema='hyoglossus-first-curved-sheet-candidate/v1',elapsed_seconds=time.monotonic()-START,
        native_deadline_seconds=50,skeleton_modified=False,candidates=all_rows,installed=False,anatomical_acceptance=False)
    (OUT/'candidate.json').write_text(json.dumps(report,indent=2)+'\n')
    spec=json.loads((GUIDES/'review-spec.json').read_text());spec['meshes']=[r for r in spec['meshes']if 'bone-origin-surface' not in r['mesh']]
    for row in spec['meshes']:
        if 'genioglossus-' in row['mesh']:
            side='left'if'left'in row['mesh']else'right';row['mesh']=str(SOURCE/f'native-{1 if side=="left"else 2:02d}.npz');row['name']=f'Native GG ({side})'
        if row['mesh'].endswith('native-03.npz')or row['mesh'].endswith('native-04.npz'):
            side='left'if row['mesh'].endswith('native-03.npz')else'right';row['mesh']=str(OUT/f'hyoglossus-{side}.npz');row['name']=f'First broad HG origin candidate ({side})'
    (OUT/'review-spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
