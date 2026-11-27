"""Partition the existing renal collecting wall by its two actual carriers."""
import numpy as np
import point_cloud_utils as pcu
from .renal_master import bound,OUT,FRAME
from .common import atomic_json,atomic_npz
from .harmonic_wall import layers
from .model import Mesh
from .audit import audit,topology


def build(version='446',*,method='distance-ratio'):
    out=OUT/('collecting-harmonic-v'+version)
    if out.exists():raise FileExistsError(out)
    _,z=bound('renal-ureter-fluid-v317');v=z['positions'];f=z['wholeWallFaces'];tri=v[f]
    def distance(points,faces):return abs(pcu.signed_distance_to_mesh(np.ascontiguousarray(points),v,faces.astype('i4'))[0])
    do=distance(tri.reshape(-1,3),z['outerClassifier']).reshape(-1,3);di=distance(tri.reshape(-1,3),z['innerClassifier']).reshape(-1,3)
    outer=(do.max(1)<2e-12)&(distance(tri.mean(1),z['outerClassifier'])<2e-12)
    inner=(di.max(1)<2e-12)&(distance(tri.mean(1),z['innerClassifier'])<2e-12)
    q=v[z['innerClassifier']];common=distance(q.reshape(-1,3),z['outerClassifier']).reshape(-1,3).max(1)<2e-12
    caps=q[common];centers=caps.mean(1);normal=np.cross(caps[:,1]-caps[:,0],caps[:,2]-caps[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    end=np.zeros(len(f),bool)
    for c,n in zip(centers,normal):
        possible=np.linalg.norm(tri.mean(1)-c,axis=1)<.006
        end|=possible&(np.max(abs((tri-c)@n),axis=1)<3e-9)
    outer&=~end;inner&=~end
    ambiguous=outer&inner;outer&=~ambiguous;inner&=~ambiguous
    shared=np.intersect1d(np.unique(f[outer]),np.unique(f[inner]))
    print('Carriers',int(outer.sum()),int(inner.sum()),'rim faces',int(end.sum()),'shared controls',len(np.intersect1d(np.unique(f[outer]),np.unique(f[inner]))),flush=True)
    original_f=f;outer_carrier=f[outer];inner_carrier=f[inner]
    # Resolve the annular attachment surface itself; a whole facet whose
    # corners lie on the common outer/inner rim cannot carry a 0/0 ratio.
    # Adding its centroid preserves its exact planar geometric boundary.
    vertices=v.tolist();refined=[];parents=[]
    for i,face in enumerate(f):
        if method.startswith('distance-')and not outer[i]and not inner[i]:
            center=len(vertices);vertices.append(v[face].mean(0).tolist())
            refined.extend([[int(a),int(b),center]for a,b in zip(face,np.roll(face,-1))]);parents.extend([i]*3)
        else:refined.append(face.tolist());parents.append(i)
    v=np.array(vertices);f=np.array(refined,int);parents=np.array(parents,int)
    r=layers(v,f,outer_carrier,inner_carrier,['adventitia','muscularis','mucosa'],origin=[.0625,1.1875,-.015625],free_source_vertices=shared,method=method)
    arrays={'positions':r['positions']};rows=[]
    for i,(sid,fs)in enumerate(r['surfaces'].items()):
        tags=r['source_faces'][sid].copy();valid=tags>=0;tags[valid]=parents[tags[valid]]
        arrays['surface_'+str(i)]=fs;arrays['source_face_'+str(i)]=tags;rows.append(dict(id=sid,array='surface_'+str(i),source_face_array='source_face_'+str(i)))
    atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'result.json',dict(surfaces=rows,cells=r['cells'],evidence=r['evidence'],carrier_outer_faces=int(outer.sum()),carrier_inner_faces=int(inner.sum()),rim_faces=int(end.sum())))
    for name,uses in r['cells'].items():
        fs=np.concatenate([r['surfaces'][u['surface']]if u['side']==1 else r['surfaces'][u['surface']][:,[0,2,1]]for u in uses]);m=Mesh(r['positions'],fs,FRAME);m.save(out/(name+'.npz'))
        t=topology(m);topology_ok=t['closed_oriented_manifold']and [c['closed_orientable_genus']for c in t['components']]==[11]
        q=audit(m,self_intersections=topology_ok,require_closed=True,max_pairs=1000000,expected_genus=[11]);atomic_json(out/(name+'-audit.json'),q)
        print(name,len(fs),q['passed'],q['failures'],[x['closed_orientable_genus']for x in q['topology']['components']],flush=True)
    return r


if __name__=='__main__':build()
