"""Read-only actual triangle contact comparison for a saved LCFA course."""
from pathlib import Path
import json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1'
FRAME='model:local-ta98:assembled-metres'

def run(side):
    start=time.monotonic();out=B/'context-proof-v1';out.mkdir(parents=True,exist_ok=True)
    cat=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(cat)=='b26860231a10b4b4382802751d6a8c15a7491e8baa83d629d73b93f5abb61707'
    atlas=json.loads(cat.read_bytes());reader=AtlasReader();parts={p['id']:p for p in atlas['parts']}
    z=np.load(B/'bilateral-course-proof-v2'/f'{side}-course.npz')
    if side=='l':old=np.load(B/'course-correction-v3/candidate.npz')['source_vertices']
    else:old=np.load(B/'source-gate-v1/r-component-0.npz')['vertices']
    new=z['vertices'];faces=z['faces'];before=Mesh(old,faces,FRAME);after=Mesh(new,faces,FRAME)
    source=reader.read_part(atlas,parts['ZA:Lateral circumflex femoral artery.'+side],frame=FRAME)
    assert np.array_equal(source.triangles,before.triangles),'Source literal triangle order changed'
    changed=np.any(np.any(before.triangles!=after.triangles,axis=2),axis=1)
    changed_ids=np.flatnonzero(changed);assert len(changed_ids)>0
    swept=np.concatenate([before.triangles[changed],after.triangles[changed]])
    lo=swept.min((0,1));hi=swept.max((0,1))
    # Root/fork ownership is assessed over full course, other context only over
    # the changed faces. Literal unchanged geometry cannot create new contacts.
    roots={'ZA:Deep femoral artery.'+side,'TA98NATIVE:male:LCFA:ascending:'+side}
    systems={'skeletal','arterial','venous','nervous','muscular','connective'}
    peers=[p for p in atlas['parts'] if p['id']!=source.source['record']['id'] and not p.get('defaultHidden') and p.get('system') in systems and (p['id'] in roots or (np.all(np.asarray(p['bounds'])[1]>=lo) and np.all(np.asarray(p['bounds'])[0]<=hi)))]
    rows=[]
    for part in peers:
        mesh=reader.read_part(atlas,part,frame=FRAME)
        if part['id'] not in roots:
            tri=mesh.triangles;mask=np.all(tri.max(1)>=lo,axis=1)&np.all(tri.min(1)<=hi,axis=1)
            if not mask.any():continue
            peer_ids=np.flatnonzero(mask);mesh=Mesh(mesh.vertices,mesh.faces[mask],FRAME)
            a=Mesh(old,faces[changed],FRAME);b=Mesh(new,faces[changed],FRAME);source_ids=changed_ids
        else:a=before;b=after;source_ids=np.arange(len(faces));peer_ids=np.arange(len(mesh.faces))
        tree=BVH(mesh)
        first,s0=pair_contacts(a,mesh,0.,second_bvh=tree,max_tests=500000,max_pairs=100000)
        second,s1=pair_contacts(b,mesh,0.,second_bvh=tree,max_tests=500000,max_pairs=100000)
        assert s0['complete_within_pair_budget'] and s1['complete_within_pair_budget']
        def meaningful(q):
            ids=np.flatnonzero(np.isin(q['relation'],['crossing','coplanar-overlap']))
            return {(int(source_ids[q['pairs'][i,0]]),int(peer_ids[q['pairs'][i,1]]),str(q['relation'][i]))for i in ids}
        oldpairs,newpairs=meaningful(first),meaningful(second)
        added=sorted(newpairs-oldpairs);removed=sorted(oldpairs-newpairs)
        row=dict(id=part['id'],system=part['system'],unchanged_interface=part['id'] in roots,before=len(oldpairs),after=len(newpairs),added_count=len(added),removed_count=len(removed),added_examples=added[:12],removed_examples=removed[:12],actual_tests=s0['candidate_triangle_pairs_tested']+s1['candidate_triangle_pairs_tested'])
        if part['id'] in roots:
            assert oldpairs==newpairs,'Root contact pairs changed'
            used=np.unique(faces[np.unique([r[0]for r in oldpairs])]) if oldpairs else np.array([],int)
            assert np.array_equal(old[used],new[used]),'Root contacting face packet moved'
        rows.append(row)
        atomic_json(out/f'{side}-peer-{len(rows):03}.json',row)
    risks=[r for r in rows if r['added_count']]
    proof=dict(schema='native-lcfa-course-context/v1',side=side,source_catalogue_sha256=digest_file(cat),candidate_sha256=digest_file(B/'bilateral-course-proof-v2'/f'{side}-course.npz'),changed_faces=len(changed_ids),broadphase_records=len(peers),actual_peer_count=len(rows),scope='Exact triangle contacts over changed course faces; complete original deep-femoral/ascending contact packets frozen. No anatomical acceptance or lumen patency claim.',rows=rows,new_contact_risks=risks,passed=not risks,elapsed_seconds=time.monotonic()-start)
    atomic_json(out/f'{side}-proof.json',proof);print(json.dumps(dict(side=side,passed=proof['passed'],elapsed_seconds=proof['elapsed_seconds'],checked=len(rows),nonzero=[r for r in rows if r['before'] or r['after']])))

if __name__=='__main__':run(sys.argv[1])
