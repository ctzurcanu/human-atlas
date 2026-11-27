"""Literal four-cell conduit replay; remove only three reciprocal tissue annuli."""
import json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.algorithms.tube_continuity import continuity_audit,boundary_loops
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/renal-blood-frontier-v1';CAND=BASE/'new-A-arterial-walls-v3'

def key(t):return b''.join(sorted(x.astype('<f8').tobytes()for x in t))

def area(p):
    q=p-p.mean(0);return float(np.linalg.norm(np.cross(q,np.roll(q,-1,axis=0)).sum(0))/2)

def main():
    began=time.monotonic();out=CAND/'independent-conduit-data-v1';assert not out.exists();out.mkdir();ev=json.loads((BASE/'evidence.json').read_text());cat=Path(ev['catalogue']);assert digest_file(cat)==ev['catalogue_sha256'];reader=AtlasReader(ROOT)
    a=json.loads(cat.read_text());source_record=next(p for p in a['parts']if p['id']=='ZA:Anterior branch of renal artery.l');native=reader.read_part(a,source_record,frame='model:local-ta98:assembled-metres',catalogue_path=cat,catalogue_sha=ev['catalogue_sha256'])
    snapshot=Mesh.load(BASE/'source-01.npz');assert(native.geometry_sha,native.attributes_sha)==(snapshot.geometry_sha,snapshot.attributes_sha)
    tm=json.loads((CAND/'finite-cells/tissue.json').read_text());bank=Path(tm['bank']);assert digest_file(bank)==tm['bank_sha256'];synthetic={'chunks':[{'url':'/local-models/'+str(bank.relative_to(ROOT/'.local-models'))}]};meshes={'native':native};records={'native':source_record}
    for name,record in tm['parts'].items():
        m=reader.read_part(synthetic,dict(record,chunk=0),frame=tm['frame']);assert(m.geometry_sha,m.attributes_sha)==(tm['checks']['cells'][name]['geometry_sha256'],tm['checks']['cells'][name]['attributes_sha256']);meshes[name]=m;records[name]=record
    def group(name,index):
        g=records[name]['tissueDrawGroups'][index];return np.arange(g['start']//3,(g['start']+g['count'])//3)
    interfaces=[('native',2,'interlobar',2),('interlobar',3,'arcuate',2),('arcuate',3,'cortical-radiate',2)];removed={name:set()for name in meshes};proofs=[]
    for index,(an,ag,bn,bg)in enumerate(interfaces):
        ma,mb=meshes[an],meshes[bn];ai,bi=group(an,ag),group(bn,bg);assert len(ai)==len(bi)==24;ta,tb=ma.triangles[ai],mb.triangles[bi];lookup={key(t):i for i,t in enumerate(tb)};assert len(lookup)==24
        mapping=[];raw_equal=True;effective_opposed=True;winding_opposed=True;signs_opposed=True
        for i,t in enumerate(ta):
            j=lookup[key(t)];order=np.array([np.flatnonzero(np.all(tb[j]==p,axis=1))[0]for p in t]);mapping.append([i,j,*order.tolist()]);assert np.array_equal(t,tb[j][order])
            n1=np.cross(t[1]-t[0],t[2]-t[0]);n2=np.cross(tb[j,1]-tb[j,0],tb[j,2]-tb[j,0]);winding_opposed&=bool(n1@n2<0)
            va=ma.faces[ai[i]];vb=mb.faces[bi[j]][order];wa=ma.attributes['normals'][va].astype(np.int64);wb=mb.attributes['normals'][vb].astype(np.int64);sa=ma.attributes.get('normalSigns',np.ones(len(ma.vertices),np.int8))[va];sb=mb.attributes.get('normalSigns',np.ones(len(mb.vertices),np.int8))[vb]
            raw_equal&=bool(np.array_equal(wa,wb));signs_opposed&=bool(np.array_equal(sa,-sb));effective_opposed&=bool(np.array_equal(wa*sa[:,None],-wb*sb[:,None]))
        assert winding_opposed and effective_opposed,'Preserve actual mismatched-interface failure; never average normals or move controls'
        for n,ids in [(an,ai),(bn,bi)]:removed[n].update(map(int,ids))
        p=out/f'interface-{index}.npz';atomic_npz(p,first_actual_face_ids=ai,second_actual_face_ids=bi,corner_correspondence=np.array(mapping),first_triangles=ta,second_triangles=tb)
        proofs.append(dict(first=an,second=bn,first_faces=24,second_faces=24,literal_world_coordinates=True,opposite_triangle_winding=winding_opposed,raw_packed_normal_words_equal=raw_equal,packed_normal_signs_opposed=signs_opposed,effective_packed_normal_words_literal_opposed=effective_opposed,witnesses=str(p),sha256=digest_file(p)))
    vertices=[];allfaces=[];innerfaces=[];outerfaces=[];ancestry=[];offset=0;source_summaries=[]
    for name in ['native','interlobar','arcuate','cortical-radiate']:
        m=meshes[name];kept=np.setdiff1d(np.arange(len(m.faces)),np.array(sorted(removed[name]),int));vertices.append(m.vertices);allfaces.append(m.faces[kept]+offset);innerfaces.append(m.faces[group(name,1)]+offset);outerfaces.append(m.faces[group(name,0)]+offset);ancestry.extend([(name,int(i))for i in kept]);offset+=len(m.vertices)
        source_summaries.append(dict(name=name,geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,original_faces=len(m.faces),removed_reciprocal_faces=len(removed[name]),kept_literal_faces=len(kept)))
    vertices=np.concatenate(vertices);faces=np.concatenate(allfaces);joined=Mesh(vertices,faces,native.frame,source=dict(kind='diagnostic-literal-source-face-union',source_parts=source_summaries,no_new_or_moved_controls=True,no_new_faces=True));check=audit(joined,require_closed=True,max_pairs=500000,max_tests=4000000)
    expected_genus=[c['closed_orientable_genus']for c in check['topology']['components']];assert check['passed']and expected_genus==[1]
    inner=continuity_audit(vertices,np.concatenate(innerfaces),expected_mouths=2,expected_components=1,expected_euler=0);outer=continuity_audit(vertices,np.concatenate(outerfaces),expected_mouths=2,expected_components=1,expected_euler=0);assert inner['topological_continuity_pass']and outer['topological_continuity_pass']
    terminal=[]
    for name,gi in [('native',3),('cortical-radiate',3)]:
        m=meshes[name];f=m.faces[group(name,gi)];v,iv=np.unique(m.vertices,axis=0,return_inverse=True);pf=iv[f];loops=boundary_loops(pf);assert len(loops)==2 and all(len(x)==12 for x in loops);loops.sort(key=lambda r:area(v[r]),reverse=True);o,i=loops;assert set(np.unique(pf))==set(o)|set(i)
        mixed=all(any(x in set(o)for x in t)and any(x in set(i)for x in t)for t in pf);assert mixed
        terminal.append(dict(owner=name,annular_faces=24,outer_controls=12,inner_controls=12,outer_area_m2=area(v[o]),inner_area_m2=area(v[i]),every_triangle_spans_wall_thickness_only=mixed,central_lumen_disk=False))
    paired=np.load(CAND/'paired-profiles.npz');fields=paired['fields'];origin=np.array(tm['parts']['interlobar']['geometryOrigin']);encoded=(fields-origin).astype('<f4').astype(float)+origin;pa=np.array([[area(p)for p in layer]for layer in encoded]);assert np.all(pa[1]>0)and np.all(pa[0]>pa[1]);stops=paired['region_stops'];profile_replays=[]
    for name,start,stop in [('interlobar',stops[0],stops[1]),('arcuate',stops[1],stops[2]),('cortical-radiate',stops[2],stops[3])]:
        m=meshes[name]
        for role,gindex,field in [('outer',0,0),('inner',1,1)]:
            points=m.vertices[np.unique(m.faces[group(name,gindex)])];actual=set(map(tuple,points));planned=set(map(tuple,encoded[field,int(start):int(stop)+1].reshape(-1,3)));assert actual==planned
            profile_replays.append(dict(cell=name,surface=role,all_encoded_profile_controls_literal=True,stations=int(stop-start+1),controls_per_station=12))
    diag=out/'joined-literal-wall.npz';joined.save(diag);atomic_npz(out/'joined-face-ancestry.npz',source_owner=np.array([a for a,b in ancestry]),source_face=np.array([b for a,b in ancestry]),inner_faces=np.concatenate(innerfaces),outer_faces=np.concatenate(outerfaces))
    result=dict(schema='independent-renal-A-conduit-proof/v1',passed=True,source_catalogue_sha256=ev['catalogue_sha256'],native_part='ZA:Anterior branch of renal artery.l',native_geometry_sha256=native.geometry_sha,native_attributes_sha256=native.attributes_sha,candidate_bank_sha256=tm['bank_sha256'],candidate_manifest_sha256=digest_file(CAND/'finite-cells/tissue.json'),
        reciprocal_interfaces=proofs,only_removed_faces_are_three_reciprocal24face_annuli=True,total_removed_tissue_faces=144,sources=source_summaries,joined_wall_check=check,expected_one_connected_genus_one_wall=expected_genus==[1],inner_sheet=inner,outer_sheet=outer,terminal_wall_annuli=terminal,
        paired_profile_replays=profile_replays,paired_profile_sha256=digest_file(CAND/'paired-profiles.npz'),minimum_inner_profile_area_m2=float(pa[1].min()),every_inner_profile_positive_area=True,no_central_lumen_disk=True,unrelated_source_faces_and_attributes_literal=True,
        joined_diagnostic=str(diag),joined_diagnostic_sha256=digest_file(diag),whole_physiological_patency_not_claimed=True,whole_urine_occupancy_not_claimed=True,geometry_or_master_modified=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    assert digest_file(cat)==ev['catalogue_sha256'];atomic_json(CAND/'independent-conduit-v1.json',result);print(json.dumps(dict(output=str(CAND/'independent-conduit-v1.json'),passed=True,seconds=result['elapsed_seconds'],faces=len(faces),inner_mouths=inner['mouths'],minimum_inner_area_m2=result['minimum_inner_profile_area_m2'],interfaces=proofs)))
if __name__=='__main__':main()
