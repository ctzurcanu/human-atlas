"""Read-only native Tongue retirement and literal Mouth guide clearance proof.

This never edits, rebinds or certifies the existing alimentary route. Segment
queries use the two saved Mouth guide points, not a fabricated tube or lumen.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,sys,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.distance import triangle_pairs
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
ROUTE=ROOT/'.local-models/ta98-review/visits/female.json'
FRAME='model:local-female-ta98:assembled-metres'
PID='LOCAL:female:tongue:1089'


def main(old,candidate,output):
    started=time.monotonic()
    assert candidate.resolve().is_relative_to(WORK) and output.resolve().is_relative_to(WORK)
    assert not output.exists(),'Use a new evidence path'
    hashes={str(p):digest_file(p) for p in (old,candidate,ROUTE)}
    raw=old.read_bytes();oldraw=gzip.decompress(raw)if raw[:2]==b'\x1f\x8b'else raw
    oldsha=hashlib.sha256(oldraw).hexdigest()
    assert oldsha=='770f47f802c79cf6f3e4336bf2fa7c42cf710998069221f533a1ea3f8322e470','Bind the exact original female source'
    before=json.loads(oldraw);after=json.loads(candidate.read_bytes());route=json.loads(ROUTE.read_bytes())
    bp={p['id']:p for p in before['parts']};ap={p['id']:p for p in after['parts']}
    p,q=bp[PID],ap[PID]
    keys=set(p)|set(q);changed=sorted(k for k in keys if p.get(k)!=q.get(k))
    assert changed==['suppressed'] and not p.get('suppressed') and q.get('suppressed') is True
    assert {k:v for k,v in p.items()if k!='suppressed'}=={k:v for k,v in q.items()if k!='suppressed'}
    assert before['chunks'][p['chunk']]==after['chunks'][q['chunk']]
    marker=after['finiteTongueIntegration'];assert marker['partialMode']=='female-core-only'
    assert marker['sourcePartId']==PID and marker['legacyWholeAlias'] is True
    assert PID not in marker['partIds'] and len(marker['partIds'])==6
    reader=AtlasReader();oldbank=reader.asset(before['chunks'][p['chunk']]['url'])
    newbank=reader.asset(after['chunks'][q['chunk']]['url'])
    assert digest_file(oldbank)==digest_file(newbank)
    oldmesh=reader.read_part(before,p,frame=FRAME);retired=reader.read_part(after,q,frame=FRAME)
    assert np.array_equal(oldmesh.vertices,retired.vertices) and np.array_equal(oldmesh.faces,retired.faces)
    assert set(oldmesh.attributes)==set(retired.attributes)
    assert all(np.array_equal(a,retired.attributes[k]) for k,a in oldmesh.attributes.items())
    mouth=next(s for s in route['segments']if s['name']=='Mouth')
    points=np.asarray(mouth['points'],float);assert points.shape==(2,3) and mouth['parts']==[PID]
    start=mouth['startIndex'];assert np.array_equal(points,np.asarray(route['points'][start:start+2]))
    mucosa=ap['TA98TONGUE:female:mucosa'];mesh=reader.read_part(after,mucosa,frame=FRAME)
    bvh=BVH(mesh);ends=[bvh.closest_point(v)for v in points]
    # The distance kernel explicitly retains segment distance for degenerate
    # inputs. [a,b,b] is only a numerical query, never anatomical geometry.
    query=np.broadcast_to(np.asarray([points[0],points[1],points[1]]),(len(mesh.faces),3,3))
    distances=triangle_pairs(query,mesh.triangles);i=int(np.argmin(distances['distance_m']))
    segment_t=float(np.clip(np.dot(distances['first'][i]-points[0],points[1]-points[0])/np.dot(points[1]-points[0],points[1]-points[0]),0,1))
    saved_binding=next(b for b in route['bindings']if b['id']==PID)
    unchanged_fields={k:saved_binding.get(k)==p.get(k)for k in ['positions','indices','vertexCount','bounds']}
    unchanged_fields['chunkUrl']=saved_binding['chunkUrl']==before['chunks'][p['chunk']]['url']
    stale_bindings=[]
    for binding in route['bindings']:
        source=bp[binding['id']]
        fields=[k for k in ['positions','indices','vertexCount','bounds']if binding.get(k)!=source.get(k)]
        if binding.get('chunkUrl')!=before['chunks'][source['chunk']]['url']:fields.append('chunkUrl')
        if fields:stale_bindings.append(dict(part_id=binding['id'],changed_saved_fields=fields))
    assert hashes=={str(path):digest_file(path)for path in (old,candidate,ROUTE)},'Read-only proof changed a source'
    result=dict(schema='female-tongue-retired-route-source-literal/v1',candidate=str(candidate),candidate_sha256=hashes[str(candidate)],
        old_catalogue_path=str(old),old_catalogue_sha256=oldsha,old_file_sha256=hashes[str(old)],route_sha256=hashes[str(ROUTE)],route_source_catalogue_sha256=route['sourceCatalogueSha256'],
        source_part_id=PID,changed_record_fields=changed,old_suppressed=p.get('suppressed'),new_suppressed=q['suppressed'],
        physical_record_except_suppressed_literal=True,chunk_record_literal=True,bank_path=str(oldbank),bank_sha256=digest_file(oldbank),
        decoded_geometry_and_all_attributes_literal=True,source_geometry_sha256=oldmesh.geometry_sha,
        saved_tongue_binding_field_matches=unchanged_fields,mouth_guide_points=points.tolist(),
        replacement_mucosa_id=mucosa['id'],replacement_mucosa_geometry_sha256=mesh.geometry_sha,
        endpoint_surface_witnesses=[dict(point=points[j].tolist(),distance_m=r['distance_m'],face_id=r['face_id'],barycentric=r['barycentric'].tolist(),surface_point=r['point'].tolist())for j,r in enumerate(ends)],
        complete_segment_surface_witness=dict(distance_m=float(distances['distance_m'][i]),face_id=i,segment_parameter=segment_t,
            guide_point=distances['first'][i].tolist(),surface_point=distances['second'][i].tolist(),barycentric=distances['second_barycentric'][i].tolist(),
            every_actual_triangle_tested=len(mesh.faces),clear=float(distances['distance_m'][i])>float(distances['numerical_roundoff_m'][i])),
        existing_stale_binding_fields=stale_bindings,route_bytes_unchanged=True,route_rebound=False,route_certified=False,
        physiological_lumen_claim=False,elapsed_seconds=time.monotonic()-started)
    atomic_json(output,result)
    print(json.dumps({k:result[k]for k in ['candidate_sha256','changed_record_fields','physical_record_except_suppressed_literal','decoded_geometry_and_all_attributes_literal','saved_tongue_binding_field_matches','complete_segment_surface_witness','existing_stale_binding_fields','route_bytes_unchanged','route_certified','elapsed_seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--before',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.before.resolve(),a.candidate.resolve(),a.output.resolve())
