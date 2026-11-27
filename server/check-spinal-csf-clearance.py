"""Literal emitted arachnoid–pia clearance; no inferred CSF solid or join."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts

ap=argparse.ArgumentParser();ap.add_argument('--pia',required=True);ap.add_argument('--coats',required=True);ap.add_argument('--output',required=True);ap.add_argument('--begin',type=int,default=0);ap.add_argument('--count',type=int,default=1);args=ap.parse_args();start=time.monotonic()
def free_surface(bank,family,segment=None):
    bank=Path(bank).resolve();d=json.loads((bank/'tissue.json').read_text());assert d['passed'];raw=(bank/'cells.bin').read_bytes();selected=[]
    for name,p in d['parts'].items():
        if name.split('@')[0]!=family:continue
        if segment and name!='pia@'+segment:continue
        vertices=np.ndarray((p['vertexCount'],3),'<f4',buffer=raw,offset=p['positions']).astype(float)+[0.,1.35,-.05]
        f=np.ndarray((p['indexCount']//3,3),'<u4',buffer=raw,offset=p['indices'])
        for g in p['tissueDrawGroups']:
            if not g.get('neighbor'):selected.append(f[g['start']//3:(g['start']+g['count'])//3])
    assert selected
    return Mesh(vertices,np.concatenate(selected),d['frame']),digest_file(bank/'cells.bin')
manifest=json.loads((Path(args.pia)/'tissue.json').read_text());segments=sorted(name.split('@')[1]for name in manifest['parts']if name.startswith('pia@'))
arachnoid,asha=free_surface(args.coats,'arachnoid');y=arachnoid.triangles[:,:,1];low,high=y.min(1),y.max(1);results=[]
for segment in segments[args.begin:args.begin+args.count]:
    pia,psha=free_surface(args.pia,'pia',segment);margin=.005
    nearby=(high>=pia.bounds[0,1]-margin)&(low<=pia.bounds[1,1]+margin)
    local=Mesh(arachnoid.vertices,arachnoid.faces[nearby],arachnoid.frame);pt,at=BVH(pia),BVH(local)
    r,s=pair_contacts(pia,local,0.,first_bvh=pt,second_bvh=at,max_pairs=200000,max_tests=2000000);assert not s['truncated']
    nearest=pt.closest_mesh(at);distance=float(nearest['distance_m'])
    # Any omitted triangle is at least margin away in Y. A smaller witnessed
    # minimum therefore proves the same minimum against the whole carrier.
    assert distance<margin
    passed=s['crossing_face_pairs']==0 and s['coplanar_overlap_face_pairs']==0 and distance>0
    result=dict(segment=segment,passed=passed,pial_bank_sha256=psha,meningeal_bank_sha256=asha,crossings=s['crossing_face_pairs'],coplanar_overlaps=s['coplanar_overlap_face_pairs'],
        nearest_surface_distance_m=distance,omitted_triangles_distance_lower_bound_m=margin,global_minimum_certified_from_axial_bound=True,
        actual_float32_source_surfaces=True,spinal_csf_space_not_filled=True,skeleton_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start)
    atomic_json(Path(args.output).with_name(Path(args.output).stem+'-'+segment+'.json'),result);print(json.dumps(result),flush=True);assert passed
    results.append(result)
atomic_json(args.output,dict(passed=all(r['passed']for r in results),segments=results,whole_assembly_complete=len(results)==len(segments),elapsed_seconds=time.monotonic()-start))
