"""Bounded retained-territory audit of a saved, original tongue candidate.

An overlap with an external-muscle envelope is recorded as unresolved ownership;
neither outer-tongue containment nor anatomical interdigitation excuses it.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts


def main(candidate):
    began=time.monotonic()
    construction=json.loads((candidate/'construction.json').read_text())
    checked=json.loads((candidate/'checked-geometry.json').read_text())
    assert checked['passed'] and checked['construction_sha256']==digest_file(candidate/'construction.json')
    source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    assert digest_file(source)==construction['source_catalogue_sha256']
    work=ROOT/'.local-models/geometry-workbench/tongue-completion'
    source_audit=work/'source-v1/source-audit.json'
    audit=json.loads(source_audit.read_text())
    context=[]
    for p in audit['contextParts'][1:5]:
        path=Path(p['meshSavePath'])
        assert digest_file(path)==p['meshSaveSha256']
        context.append(dict(key=p['partId'],name=p['name'],mesh=str(path),sha256=digest_file(path)))
    stylo=work/'styloglossus-v2/candidate.json'
    stylo_data=json.loads(stylo.read_text())
    assert stylo_data['source_catalogue_sha256']==digest_file(source)
    for side,p in stylo_data['sides'].items():
        assert digest_file(p['mesh'])==p['sha256']
        context.append(dict(key='styloglossus-v2-'+side,name='Original styloglossus candidate '+side,
                            mesh=p['mesh'],sha256=p['sha256']))
    out=candidate/'retained-context-audit'
    out.mkdir(parents=True,exist_ok=True)
    context_meshes=[]
    for row in context:
        mesh=Mesh.load(row['mesh']); row['geometry_sha256']=mesh.geometry_sha
        row['topology']=topology(mesh)
        context_meshes.append(mesh)
    rows=[]
    for intrinsic_row in construction['meshes']:
        assert digest_file(intrinsic_row['mesh'])==intrinsic_row['sha256']
        intrinsic=Mesh.load(intrinsic_row['mesh'])
        q=np.concatenate([intrinsic.vertices,intrinsic.triangles.mean(1)])
        for context_row,retained in zip(context,context_meshes):
            stage=time.monotonic()
            arrays,summary=pair_contacts(intrinsic,retained,0.,max_pairs=50000,max_tests=1000000)
            bad=np.isin(arrays['relation'],['crossing','coplanar-overlap'])
            pairs=arrays['pairs'][bad]
            inside_count=None; max_winding=None
            if context_row['topology']['closed_oriented_manifold']:
                winding=np.abs(igl.winding_number(retained.vertices,retained.faces,q))
                inside_count=int(np.count_nonzero(winding>.99)); max_winding=float(winding.max())
            reverse=np.abs(igl.winding_number(intrinsic.vertices,intrinsic.faces,retained.vertices))
            witnesses=[]
            for a,b in pairs[:20]:
                witnesses.append(dict(intrinsic_face=int(a),retained_face=int(b),
                    intrinsic_triangle=intrinsic.triangles[a].tolist(),retained_triangle=retained.triangles[b].tolist()))
            row=dict(intrinsic=intrinsic_row['key'],retained=context_row['key'],
                intrinsic_geometry_sha256=intrinsic.geometry_sha,retained_geometry_sha256=retained.geometry_sha,
                contact_summary=summary,defect_face_pairs=int(len(pairs)),exact_witnesses=witnesses,
                intrinsic_sample_queries=len(q),intrinsic_samples_inside_retained_envelope=inside_count,
                maximum_absolute_retained_winding=max_winding,
                retained_vertices_inside_intrinsic=int(np.count_nonzero(reverse>.99)),
                interpretation=('Unresolved tissue ownership: literal envelopes cross or contain each other; anatomical interdigitation is not proved'
                    if len(pairs) or inside_count or np.any(reverse>.99)
                    else 'No crossing or sampled containment detected; this does not establish an attachment or complete anatomy'),
                elapsed_seconds=time.monotonic()-stage)
            row_file=out/(intrinsic_row['key']+'--'+context_row['key'].replace(':','_').replace('/','_')+'.json')
            atomic_json(row_file,row)
            np.savez_compressed(row_file.with_suffix('.npz'),**arrays)
            rows.append(row)
    unresolved=[r for r in rows if r['defect_face_pairs'] or r['intrinsic_samples_inside_retained_envelope'] or r['retained_vertices_inside_intrinsic']]
    complete=all(r['contact_summary']['complete_within_pair_budget'] for r in rows)
    report=dict(schema='tongue-intrinsic-retained-context-audit/v1',
        source_catalogue_sha256=digest_file(source),source_audit_sha256=digest_file(source_audit),
        construction_sha256=digest_file(candidate/'construction.json'),
        candidate_check_sha256=digest_file(candidate/'checked-geometry.json'),
        styloglossus_candidate_sha256=digest_file(stylo),context=context,pairs=rows,
        pair_count=len(rows),expected_pair_count=7*6,complete_within_pair_budget=complete,
        crossing_face_pairs=sum(r['contact_summary']['crossing_face_pairs'] for r in rows),
        coplanar_overlap_face_pairs=sum(r['contact_summary']['coplanar_overlap_face_pairs'] for r in rows),
        unresolved_ownership_pairs=len(unresolved),anatomical_acceptance=False,
        assembly_exclusive_tissue_verified=False,source_geometry_modified=False,installed=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'summary.json',report)
    print(json.dumps({k:report[k] for k in ['pair_count','complete_within_pair_budget','crossing_face_pairs',
        'coplanar_overlap_face_pairs','unresolved_ownership_pairs','elapsed_seconds']}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True,type=Path)
    args=ap.parse_args();candidate=args.candidate.resolve()
    if not candidate.is_relative_to(ROOT/'.local-models/geometry-workbench/tongue-completion'):
        raise ValueError('Audit only saved local tongue candidates')
    main(candidate)
