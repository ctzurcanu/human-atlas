"""Bounded source-network, section gauge and fixed bone/skin proposal checks.

Actual triangle crossing counts and section area ranges are diagnostics, not
proof of anatomical clearance, finite vascular walls or perfusion.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.audit import topology
from geometry_tools.contact import self_contact_summary
from geometry_tools.course_io import sha, write_json
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.commands.filter_hand_components import crossings
from geometry_tools.algorithms.source_tube_topology import recover_regular_source_annuli


def combine(meshes):
    vertices, faces, offset = [], [], 0
    for m in meshes:
        vertices.append(m.vertices)
        faces.append(m.faces+offset)
        offset += len(m.vertices)
    return Mesh(np.concatenate(vertices), np.concatenate(faces), meshes[0].frame)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['candidate', 'bones', 'skin-reference', 'output']:
        p.add_argument('--'+key, required=True)
    p.add_argument('--parts',help='Scoped diagnostic only; comma-separated source part IDs, never whole-assembly acceptance')
    p.add_argument('--profiles',help='Explicit current source profile bank; no row-order pairing fallback is allowed')
    args = p.parse_args()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve()
             for k in ['candidate', 'bones', 'skin-reference']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    banks={}
    if args.profiles:
        paths['profiles']=Path(args.profiles).resolve()
        banks={r['partId']:r for r in json.loads(paths['profiles'].read_text())['rows']}
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new diagnostic evidence file')
    started = time.monotonic()
    scope=set(args.parts.split(','))if args.parts else None
    if scope and not scope.issubset({r['partId']for r in data['candidate']['rows']}):raise ValueError('Unknown source part in diagnostic scope')
    skin = data['skin-reference']
    if not skin['defaultHiddenExcluded'] or sha(skin['path']) != skin['sha256']:
        raise ValueError('Active skin source changed')
    skin_tree = BVH(Mesh.load(skin['path']))
    rows, whole, scoped = [], [], []
    for side in ['left', 'right']:
        selected = [r for r in data['candidate']['rows'] if r['side'] == side]
        if scope:selected=[r for r in selected if r['partId']in scope]
        if not selected:continue
        m = combine([Mesh.load(r['path']) for r in selected])
        summary, _ = self_contact_summary(m)
        topo = topology(m)
        assembly=dict(side=side, selfPairs=summary['defect_pairs'],
                          selfComplete=summary['complete_within_pair_budget'],
                          components=len(topo['components']), nonmanifoldEdges=topo['nonmanifold_edges'],
                          invalidVertexLinks=len(topo['invalid_vertex_links']))
        (scoped if scope else whole).append(assembly)
        bone_rows = [r for r in data['bones']['rows']
                     if r['partId'].endswith('.l' if side == 'left' else '.r')]
        for row in bone_rows:
            if sha(row['path']) != row['sha256']:
                raise ValueError('Fixed bone reference changed')
        bone_tree = BVH(combine([Mesh.load(r['path']) for r in bone_rows]))
        for row in selected:
            if sha(row['path']) != row['sha256'] or sha(row['sourcePath']) != row['sourceSha256']:
                raise ValueError('Source/candidate array changed')
            a, b = Mesh.load(row['sourcePath']), Mesh.load(row['path'])
            if not np.array_equal(a.faces, b.faces):
                raise ValueError('This checker requires literal source face ancestry')
            pairing = load_source_pairing(banks.get(row['partId'],row.get('sourceTrianglePairing')), a.faces)
            rec = recover_regular_source_annuli(a.vertices, a.faces, pairing)
            ratios = []
            for ring in rec['rings']:
                ids = rec['representative_raw_ids'][ring]
                def area(points):
                    q = points-points.mean(0)
                    return np.linalg.norm(np.cross(q, np.roll(q, -1, axis=0)).sum(0))/2
                ratios.append(area(b.vertices[ids])/area(a.vertices[ids]))
            before = [crossings(a, tree, started) for tree in [bone_tree, skin_tree]]
            after = [crossings(b, tree, started) for tree in [bone_tree, skin_tree]]
            record=dict(partId=row['partId'], side=side, before=before, after=after,
                             profiles=len(ratios), areaRange=[min(ratios, default=1), max(ratios, default=1)])
            if scope:
                own,_=self_contact_summary(b);record.update(ownSelfPairs=own['defect_pairs'],ownSelfComplete=own['complete_within_pair_budget'])
                adjacent=set()
                for pair in data['candidate']['pairs']:
                    if pair['firstPart']==row['partId']:adjacent.add(pair['secondPart'])
                    if pair['secondPart']==row['partId']:adjacent.add(pair['firstPart'])
                other=[];box=np.array([b.vertices.min(0),b.vertices.max(0)])
                for peer in data['candidate']['rows']:
                    if peer['side']!=side or peer['partId']==row['partId']or peer['partId']in adjacent:continue
                    if sha(peer['path'])!=peer['sha256']:raise ValueError('Unjoined source artery changed')
                    pm=Mesh.load(peer['path'])
                    if np.any(pm.vertices.max(0)<box[0])or np.any(pm.vertices.min(0)>box[1]):continue
                    count=crossings(b,BVH(pm),started)
                    if count:other.append(dict(partId=peer['partId'],pairs=count))
                record.update(unjoinedArteryCrossings=other,adjoiningPortsExcludedPendingRecoupling=sorted(adjacent))
            rows.append(record)
        if time.monotonic()-started > 22:
            raise TimeoutError('22-second complete network diagnostic')
    totals = {side: np.sum([[r['before'], r['after']] for r in rows if r['side'] == side], axis=0).tolist()
              for side in ['left', 'right']}
    report = dict(candidateSha256=sha(paths['candidate']), inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  rows=rows, wholeAssembly=whole, boneSkinTotalsBeforeAfter=totals,
                  sectionAreaRange=[min(r['areaRange'][0] for r in rows), max(r['areaRange'][1] for r in rows)],
                  complete=True,fullAssemblyChecked=scope is None,scopedSourceParts=None if scope is None else sorted(scope),scopedSurfaces=scoped,
                  anatomicalAcceptance=False, sourceCodeSha256=sha(__file__),
                  elapsedSeconds=time.monotonic()-started)
    write_json(out, report)
    print(json.dumps({k: report[k] for k in ['wholeAssembly', 'boneSkinTotalsBeforeAfter', 'sectionAreaRange', 'elapsedSeconds']}))


if __name__ == '__main__':
    main()
