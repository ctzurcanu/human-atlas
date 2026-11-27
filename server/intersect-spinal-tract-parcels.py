"""Replace eligible white tissue with disjoint owned tract/remainder cells.

This stage does not add another superimposed white-matter model. The existing
grey, pial, outer-cord and canal boundaries are replayed in one Boolean graph.
Native funiculus membership remains attached to every resulting tissue piece.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import manifold3d as md
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.boolean_cells import BooleanCells
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['intersect','assemble']);ap.add_argument('--parcels',required=True)
ap.add_argument('--source');ap.add_argument('--output',required=True);args=ap.parse_args();started=time.monotonic()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);parcels=Path(args.parcels).resolve()
pm=json.loads((parcels/'masks.json').read_text());segment=pm['segment']
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
native=base/'female-31-segments-v2';meta=json.loads((native/'cells.json').read_text())
if args.stage=='intersect':
    z=np.load(native/'cells.npz');v=z['positions'];bank=BooleanCells(tolerance_m=2e-12);original={};solids={};attributes={}
    for name,uses in meta['cells'].items():
        if not name.endswith('@'+segment):continue
        f=np.concatenate([z[u['surface']]if u['side']==1 else z[u['surface']][:,[0,2,1]]for u in uses])
        used,inv=np.unique(f,return_inverse=True);p=v[used];f=inv.reshape(-1,3)
        if name.startswith('white-'):original[name]=bank.source(p,f,'native:'+name)
        else:solids[name]=bank.literal_cell(p,f,'native:'+name)
        attributes[name]=meta['attributes'][name]
    assert len(original)==6 and len(solids)==7
    mz=np.load(parcels/'masks.npz');masks={}
    for name,uses in pm['cells'].items():
        f=mz[uses[0]['surface']];used,inv=np.unique(f,return_inverse=True)
        masks[name]=bank.source(mz['positions'][used],inv.reshape(-1,3),'tract-mask:'+name)
    before=sum(float(s.volume())for s in original.values());tract_volumes={};empty=[]
    # Parcels already have mutually exclusive ownership. Intersections all
    # start from the literal parent; subtraction is only for its remainder.
    for parent,solid in original.items():
        remainder=solid
        for name,mask in masks.items():
            # The atlas probability centroid can straddle the native median
            # plane by a fraction of a voxel. A named ipsilateral tract does
            # not acquire the opposite cord half from that sampling error.
            if name.rsplit('-',1)[1]!=parent.split('@')[0].rsplit('-',1)[1]:continue
            tract=name.rsplit('-',1)[0]
            funiculus=('posterior'if tract in ['gracile','cuneate']else
                'lateral'if tract in ['lateral-corticospinal','rubrospinal','anterior-spinocerebellar']else'anterior')
            if parent.split('-')[1]!=funiculus:continue
            q=solid^mask
            if q.is_empty():continue
            assert str(q.status())=='Error.NoError'
            cell='tract-'+name+'~'+parent
            solids[cell]=q;remainder=remainder-mask
            attributes[cell]=dict(meta['attributes'][parent],tract=name,parent_white_cell=parent,
                reference_status='Registered probabilistic territory; not donor histology')
            tract_volumes[name]=tract_volumes.get(name,0.)+float(q.volume())
        if not remainder.is_empty():solids[parent]=remainder
        else:empty.append(parent)
    assert all(str(s.status())=='Error.NoError'for s in solids.values()if hasattr(s,'status'))
    after=sum(float(s.volume())for n,s in solids.items()if n.startswith('white-')or n.startswith('tract-'))
    assert abs(before-after)<max(1e-15,before*1e-8),(before,after)
    bank.save(solids,out/'boolean-source.npz')
    atomic_json(out/'stage.json',dict(segment=segment,attributes=attributes,empty_remainders=empty,
        white_volume_before_m3=before,white_volume_after_m3=after,tract_volumes_m3=tract_volumes,
        source_native_sha256=digest_file(native/'cells.npz'),source_masks_sha256=digest_file(parcels/'masks.npz'),
        construction_only=True,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started))
    print(json.dumps(dict(segment=segment,cells=len(solids),tracts=len(tract_volumes),white_volume_difference_m3=after-before,
        elapsed_seconds=time.monotonic()-started)),flush=True)
else:
    source=Path(args.source).resolve();metadata=json.loads((source/'stage.json').read_text());assert metadata['segment']==segment
    bank,solids=BooleanCells.load(source/'boolean-source.npz')
    result=bank.assemble(solids,conform_collinear_cycles=True,
        progress=lambda n,total,overlays:print(json.dumps(dict(stage='reciprocal-charts',groups=n,total=total,overlays=overlays)),flush=True))
    checks={}
    for name,uses in result['cells'].items():
        f=np.concatenate([result['surfaces'][u['surface']]if u['side']==1 else result['surfaces'][u['surface']][:,[0,2,1]]for u in uses])
        t=topology(Mesh(result['positions'],f,'model:local-female-ta98:assembled-metres'))
        checks[name]={k:t[k]for k in ['closed_oriented_manifold','open_edges','nonmanifold_edges','invalid_vertex_links','zero_area_faces','components']}
    failed=[n for n,t in checks.items()if not t['closed_oriented_manifold']]
    atomic_npz(out/'cells.npz',positions=result['positions'],**result['surfaces'])
    atomic_json(out/'cells.json',dict(cells=result['cells'],surfaces={k:len(f)for k,f in result['surfaces'].items()},attributes=metadata['attributes'],
        construction={**result['evidence'],'source_native_sha256':metadata['source_native_sha256'],
            'source_masks_sha256':metadata['source_masks_sha256'],'native_envelope_and_grey_replayed':True,
            'probabilistic_reference_regions':True,'not_donor_histological_boundaries':True},
        topology=checks,segment=segment,passed=not failed,installed=False,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-started))
    assert not failed,failed
    print(json.dumps(dict(segment=segment,passed=True,cells=len(result['cells']),interfaces=len(result['surfaces']),
        triangles=sum(len(f)for f in result['surfaces'].values()),elapsed_seconds=time.monotonic()-started)),flush=True)
