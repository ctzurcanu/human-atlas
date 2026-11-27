"""Partition the fitted shared tissues into the 31 registered cord segments.

All cuts use one position bank and reciprocal caps. Existing macro-region
cuts are removed before segmentation, including the pia. Segment boundaries
are root-registered estimates, not claims of 31 additional TA98 terms.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.plane_tissue_partition import partition
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion';registration_path=base/'female-root-regions-v2/root-regions.json'
registration=json.loads(registration_path.read_text());d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz')
owners=defaultdict(list)
for name,uu in d['cells'].items():
    for u in uu:owners[u['surface']].append((name,u))
cells=defaultdict(list);surfaces={};removed=0
for sid,oo in owners.items():
    names=[n.split('@')[0]for n,_ in oo]
    if len(oo)==2 and names[0]==names[1]:removed+=1;continue
    surfaces[sid]=z[sid]
    for (name,u),macro in zip(oo,names):
        peer=u.get('neighbor');peer=peer.split('@')[0]if peer else None
        cells[macro].append(dict(surface=sid,side=u['side'],neighbor=peer))
assert len(cells)==13 and 'pia'in cells
bands=[(r['name'],*r['y_range_m'])for r in registration['segments'][::-1]]
print(json.dumps(dict(stage='shared-segment-cuts',segments=len(bands),macro_tissues=len(cells),removed_internal_charts=removed)),flush=True)
result=partition(z['positions'],surfaces,cells,bands,origin=[0.,1.35,-.05])
attributes={};proof={}
for i,(name,uu)in enumerate(result['cells'].items()):
    macro,segment=name.split('@');band=next(b for b in bands if b[0]==segment)
    attributes[name]=dict(native_region=macro,axial_region=segment,y_range_m=list(band[1:]),segment_boundary_status=next(r['boundary_status']for r in registration['segments']if r['name']==segment))
    faces=np.concatenate([result['surfaces'][u['surface']]if u['side']==1 else result['surfaces'][u['surface']][:,[0,2,1]]for u in uu])
    t=topology(Mesh(result['positions'],faces,'model:local-female-ta98:assembled-metres'))
    proof[name]={k:t[k]for k in ['closed_oriented_manifold','open_edges','nonmanifold_edges','invalid_vertex_links','zero_area_faces']}
    if(i+1)%50==0:print(json.dumps(dict(stage='topology',cells=i+1,total=len(result['cells']),elapsed_seconds=time.monotonic()-start)),flush=True)
bad=[n for n,t in proof.items()if not t['closed_oriented_manifold']]
atomic_npz(out/'cells.npz',positions=result['positions'],**result['surfaces'])
atomic_json(out/'cells.json',dict(cells=result['cells'],surfaces={k:len(f)for k,f in result['surfaces'].items()},attributes=attributes,bands=bands,
    construction={**d['construction'],'source_cell_bank_sha256':digest_file(source/'cells.npz'),'registered_segments':31,'segment_partition':result['evidence'],
        'regional_boundaries_are_estimates':True,'segment_names_are_not_new_ta98_terms':True},topology=proof,registration_sha256=digest_file(registration_path),
    anatomical_acceptance=False,installed=False,passed=not bad,elapsed_seconds=time.monotonic()-start))
assert not bad,bad
print(json.dumps(dict(passed=True,cells=len(proof),segments=31,triangles=sum(len(f)for f in result['surfaces'].values()),elapsed_seconds=time.monotonic()-start)),flush=True)
