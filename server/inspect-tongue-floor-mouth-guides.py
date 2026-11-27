"""Bound source guide observations only; no tissue or attachment is inferred."""
from pathlib import Path
import argparse, json, sys, time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
import numpy as np
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from tongue_whole_input import WORK, CATALOGUE

def bound(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':digest_file(p)}

def main(output):
    began=time.monotonic()
    if output.exists() or not output.is_relative_to(WORK):raise ValueError('Use a new local stage')
    base=WORK/'floor-of-mouth-contract-v1/contract.json';contract=json.loads(base.read_text())
    if digest_file(CATALOGUE)!=contract['source_catalogue']['sha256']:raise ValueError('Current catalogue changed: no source mixing')
    output.mkdir();atlas=json.loads(CATALOGUE.read_text());part=next(p for p in atlas['parts']if p['id']=='ZA:Gingiva')
    gingiva=AtlasReader().read_part(atlas,part,frame=contract['frame'],catalogue_path=CATALOGUE,catalogue_sha=digest_file(CATALOGUE))
    gp=output/'actual-gingiva.npz';gingiva.save(gp)
    loops=[]
    for row in contract['duct_rim_candidates']:
        p=Path(row['witness']['path'])
        if digest_file(p)!=row['witness']['sha256']:raise ValueError('Source duct guide changed')
        with np.load(p,allow_pickle=False)as raw:vertices=raw['positions'].copy();edges=raw['edges'].copy()
        adjacency={}
        for a,b in edges:
            adjacency.setdefault(int(a),[]).append(int(b));adjacency.setdefault(int(b),[]).append(int(a))
        if any(len(x)!=2 for x in adjacency.values()):raise ValueError('Boundary is not disjoint simple loops')
        remaining=set(adjacency);components=[]
        while remaining:
            loop=[min(remaining)];last=None
            while True:
                now=loop[-1];nxt=next(q for q in adjacency[now]if q!=last)
                if nxt==loop[0]:break
                if nxt in loop:raise ValueError('Boundary cycle repeats')
                loop.append(nxt);last=now
            remaining.difference_update(loop);components.append(loop)
        if len(components)!=2:raise ValueError('Expected exactly two source duct endpoints')
        # Infer only an endpoint vicinity from the actual course, not a built ostium.
        means=np.array([vertices[x].mean(0)for x in components])
        anterior=int(np.argmax(means[:,2]));medial=int(np.argmin(np.abs(means[:,0])))
        if anterior!=medial:raise ValueError('Anterior and medial endpoint rankings disagree')
        for i,indices in enumerate(components):
            pts=vertices[indices];label='anterior-medial-oral-candidate'if i==anterior else'posterior-lateral-gland-candidate'
            path=output/(row['part_id'].replace(':','_').replace('.','_')+'-'+label+'.npz')
            atomic_npz(path,positions=pts,source_position_indices=np.array(indices),frame=np.array(contract['frame']))
            loops.append(dict(part_id=row['part_id'],loop=bound(path),source=bound(p),vertices=len(indices),
                center_m=pts.mean(0).tolist(),bounds_m=np.array([pts.min(0),pts.max(0)]).tolist(),
                perimeter_m=float(np.linalg.norm(pts-np.roll(pts,1,axis=0),axis=1).sum()),
                classification=label,anatomical_ostium_accepted=False,finite_wall_built=False,
                qualification='Two source endpoint loops. Anterior + medial ordering supports oral vicinity; no current caruncle or through-mucosa ostium exists.'))
    reviewdir=ROOT/'.local-models/geometry-workbench/review/0e0855785f305d7323bc269ba646bdffbed2f31faa3429f9e9c94ea004a487b2'
    views=[]
    for i,name in enumerate(['superior-context','right-context','sagittal-context']):
        image=reviewdir/f'view-{i}.png';views.append(dict(name=name,**bound(image),
            inspected=True,diagnostic_section_uncapped=i==2))
    result=dict(schema='tongue-fom-source-guides/v1',contract=bound(base),source_catalogue=bound(CATALOGUE),
        gingiva=dict(part_id=part['id'],mesh=bound(gp),source_record=part,geometry_sha256=gingiva.geometry_sha,
            topology=topology(gingiva),bounds_m=gingiva.bounds.tolist(),
            interpretation='Actual upper/lower gingiva geometry; lingual lower margin and fascial insertion are not yet annotated.'),
        duct_endpoint_loops=loops,
        actual_blender_context_review=dict(request_id=reviewdir.name,views=views,original_objects_preserved=True,
            observed='The selected tongue lies above the bilateral sublingual glands and mylohyoid/geniohyoid support, inside the native mandibular arc. Side and sagittal views show the unconstructed FOM space beneath the tongue. Context is translucent/dim and not proof of tissue interfaces.'),
        next_geometry_required=[
            'Select the actual lower lingual gingival margin on the exported gingiva and pair it to inner mandibular cortex at the broad fascial insertion arc; store face/barycentric attachment witnesses.',
            'Join that witnessed horseshoe arc to the literal ventral tongue reflection loop with curved finite mucosal/fascial fields, retaining source glands, ducts, muscles and nerves as fixed exclusions.',
            'Provide paired finite duct-wall collars through the new field at the exported anterior-medial loops. Do not cap the two mouths. Caruncular relief belongs to the mucosal field around those channels.',
            'Use a broad variable-thickness connective roof, with central elevation/reflection producing the frenulum as a tissue fold; no independent cord or generic triangular slab.',
            'Replace the current tongue reflection closure only by explicit reciprocal charts and shared positions/normals; keep posterolateral boundaries open to continued oral mucosa until source attachment is defined.'],
        no_existing_mesh_modified=True,no_new_tissue_constructed=True,installed=False,anatomical_acceptance=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'guides.json',result)
    print(json.dumps({'output':str(output/'guides.json'),'gingiva_triangles':len(gingiva.faces),
        'duct_loops':[{k:r[k]for k in ['part_id','vertices','classification','center_m']}for r in loops],
        'elapsed_seconds':time.monotonic()-began}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();main(args.output.resolve())
