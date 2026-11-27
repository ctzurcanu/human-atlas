"""Source-bound mucosa/aponeurosis cells with one reciprocal interface bank."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.tissue import layer_arrangement,normal_words,pack
from geometry_tools.audit import topology
from tongue_whole_input import whole_input

WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
REFERENCE='https://onlinelibrary.wiley.com/doi/10.1111/joa.13637'


def main(output,whole=None,fit_result=None):
    began=time.monotonic();output=output.resolve()
    if output.exists() or not output.is_relative_to(WORK) or output==WORK:
        raise ValueError('Use a new local tongue stage')
    m,source,binding=whole_input(whole,fit_result)
    if digest_file(m.source['catalogue_path'])!=m.source['catalogue_sha256']:
        raise ValueError('Stale source catalogue')
    v,f,n=m.vertices,m.faces,m.attributes['normals'].astype(float)
    n/=np.linalg.norm(n,axis=1)[:,None]
    thickness=.00035+.00007*np.clip(n[:,1],0,1)-.00005*np.clip(-n[:,1],0,1)
    inner=v-thickness[:,None]*n
    tri=v[f];centres=tri.mean(1);fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);fn/=np.linalg.norm(fn,axis=1)[:,None]
    # Aponeurosis is a dorsal plate, fading at the margins and behind the
    # terminal sulcus. The V is an explicit approximate source landmark;
    # it is not a newly established donor sulcus or a ventral whole shell.
    z0,z1=v[:,2].min(),v[:,2].max()
    terminal=z1-(z1-z0)*2/3+.4*np.abs(centres[:,0])
    mask=(fn[:,1]>.45)&(centres[:,2]>terminal)
    selected=np.flatnonzero(mask);patch=f[selected]
    patch_topology=topology(Mesh(v,patch,m.frame))
    if len(patch_topology['components'])!=1:
        raise ValueError('Dorsal aponeurosis patch must be one connected region')
    # Primary histology bounds: central 280–320 µm; lateral 100–150 µm.
    # The explicit lateral taper is illustrative, not an individual measure.
    lateral=np.clip(np.abs(v[:,0])/.016,0,1)
    apo_thickness=.00030-.000175*lateral*lateral*(3-2*lateral)
    deep=inner-apo_thickness[:,None]*n
    output.mkdir(parents=True)
    plan=dict(schema='tongue-finite-coverings/v1',source=str(source),source_sha256=digest_file(source),whole_input=binding,
        source_geometry_sha256=m.geometry_sha,source_catalogue_sha256=m.source['catalogue_sha256'],
        terms=['A05.1.04.011','A05.1.04.025'],reference=REFERENCE,
        reference_scope='Human central/lateral/inferior histology; a dorsal aponeurotic plate rather than a whole-organ shell',
        modeling_field_not_donor_measurement=True,aponeurosis_source_faces=selected.tolist(),
        approximate_terminal_boundary=True,aponeurosis_not_ventral=True,
        mucosa_range_m=[float(thickness.min()),float(thickness.max())],
        aponeurosis_range_m=[float(apo_thickness[np.unique(patch)].min()),float(apo_thickness[np.unique(patch)].max())],
        fixed=['All source coordinates and source skeleton'],
        adjoining_interfaces=[dict(first='mucosa',second='aponeurosis',kind='literal-opposed-shared-chart'),
            dict(first='aponeurosis',second='underlying-muscle-region',kind='unjoined; final muscle placement check pending')],
        source_geometry_modified=False,anatomical_acceptance=False,installed=False,wall_seconds=40)
    atomic_json(output/'construction-input.json',plan)
    count=len(v)
    _,apo_surfaces,apo_cells,prism=layer_arrangement(np.array([inner,deep]),patch,['aponeurosis'])
    surfaces={'mucosa:outer':f,'mucosa:inner-rest':f[~mask]+count,'mucosa:aponeurosis':patch+count,'aponeurosis:inner':patch+2*count}
    cells={'mucosa':[dict(surface='mucosa:outer',side=1,neighbor=None),dict(surface='mucosa:inner-rest',side=-1,neighbor=None),dict(surface='mucosa:aponeurosis',side=-1,neighbor='aponeurosis')],
        'aponeurosis':[dict(surface='mucosa:aponeurosis',side=1,neighbor='mucosa'),dict(surface='aponeurosis:inner',side=-1,neighbor=None)]}
    for name,faces in apo_surfaces.items():
        if name.startswith('rim:'):
            key='aponeurosis:'+name;surfaces[key]=faces+count
            cells['aponeurosis'].append(dict(surface=key,side=1,neighbor=None))
    positions=np.r_[v,inner,deep]
    words=normal_words(positions,surfaces)
    words['mucosa:outer']=m.attributes['normals'][f].copy()
    result=pack(positions,surfaces,cells,{'mucosa':'TA98TONGUE:male:mucosa','aponeurosis':'TA98TONGUE:male:aponeurosis'},
        m.frame,output/'packed',origin=(0,1.5,0),corner_words=words,construction=prism,max_pairs=250000,max_tests=2000000)
    if time.monotonic()-began>=40:raise TimeoutError('Coverings exceeded 40-second worker budget')
    q=np.r_[inner,deep[np.unique(patch)],deep[patch].mean(1)]
    winding=np.abs(igl.winding_number(v,f,q));contained=bool(np.all(winding>.99))
    report=dict(**plan,passed=bool(result['checks']['passed'] and contained),packed_manifest=str(output/'packed/tissue.json'),
        actual_emitted_checks=result['checks'],aponeurosis_patch_topology=patch_topology,inner_samples_contained=contained,
        shared_interface_stored_once=True,reciprocal_normals_same_words_opposite_sign=True,total_bytes=result['total_bytes'],new_image_assets=0,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'construction.json',report)
    print({key:report[key] for key in ['passed','total_bytes','inner_samples_contained','elapsed_seconds']},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--whole-tongue',type=Path);parser.add_argument('--fit-result',type=Path)
    args=parser.parse_args();main(args.output,args.whole_tongue,args.fit_result)
