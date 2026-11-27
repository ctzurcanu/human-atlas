"""Local-only named landmark inventory; exports no reference mesh arrays.

Language entries are distinguished from actual structures, patches and nodes.
The complete local metadata set is searched once to find usable references
missed by the current body catalogue. Coarse shape relations are deliberately
rounded to 5% of the reference whole tongue, not a registration or a template.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import struct
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,digest_file
EXT=Path('/Users/christiantzurcanu/Documents/dev/3d/ext')
WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'
TERMS=[
 dict(term='A05.1.04.002',name='body of tongue',match=[]),
 dict(term='A05.1.04.003',name='root of tongue',match=[]),
 dict(term='A05.1.04.008',name='fimbriated fold',match=[]),
 dict(term='A05.1.04.012',name='frenulum of tongue',match=[]),
 dict(term='A05.1.04.018',name='midline groove of tongue',match=[]),
 dict(term='A05.1.04.019',name='terminal sulcus of tongue',match=[]),
 dict(term='A05.1.04.020',name='foramen caecum of tongue',match=[]),
 dict(term='A05.1.04.022',name='lingual tonsil',match=[26125]),
 dict(term='A05.1.04.023',name='lymphoid nodules',match=[]),
 dict(term='A05.1.04.016',name='vallate papillae',match=[15101,26176]),
 dict(term='A05.1.04.017',name='foliate papillae',match=[15095])]
RELATED={19491,15183,26279,26125,16127,16001,15101,26176,15095,15165}


def header(path):
    with path.open('rb')as f:
        h=f.read(20);magic,version,total,length,kind=struct.unpack('<4sIIII',h)
        if magic!=b'glTF'or version!=2 or kind!=0x4e4f534a:raise ValueError('Invalid GLB')
        raw=f.read(length);g=json.loads(raw)
    return g,hashlib.sha256(raw).hexdigest()


def main(output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local stage')
    language_path=EXT/'meshes/en.json';language=json.loads(language_path.read_text())
    names={r['cid']:r for r in language['items']}
    metadata=[];records=[]
    for path in sorted((EXT/'meshes').glob('*/assets.json')):
        data=json.loads(path.read_text());found=[]
        for collection in ['structures','clusters','shapes','spots','compositions']:
            for index,row in enumerate(data.get(collection,[])):
                if not isinstance(row,dict)or row.get('cid')not in RELATED:continue
                value=dict(asset=path.parent.name,collection=collection,index=index,
                    cid=row['cid'],name=names[row['cid']]['name'],source_id=row.get('id'),
                    nodes=[s['id']for s in row.get('sni',[])],
                    linked_structure_indices=row.get('linked',[]),
                    parent_or_child_group_only=collection=='clusters')
                found.append(value);records.append(value)
        if found:
            metadata.append(dict(path=str(path),sha256=digest_file(path),matching_records=found,
                description=json.loads((path.parent/'model_description.json').read_text())if(path.parent/'model_description.json').exists()else None))
    models=[];model_json={}
    for relative in ['models/Male_2024.glb','models/Female_2024.glb',
                     'models2/Male_2024.glb','models2/Female_2024.glb',
                     'models5/Male_2022.glb','models5/Tongue.glb']:
        path=EXT/relative
        if not path.exists():continue
        g,sha=header(path);by={n.get('name'):i for i,n in enumerate(g['nodes'])};items=[]
        for row in records:
            for name in row['nodes']:
                if name not in by:continue
                node=g['nodes'][by[name]]
                if 'mesh'not in node:continue
                mesh=g['meshes'][node['mesh']]
                items.append(dict(cid=row['cid'],name=row['name'],node_name=name,node_index=by[name],
                    mesh_index=node['mesh'],source_metadata_asset=row['asset'],
                    primitives=[dict(mode=p.get('mode',4),position_accessor=p['attributes']['POSITION'],
                        positions=g['accessors'][p['attributes']['POSITION']]['count'],
                        index_accessor=p.get('indices'),indices=g['accessors'][p['indices']]['count']if'indices'in p else None)
                        for p in mesh['primitives']]))
        # Duplicate metadata versions can point to one literal source node.
        unique={r['node_name']:r for r in items}
        models.append(dict(path=str(path),file_bytes=path.stat().st_size,json_sha256=sha,
            asset=g.get('asset'),matched_geometry=list(unique.values()),
            no_named_geometry_for_terms=[r['term']for r in TERMS if not any(x['cid']in r['match']for x in unique.values())]))
        model_json[relative]=g
    # Only the original local reader performs the selected in-memory mapping.
    spec=importlib.util.spec_from_file_location('tongue_reference_reader',ROOT/'server/analyse-tongue-reference-shapes.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    coarse=[]
    for relative in ['models5/Male_2022.glb']:
        path=EXT/relative;reader=module.SelectedGlb(path)
        target={r['cid']:r for r in next(x for x in models if x['path']==str(path))['matched_geometry']}
        if not all(k in target for k in [26279,26125,26176]):
            coarse.append(dict(model=str(path),supported=False,reason='Named posterior landmark geometry absent from this GLB'))
            continue
        whole,_,_=reader.mesh(target[26279]['node_name']);lower=whole.min(0);extent=np.ptp(whole,axis=0)
        local=[];centers={}
        for cid in [26125,26176]:
            v,f,_=reader.mesh(target[cid]['node_name']);used=v[np.unique(f)];centers[cid]=used.mean(0)
            # No exact reference position, triangle, UV or normal is written.
            bounds=(np.array([used.min(0),used.max(0)])-lower)/extent
            local.append(dict(cid=cid,name=target[cid]['name'],triangles=len(f),
                normalized_whole_bounds_in_5_percent_steps=(np.round(bounds/.05)*.05).tolist(),
                coarse_only=True))
        coarse.append(dict(model=str(path),supported=True,world_axes=['left-right','superior-inferior','anterior-posterior'],
            items=local,tonsil_center_posterior_to_vallate=bool(centers[26125][2]<centers[26176][2]),
            reference_registration_to_current_atlas=False,
            transfer_exact_curve_or_profile_permitted=False))
    declarations=[]
    p=EXT/'meshes/07cf4f8f5d8b76282917320715dda2ad/assets.gltf'
    g=json.loads(p.read_text());by={n.get('name'):i for i,n in enumerate(g['nodes'])}
    ids={26125:'71a7ce8f0b3da33646612ec1a36f6e3c',
         26176:'5682d31c4999a1c2c21d8eca20076561',
         26279:'5daf1ee2865580503f26b52edebb97d1'}
    bounds={}
    for cid,name in ids.items():
        node=g['nodes'][by[name]];rows=[]
        for primitive in g['meshes'][node['mesh']]['primitives']:
            a=g['accessors'][primitive['attributes']['POSITION']]
            rows.append(dict(position_accessor=primitive['attributes']['POSITION'],
                positions_declared=a['count'],indices_declared=g['accessors'][primitive['indices']]['count'],
                primitive_mode=primitive.get('mode',4)))
            bounds[cid]=np.asarray([a['min'],a['max']])
        declarations.append(dict(cid=cid,name=names[cid]['name'],node_name=name,
            node_index=by[name],mesh_index=node['mesh'],primitives=rows,
            source_gltf=dict(path=str(p),sha256=digest_file(p)),
            vertex_index_binary_available=False,actual_surface_inspected=False))
    lo=bounds[26279][0];extent=np.ptp(bounds[26279],axis=0)
    declared_coarse=[dict(cid=cid,name=names[cid]['name'],
        normalized_whole_accessor_bounds_in_5_percent_steps=(np.round((b-lo)/extent/.05)*.05).tolist())
        for cid,b in bounds.items()if cid!=26279]
    missing_buffers=[dict(uri=r['uri'],declared_bytes=r['byteLength'],
        exists_at_declared_path=(p.parent/r['uri']).exists())for r in g['buffers'][:2]]
    summaries=[]
    for row in TERMS:
        item=dict(row);item.pop('match')
        item['metadata_records']=[r for r in records if r['cid']in row['match']]
        item['actual_models']=[dict(model=m['path'],records=[r for r in m['matched_geometry']if r['cid']in row['match']])
            for m in models if any(r['cid']in row['match']for r in m['matched_geometry'])]
        item['named_reference_geometry_found']=bool(item['actual_models'])
        item['gltf_geometry_declarations']=[r for r in declarations if r['cid']in row['match']]
        item['named_geometry_declared_but_binary_unavailable']=bool(item['gltf_geometry_declarations'])and not item['named_reference_geometry_found']
        item['missing_named_reference_does_not_prove_absent_from_carrier_texture']=not item['named_reference_geometry_found']
        summaries.append(item)
    output.mkdir(parents=True)
    report=dict(schema='tongue-macroscopic-reference-inventory/v1',
        language=dict(path=str(language_path),sha256=digest_file(language_path)),
        metadata=metadata,models=models,terms=summaries,coarse_posterior_relations=coarse,
        older_male_posterior_geometry_declarations=dict(items=declarations,buffers=missing_buffers,
            normalized_accessor_bounds=declared_coarse,
            tonsil_accessor_bound_center_posterior_to_vallate=bool(bounds[26125][:,2].mean()<bounds[26176][:,2].mean()),
            source_parent_transforms_and_surface_shape_not_reconstructed=True,
            declaration_not_actual_surface_review=True),
        description_paraphrases=dict(
            mucosa='Supplied description distinguishes dorsal and ventral mucosa, describes the median and terminal grooves, posterior tonsillar covering and glossoepiglottic continuity. This is text evidence, not a separately segmented groove.',
            tongue='Supplied description identifies a V-shaped anterior/posterior division and a posterior apical blind pit; the named main 2024 body geometry contains no separate landmark identifiers.',
            tonsil='Supplied definition describes collective lymphoid nodules at the tongue root; no individual nodule source IDs were found.'),
        primary_anatomical_constraints=[
            dict(terms=['A05.1.04.012'],reference='https://onlinelibrary.wiley.com/doi/10.1002/ca.23343',
                source='Nine fresh adult cadaver dissections',
                constraints=['Frenulum is a mobilized midline fold of the floor-of-mouth fascia and overlying oral mucosa, not a discrete cord.',
                    'Fascia joins the inner mandibular periosteum broadly and the anterior/ventral tongue connective tissue.',
                    'GG and sublingual structures lie beneath/suspended from this fascia.'],
                construction_implication='Establish actual ventral tongue/floor-of-mouth carrier and mandibular rim first; author a continuous finite fold, preserve glands/ducts/lingual nerves.'),
            dict(terms=['A05.1.04.012'],reference='https://onlinelibrary.wiley.com/doi/full/10.1155/2020/1820978',
                source='Histology of thirteen adult cadavers',
                constraints=['Fascial thickness and incorporation into the fold vary; no universal thin midline band is justified.'],
                construction_implication='Do not infer a donor-specific thickness or a midline GG partition from a schematic.'),
            dict(terms=['A05.1.04.008'],reference='https://ifaa.unifr.ch/Public/TNAEntryPage/auto/unit/LAES/TAH2410%20Unit%20EN.htm',
                source='FIPAT/IFAA taxonomic definition',
                constraints=['Fimbriated fold is a lobulated zone of the inferior surface extending posteriorly and laterally from the anterior frenulum extremity.'],
                construction_implication='First locate the actual frenulum extremity; no arbitrary paired ventral ribbons.'),
            dict(terms=['A05.1.04.020'],reference='https://pubmed.ncbi.nlm.nih.gov/1746229/',
                source='Human adult cadaver histology, 35 tongues',
                constraints=['A foramen caecum was observed in 18 of 35 adult specimens; presence and depth vary.'],
                construction_implication='A supported blind surface pit is a possible ordinary adult representative, but not a mandatory universal deep tube or a separate solid.'),
            dict(terms=['A05.1.04.022','A05.1.04.023'],reference='https://pubmed.ncbi.nlm.nih.gov/1634989/',
                source='Human lingual-tonsil gross/histological autopsy investigation',
                constraints=['Gross tonsillar prominence and lymphoid/crypt tissue organization vary among adults and with age.'],
                construction_implication='Build a qualified finite posterior lymphoid field with actual mucosal crypt invaginations; not repeated disconnected balls or a universal hypertrophic bulk.')],
        next_evidence_actions=[
            'Inspect the old male whole/tonsil/vallate qualitative arrangement locally and establish independently authored coarse posterior landmark zones; current main 2024 models lack these named nodes.',
            'Use the posterior V arrangement as a reference constraint, never transfer exact mesh coordinates/curves. Locate terminal groove relative to current PG entry, actual dorsal slope and lingual tonsil/epiglottic context.',
            'Identify actual current floor-of-mouth mucosa and sublingual duct/nerve interfaces before building frenulum/fimbriated folds.',
            'Keep body/root volume division pending until terminal landmark and adjoining posterior interfaces are established.'],
        rights=dict(local_reference_only=True,redistribution_license_established=False,
            vendor_position_normal_index_uv_texture_arrays_exported=False,
            exact_vendor_dimensions_or_boundary_curves_exported=False),
        installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'inventory.json',report)
    print(json.dumps(dict(report=str(output/'inventory.json'),metadata_catalogues=len(metadata),
        terms=[dict(term=r['term'],name=r['name'],actual_reference_geometry=r['named_reference_geometry_found'])for r in summaries],
        coarse_posterior_relations=coarse,elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();main(a.output.resolve())
