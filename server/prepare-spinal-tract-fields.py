"""Register probabilistic tract fields, with explicit anatomical level limits.

These are reference-guided territories, not donor histology or axon models.
The fields are construction classifiers. Only their disjoint intersection
with the existing finite native white matter may become displayed tissue.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True)
ap.add_argument('--reticular',action='store_true',help='Include source-labelled medial, lateral and ventral reticulospinal fields with explicit TA98 parentage')
args=ap.parse_args();started=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
reference=base/'pam50-white-matter-v1';template=base/'pam50-reference-v1'
registration=json.loads((base/'female-root-regions-v2/root-regions.json').read_text())
segments=registration['segments'];level_names=[s['name']for s in segments]
assert level_names==[*(f'C{i}'for i in range(1,9)),*(f'T{i}'for i in range(1,13)),*(f'L{i}'for i in range(1,6)),*(f'S{i}'for i in range(1,6)),'Co1']
# Only direct atlas identities are used. In particular, a combined spinal
# lemniscus field is not renamed as two different anatomical tracts.
specs=[
    ('gracile','A14.1.02.243',0,'C1','S1'),
    ('cuneate','A14.1.02.244',2,'C1','T6'),
    ('lateral-corticospinal','A14.1.02.219',4,'C1','S1'),
    ('anterior-spinocerebellar','A14.1.02.226',6,'C1','S1'),
    ('rubrospinal','A14.1.02.220',8,'C1','C8'),
    ('spino-olivary','A14.1.02.230',14,'C1','S1'),
    ('lateral-vestibulospinal','A14.1.02.206',18,'C1','S1'),
    ('anterior-corticospinal','A14.1.02.205',22,'C1','L3'),
    ('tectospinal','A14.1.02.211',24,'C1','C4'),
]
if args.reticular:
    specs.extend([
        ('bulboreticulospinal','A14.1.02.221',10,'C1','S1'),
        ('reticulospinal-fibres','A14.1.02.208',20,'C1','S1'),
        ('pontoreticulospinal','A14.1.02.209',26,'C1','S1'),
    ])
authority=json.loads((ROOT/'reports/ta98-spinal-cord-inventory.json').read_text())
terms={r['authority']['id']:r['authority']for r in authority['terms']}
plan=[]
for name,term,number,first,last in specs:
    parents=[];parent=terms[term]['parent_id']
    while parent in terms and parent not in ['A14.1.02.202','A14.1.02.215','A14.1.02.239']:
        parents.append(parent);parent=terms[parent]['parent_id']
    funiculus={'A14.1.02.202':'anterior','A14.1.02.215':'lateral','A14.1.02.239':'posterior'}.get(parent)
    assert funiculus,'Every tract needs its authoritative funiculus'
    for side,offset in [('left',0),('right',1)]:
        plan.append(dict(name=name+'-'+side,term=term,authority=terms[term],atlas_number=number+offset,
            first_level=first,last_level=last,representation='Finite longitudinal white-matter territory',
            anatomical_limits_estimated=True,termination_is_not_a_reconstructed_axon_terminal=True,
            parent_terms=parents,funiculus=funiculus,
            source_identity_note=('PAM50 lateral reticulospinal maps to medullary/bulbar; medial reticulospinal maps to pontine. The ventral field remains the collective reticulospinal identity.'if name in ['bulboreticulospinal','pontoreticulospinal','reticulospinal-fibres']else'Direct source atlas identity'),
            parent='native white matter, with original funiculus ownership retained'))
atomic_json(out/'plan.json',dict(schema='native-spinal-tract-plan/v1',parts=plan,
    fixed=['All skeleton geometry','Native cord outer surface','Native grey matter','Pia','Central canal'],
    interfaces=['Every tract/remainder interface must be reciprocal','Tracts must not overlap grey matter or each other',
                'Segment boundaries remain shared','No masks are anatomy before the native intersection'],
    references=['https://spinalcordtoolbox.com/stable/overview/concepts/pam50.html',
        'https://nba.uth.tmc.edu/NEUROANATOMY/L2/Lab02p11_index.html',
        'https://nba.uth.tmc.edu/NEUROANATOMY/L6/Lab06p20_index.html',
        'https://nba.uth.tmc.edu/Assets/pdf/courses/ns_2013_lab_guide.v13.1.0_SCR.pdf',
        'https://learninglink.oup.com/access/content/blumenfeld3e-student-resources/blumenfeld3e-interactive-table-6-3-lateral-and-medial-descending-motor-systems'],
    atlas_reference_only=True,redistribution_rights_not_verified=True,anatomical_acceptance=False,installed=False))

cord=nib.load(template/'template/PAM50_cord.nii.gz');cord_data=np.asarray(cord.dataobj)
level_img=nib.load(template/'template/PAM50_spinal_levels.nii.gz');level_data=np.asarray(level_img.dataobj)
assert nib.aff2axcodes(cord.affine)==('L','A','S')
axis=np.linspace(-1.08,1.08,25);uv=np.stack(np.meshgrid(axis,axis,indexing='ij'),axis=-1).reshape(-1,2)
sections=[]
for number in range(1,31):
    zs=np.flatnonzero(np.any(level_data==number,axis=(0,1)));z=int(round(float(zs.mean())))
    ij=np.argwhere(cord_data[:,:,z]>0);low,high=ij.min(0)-.5,ij.max(0)+.5
    sections.append(dict(name=level_names[number-1],z=z,center=(low+high)/2,radii=(high-low)/2))
scores=np.zeros((30,len(uv),len(plan)),np.float64);files=[]
for k,row in enumerate(plan):
    path=reference/f"atlas/PAM50_atlas_{row['atlas_number']:02}.nii.gz";image=nib.load(path)
    assert image.shape==cord.shape and np.array_equal(image.affine,cord.affine)
    data=np.asarray(image.dataobj)
    for i,s in enumerate(sections):
        if not level_names.index(row['first_level'])<=i<=level_names.index(row['last_level']):continue
        xy=s['center']+uv*s['radii'];prob=map_coordinates(data[:,:,s['z']],xy.T,order=1,mode='constant',cval=0.)
        scores[i,:,k]=np.clip(prob,0.,1.)
    files.append(dict(name=row['name'],source_sha256=digest_file(path)))
atomic_npz(out/'fields.npz',normalized_controls=uv,axis=axis,scores=scores,
    segment_y=np.array([s['y_range_m']for s in segments]),template_z=np.array([s['z']for s in sections]))
atomic_json(out/'fields.json',dict(passed=True,construction_only=True,names=[r['name']for r in plan],
    score_axis='left, anterior; native +X,+Z',coordinate_grid=25,source_files=files,
    reference_commit=json.loads((reference/'manifest.json').read_text())['commit'],
    plan_sha256=digest_file(out/'plan.json'),probability_floor=.035,
    missing_sacral_samples=['S2','S3','S4','S5','Co1'],missing_samples_extrapolated=False,
    scope='Probabilistic template parcels; registration and anatomy still require review',
    elapsed_seconds=time.monotonic()-started,anatomical_acceptance=False,installed=False))
print(json.dumps(dict(passed=True,tract_identities=len(specs),bilateral_fields=len(plan),levels=30,
    elapsed_seconds=time.monotonic()-started)),flush=True)
