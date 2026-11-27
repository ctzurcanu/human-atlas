"""Small dimensionless spatial priors, not copied atlas surface geometry.

Per-level probability moments guide original regional construction. A soft
atlas mask does not establish a tract's anatomical longitudinal termination.
Native grey/white boundaries and accepted level ranges still constrain use.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import nibabel as nib
from geometry_tools.common import atomic_json,digest_file

ap=argparse.ArgumentParser();ap.add_argument('--reference',required=True);ap.add_argument('--template',required=True);ap.add_argument('--output',required=True);ap.add_argument('--begin',type=int,default=0);ap.add_argument('--count',type=int,default=4);args=ap.parse_args();start=time.monotonic()
reference=Path(args.reference).resolve();template=Path(args.template).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
manifest=json.loads((reference/'manifest.json').read_text());cord_image=nib.load(template/'template/PAM50_cord.nii.gz');levels_image=nib.load(template/'template/PAM50_spinal_levels.nii.gz')
assert nib.aff2axcodes(cord_image.affine)==('L','A','S');assert np.array_equal(cord_image.affine,levels_image.affine)
cord=np.asarray(cord_image.dataobj);levels=np.asarray(levels_image.dataobj);samples=[]
for label in range(1,31):
    zs=np.flatnonzero(np.any(levels==label,axis=(0,1)));z=int(round(float(zs.mean())));mask=cord[:,:,z]>0;ij=np.argwhere(mask)
    low,high=ij.min(0)-.5,ij.max(0)+.5;center=(low+high)/2;radius=(high-low)/2
    uv=(ij-center)/radius;name='C'+str(label)if label<=8 else'T'+str(label-8)if label<=20 else'L'+str(label-20)if label<=25 else'S'+str(label-25)
    samples.append(dict(level=name,z=z,ij=ij,uv=uv,template_center_voxels=center.tolist(),template_radii_voxels=radius.tolist()))
names={}
for line in(reference/'atlas/info_label.txt').read_text().splitlines():
    fields=[s.strip()for s in line.split(',')]
    if len(fields)==3 and fields[0].isdigit():names[int(fields[0])]=fields[1]
rows=[]
for number in range(args.begin,min(30,args.begin+args.count)):
    path=reference/('atlas/PAM50_atlas_'+str(number).zfill(2)+'.nii.gz');image=nib.load(path)
    assert image.shape==cord.shape and np.array_equal(image.affine,cord_image.affine)
    values=np.asarray(image.dataobj);records=[]
    for sample in samples:
        ij=sample['ij'];uv=sample['uv'];weights=np.asarray(values[ij[:,0],ij[:,1],sample['z']],float);assert np.all((weights>=0)&(weights<=1))
        mass=float(weights.sum());record={k:v for k,v in sample.items()if k not in['ij','uv']}
        if mass>1e-10:
            mean=(uv*weights[:,None]).sum(0)/mass;delta=uv-mean;cov=(delta.T*weights)@delta/mass
            record.update(probability_mass_voxels=mass,normalized_center=mean.tolist(),normalized_covariance=cov.tolist())
        else:record.update(probability_mass_voxels=0.,normalized_center=None,normalized_covariance=None)
        records.append(record)
    row=dict(source_atlas_number=number,name=names[number],source_file_sha256=digest_file(path),reference_commit=manifest['commit'],
        coordinate_axes=['left','anterior'],native_female_axis_correspondence=['+X','+Z'],coordinates_dimensionless=True,levels=records,
        anatomical_longitudinal_extent_not_certified=True,donor_specific_territory_not_established=True,geometry_copied=False,installed=False,anatomical_acceptance=False)
    atomic_json(out/('tract-'+str(number).zfill(2)+'.json'),row);rows.append(dict(number=number,name=names[number],supported_samples=sum(r['probability_mass_voxels']>0 for r in records)))
    print(json.dumps(dict(tract=number,name=names[number],sampled_levels=30,elapsed_seconds=time.monotonic()-start)),flush=True)
atomic_json(out/'priors.json',dict(reference_manifest_sha256=digest_file(reference/'manifest.json'),cord_template_sha256=digest_file(template/'template/PAM50_cord.nii.gz'),tracts=rows,
    scope='Spatial construction priors only; not installed anatomy or TA98 coverage',elapsed_seconds=time.monotonic()-start))
