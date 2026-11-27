"""Save an assembled actual-bank surface/section review spec; never queue it."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file


def main(bank,intrinsics,output):
    bank,intrinsics,output=(p.resolve()for p in (bank,intrinsics,output))
    if output.exists():
        raise FileExistsError(output)
    b=json.loads(bank.read_text());c=json.loads(intrinsics.read_text())
    if not b['passed']or set(b['parts'])!={'mucosa','aponeurosis','superior-longitudinal'}:
        raise ValueError('Use the actual passed three-cell bank')
    if digest_file(b['bank'])!=b['bank_sha256']:
        raise ValueError('Actual bank changed')
    colors={'mucosa':[.76,.39,.40,1.],'aponeurosis':[.76,.66,.61,1.],
            'superior-longitudinal':[.72,.30,.29,1.]}
    entries=[]
    for key,record in b['parts'].items():
        path=bank.with_name(record['id'].replace(':','_')+'.npz')
        entries.append(dict(mesh=str(path),part_id=record['id'],name=key,color=colors[key],
                            tissue_draw_groups=record['tissueDrawGroups']))
    for row in c['meshes']:
        if row['key']=='superior-longitudinal':
            continue
        path=Path(row['mesh'])
        if digest_file(path)!=row['sha256']:
            raise ValueError('Literal intrinsic source changed')
        entries.append(dict(mesh=str(path),name=row['key'],color=row['color'],
            part_id='TA98TONGUE:male:intrinsic:'+row['key'],corner_normal_source=row['corner_normal_source']))
    spec=dict(frame=b['frame'],assembled_tissues=True,meshes=entries,resolution=800,wall_seconds=60,
        views=[dict(name='front',direction=[0.,0.,1.]),dict(name='right',direction=[1.,0.,0.]),
               dict(name='dorsal',direction=[0.,1.,0.]),
               dict(name='sagittal',direction=[1.,0.,0.],
                    section=dict(normal=[1.,0.,0.],point=[0.,1.63,.055],clear_outer=True))],
        source_bank=str(bank),source_bank_sha256=digest_file(bank),
        source_intrinsics=str(intrinsics),source_intrinsics_sha256=digest_file(intrinsics),
        anatomical_acceptance=False)
    atomic_json(output,spec)
    print(json.dumps(dict(output=str(output),sha256=digest_file(output),actualMeshes=len(entries),
                          assembled_tissues=True,queued=False)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bank',type=Path,required=True);p.add_argument('--intrinsics',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.bank,a.intrinsics,a.output)
