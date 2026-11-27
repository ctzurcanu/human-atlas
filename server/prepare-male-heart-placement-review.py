"""Stage the installed male heart, its native surface vessels and fixed thorax.

No fitting, source edits, new main Blender scene or catalogue installation.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh

SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='6526654a0e412947121518636d1b8004a934093f8a4cf459afcfa55c80a1b63c'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-heart-placement-source-v1'
VESSELS=['ZA:Left coronary artery','ZA:Right coronary artery','ZA:Anterior interventricular artery',
         'ZA:Circumflex artery of heart','ZA:Right inferolateral branch of right coronary artery',
         'ZA:Coronary sinus','ZA:Great cardiac vein','ZA:Middle cardiac vein','ZA:Inferior vein of left ventricle']


def main():
    raw=SOURCE.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    assert sha==EXPECTED,'Explicit source revision changed'
    assert not OUT.exists(),'Keep prior review evidence'
    atlas=json.loads(raw);parts={p['id']:p for p in atlas['parts']}
    OUT.mkdir();reader=AtlasReader(ROOT);entries={};records=[]
    def export(pid,color):
        if pid in entries:return entries[pid]
        part=parts[pid];assert not part.get('suppressed')
        mesh=reader.read_part(atlas,part,frame='model:local-ta98:assembled-metres',catalogue_path=SOURCE,catalogue_sha=sha)
        if pid=='TA98HEART:male:epicardium':
            # Diagnostic free outer chart only; all coordinates/attributes stay literal.
            faces=[]
            for group in part['tissueDrawGroups']:
                if not group.get('neighbor'):
                    assert group['start']%3==0 and group['count']%3==0
                    faces.extend(range(group['start']//3,(group['start']+group['count'])//3))
            mesh=Mesh(mesh.vertices,mesh.faces[faces],mesh.frame,mesh.attributes,mesh.source)
        path=OUT/f'{len(entries)}.npz';mesh.save(path)
        entry=dict(mesh=str(path),name=part['name'],color=color);entries[pid]=entry
        records.append(dict(part_id=pid,source_record=part,geometry_sha256=mesh.geometry_sha,
                            bounds_m=mesh.bounds.tolist(),source_bank_sha256=mesh.source['bank_sha256']))
        return entry
    heart=export('TA98HEART:male:epicardium',[.67,.34,.31,1.])
    heart_mesh=Mesh.load(heart['mesh']);points=heart_mesh.vertices[np.unique(heart_mesh.faces)]
    # Principal spans describe this mesh only, not clinical long-axis dimensions.
    _,basis=np.linalg.eigh(np.cov(points.T));principal=points@basis
    metrics=dict(world_axis_spans_mm=(np.ptp(points,axis=0)*1000).tolist(),
                 principal_spans_mm=(np.ptp(principal,axis=0)*1000).tolist(),
                 pca_sampling='unique actual free-chart source controls; descriptive, not anatomical axes')
    vessels=[export(pid,[.94,.22,.12,1.]if parts[pid]['system']=='arterial'else[.16,.39,.72,1.])for pid in VESSELS if pid in parts]
    context=[]
    for concept in atlas['concepts']:
        if concept['id'].startswith('ZA:')and 'lobe of'in concept['name'].lower()and concept['name'].lower().endswith('lung'):
            context.extend(export(pid,[.48,.58,.62,.17])for pid in concept['elements']if pid in parts and not parts[pid].get('suppressed'))
    context.extend(export(pid,[.83,.79,.68,.6])for pid,p in parts.items()if pid.startswith('ZA:')and'sternum'in p['name'].lower()and not p.get('suppressed'))
    views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[-1,0,0]),dict(name='upper-oblique',direction=[.7,.4,1])]
    for name,meshes in [('surface-vessels',[heart]+vessels),('fixed-thorax',[heart]+vessels+context)]:
        spec=dict(frame='model:local-ta98:assembled-metres',meshes=meshes,views=views,resolution=1100,wall_seconds=45)
        (OUT/f'review-{name}-spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    result=dict(schema='male-heart-placement-source-review/v1',source_catalogue_sha256=sha,records=records,
                heart_metrics=metrics,source_unchanged=True,heart_fitted=False,skeleton_fixed=True,
                anatomical_acceptance=False,blood_only_scope=True)
    (OUT/'source-evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==sha
    print(json.dumps(dict(output=str(OUT),meshes=len(records),heart_metrics=metrics,source_unchanged=True)))


if __name__=='__main__':main()
