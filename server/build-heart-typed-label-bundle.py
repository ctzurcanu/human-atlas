"""Bind audited, dimensionally typed references to both installed TA98 hearts."""
from pathlib import Path
import argparse,collections,copy,json,runpy,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--geometry-stage',required=True);ap.add_argument('--boundary-stage',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args()
started=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';volume=base/args.geometry_stage;bounds=base/args.boundary_stage
audit=json.loads((volume/'float32-audit.json').read_text());assert audit['complete']and all(r['passed']for r in audit['geometries'].values()),'A finite buffer failed its complete topology/crossing audit'
data=json.loads((volume/'construction.json').read_text());native=json.loads((bounds/'typed-geometries.json').read_text());assert not data['remaining']and not native['issues']
typed=copy.deepcopy(native['geometries']);qualifiers={}
for term,row in data['geometries'].items():
    m=Mesh.load(volume/row['file']);typed[term]=dict(type=row['type'],vertices=m.vertices.tolist(),faces=m.faces.tolist(),**({'filled':False}if row['type']=='space-region'else {}));qualifiers[term]=row
assert len(typed)==71
types=runpy.run_path(str(ROOT/'server/heart-reference-types.py'))
assert all(g['type']==types['required_type'](term)for term,g in typed.items())
out=base/args.stage;out.mkdir(exist_ok=False)
plan=dict(schema='cardiac-typed-reference-install-plan/v1',geometries=typed,qualifiers=qualifiers,counts=dict(collections.Counter(g['type']for g in typed.values())),source_geometry_stage=args.geometry_stage,source_boundary_stage=args.boundary_stage,float32_audit_sha256=digest_file(volume/'float32-audit.json'),individual_microscopic_channels_verified=False,clinical_anatomical_acceptance=False)
atomic_json(out/'reference-plan.json',plan)
paths={sex:ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json'for sex in ['male','female']};atlases={sex:json.loads(p.read_text())for sex,p in paths.items()};male=atlases['male'];male_parts={p['id']:p for p in male['parts']};reader=AtlasReader();source_meshes={}
assert digest_file(paths['male'])==data['source_catalogue_sha256']==native['source_catalogue_sha256']
male_ep=male_parts['TA98HEART:male:epicardium'];male_R,male_scale=instance_pose(male_ep);male_origin=np.array(male_ep['geometryOrigin']);to_local=lambda q:(np.asarray(q)-male_origin)@male_R/male_scale
def canonical_anchor(term,g):
    if 'vertices'in g:
        vv=np.array(g['vertices']);tri=vv[np.array(g['faces'])].mean(1);return tri[np.argmin(np.linalg.norm(tri-vv.mean(0),axis=1))]
    if g['type']=='point-set':return np.array(g['points'][0])
    if 'paths'in g:
        pp=np.array(g['paths'][0]['points']);return pp.mean(0)if g['type']=='opening-rim'else pp[len(pp)//2]
    if 'surfaces'in g:
        s=max(g['surfaces'],key=lambda s:s['areaMetres2']);pid=s['partId']
        if pid not in source_meshes:source_meshes[pid]=reader.read_part(male,male_parts[pid],frame='assembled-metres',catalogue_path=paths['male'])
        tri=source_meshes[pid].triangles[s['faces']].mean(1);return to_local(tri[np.argmin(np.linalg.norm(tri-tri.mean(0),axis=1))])
    return None
anchors={t:canonical_anchor(t,g)for t,g in typed.items()}
for sex,a in atlases.items():
    target=out/sex;target.mkdir();bundle=copy.deepcopy(json.loads(paths[sex].with_name('cardiac-label-sites.json').read_text()));registry=bundle['integration'];rows={g['term']:g for g in registry['termGroups']};concepts={c['id']:c for c in bundle['concepts']};pmap={p['id']:p for p in a['parts']};heart=a['finiteCardiacIntegration'];ep=pmap[f'TA98HEART:{sex}:epicardium'];rot,scale=instance_pose(ep);origin=np.array(ep['geometryOrigin']);forward=lambda q:(np.asarray(q)@rot.T)*scale+origin
    def map_part(pid):
        candidate=pid.replace(':male:',f':{sex}:')
        if candidate in pmap:return candidate
        source_group=min((g for g in male['finiteCardiacIntegration']['termGroups']if pid in g['partIds']),key=lambda g:len(g['partIds']))
        return next(g for g in heart['termGroups']if g['key']==source_group['key'])['partIds'][0]
    for term,row in rows.items():
        g=copy.deepcopy(typed[term]);native_ids=[]
        if 'surfaces'in g:
            for surface in g['surfaces']:
                source_area=surface['areaMetres2'];surface['partId']=map_part(surface['partId']);surface['areaMetres2']=source_area*(scale/male_scale)**2;native_ids.append(surface['partId'])
        if 'partIds'in g:g['partIds']=[map_part(pid)for pid in g['partIds']];native_ids.extend(g['partIds'])
        if 'sourcePartIds'in g:native_ids.extend(map_part(pid)for pid in g['sourcePartIds'])
        carrier=map_part(qualifiers.get(term,{}).get('carrierPartId',rows[term]['anchor']['partId'].replace(f':{sex}:',':male:')))
        if 'surfaces'in g:carrier=max(g['surfaces'],key=lambda s:s['areaMetres2'])['partId']
        row['partIds']=list(dict.fromkeys([carrier,*([]if 'surfaces'in g else row['partIds']),*native_ids]));row['anchor']['partId']=carrier
        if anchors[term]is not None:row['anchor']['point']=forward(anchors[term]).tolist()
        row['geometry']=g;row['representation']=g['type'];row['points']=None
        c=concepts[row['conceptId']];c['elements']=row['partIds'];c['referenceGeometryType']=g['type'];kind=g['type'].replace('-',' ')
        detail=qualifiers.get(term,{}).get('method','Native source-bound '+kind+'. Anatomical subdivision boundaries remain qualified.')
        detail=detail.replace('source-sized','source-guided')
        if sex=='female'and term in ['A12.1.08.003','A12.1.08.004']:
            detail='Cardiac-similarity transfer of the male attachment reference; female sternal, carinal and diaphragmatic insertions remain provisional.'
        if g['type']=='point-set':detail='Three source-bound leaflet commissural attachment junctions, represented as points.'
        if g['type']=='opening-rim':detail='Unfilled anatomical mouth-rim reference; no disk or solid plug is added.'
        if term=='A12.1.01.004':detail='Representative endocardial foraminal rings; individual microscopic venous channels and their exact specimen locations are not established.'
        if term=='A12.1.01.007':detail='Closed adult septal surface site; no patent fetal foramen is created.'
        if term=='A12.1.01.009':detail+=' Native coronary-sinus end ring; its atrial tissue join remains unverified.'
        if g['type']=='space-region':detail+=' The displayed wireframe bounds empty space, not tissue.'
        c['coverageLimitation']='Provisional '+kind+' reference. '+detail
    # A collective owns its actual child carriers, including pre-existing core
    # tissues. This keeps selecting a collective useful in isolation mode.
    core={g['term']:g for g in heart['termGroups']if g['kind']=='exact'}
    def carriers(term,seen=None):
        seen=set()if seen is None else seen
        if term in seen:return []
        seen.add(term);g=rows.get(term)or core.get(term)
        if not g:return []
        geometry=g.get('geometry',{})
        if geometry.get('type')=='collection':return list(dict.fromkeys(pid for child in geometry['terms']for pid in carriers(child,seen)))
        return g['partIds']
    for term,g in rows.items():
        if g['geometry']['type']=='collection':g['partIds']=list(dict.fromkeys([g['anchor']['partId'],*carriers(term)]));concepts[g['conceptId']]['elements']=g['partIds'];concepts[g['conceptId']]['coverageLimitation']='Provisional collection of typed child references and existing tissue meshes; individual limitations remain attached to each child.'
    poses={pid:dict(geometryOrigin=pmap[pid].get('geometryOrigin',[0,0,0]),geometryInstance=pmap[pid].get('geometryInstance'),chunk=pmap[pid]['chunk'],positions=pmap[pid]['positions'],indices=pmap[pid]['indices'],indexCount=pmap[pid]['indexCount'])for g in rows.values()for pid in g['partIds']}
    registry.update(carrierPoses=poses,referencePlanSha256=digest_file(out/'reference-plan.json'),typedGeometryCount=71,typedGeometryCounts=plan['counts'],typedReferenceRevision=args.stage)
    atomic_json(target/'cardiac-label-sites.json',bundle)
    atomic_json(target/'label-proof.json',dict(source_catalogue_sha256=digest_file(paths[sex]),all_parts_literal=True,all_chunks_literal=True,fixed_skeleton=True,qualified_reference_terms=71,geometry_backed_terms=42,labeled_terms=113,all_anchors_finite=all(np.isfinite(g['anchor']['point']).all()for g in rows.values()),new_atlas_part_meshes=0,new_reference_meshes=36,new_meshes=36,new_textures=0,typed_geometry_count=71,float32_audit_complete=True,float32_audit_passed=True,float32_audit_sha256=digest_file(volume/'float32-audit.json'),anatomical_acceptance=False))
    print(sex,'71 typed reference controls; 113 total; body and skeleton unchanged',flush=True)
print('Built in',round(time.monotonic()-started,2),'seconds',flush=True)
