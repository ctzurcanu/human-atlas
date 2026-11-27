"""Validate non-solid references against their actual encoded host facets."""
import numpy as np

from .model import AtlasReader


def validate_catalogue(atlas, *, reader=None, frame=None):
    from .render_volumes import validate_liver_sections
    from .liver_completion import validate_liver_completion
    diaphragm_checks=validate_diaphragm_catalogue(atlas,reader=reader,frame=frame)+validate_liver_sections(atlas,reader=reader,frame=frame)+validate_liver_completion(atlas,reader=reader,frame=frame)
    marker=atlas.get('finiteUrinaryIntegration')or{}
    groups=[g for g in marker.get('termGroups',[])if g.get('pointReference')is not None]
    configurations=[g for g in marker.get('termGroups',[])if g.get('openingConfiguration')is not None]
    if not groups and not configurations:return diaphragm_checks
    reader=reader or AtlasReader();parts={p['id']:p for p in atlas['parts']};concepts={c['id']:c for c in atlas['concepts']}
    checks=diaphragm_checks
    for group in groups:
        ref=group['pointReference'];pid=ref['partId'];part=parts.get(pid);concept=concepts.get(group['conceptId'])
        if ref.get('representation')!='point-landmark'or not part or not concept or group['partIds']!=[pid]or concept['elements']!=[pid]or concept.get('ta98Term')!=group['term']:
            raise ValueError('A point reference needs its explicit semantic identity and host')
        default_model='local-ta98'if atlas.get('sex')=='male'else'local-female-ta98'
        mesh=reader.read_part(atlas,part,frame=frame or 'model:'+default_model+':assembled-metres')
        if mesh.source['bank_sha256']!=ref.get('sourceBankSha256'):
            raise ValueError('The point reference source bank changed')
        face=ref.get('face');weights=np.asarray(ref.get('barycentric'),float);point=np.asarray(ref.get('point'),float)
        if type(face)is not int or not 0<=face<len(mesh.faces)or weights.shape!=(3,)or point.shape!=(3,)or not np.isfinite(np.r_[weights,point]).all():
            raise ValueError('A point reference needs a finite actual facet and barycentric coordinates')
        bound=64*np.finfo(float).eps*max(1.,float(np.abs(mesh.triangles[face]).max()))
        if np.any(weights<0)or np.any(weights>1)or abs(weights.sum()-1)>64*np.finfo(float).eps:
            raise ValueError('Point reference weights must remain inside the source facet')
        actual=weights@mesh.triangles[face];error=float(np.linalg.norm(actual-point))
        if error>bound:raise ValueError('The point reference no longer lies on its actual encoded host facet')
        if group['term']=='A08.3.01.002':
            attachment=marker.get('apicalAttachment')or{}
            transfer=marker.get('apexTransfer')or{}
            if attachment:
                valid=pid==attachment.get('wallPartId')and ref.get('sourcePlanSha256')==attachment.get('sourcePlanSha256')and any(g.get('neighbor')==attachment.get('partId')and g['start']//3<=face<(g['start']+g['count'])//3 for g in part.get('tissueDrawGroups',[]))
            elif transfer:
                from pathlib import Path
                import json
                from .common import digest_file
                plan_path=Path(transfer['transferPlanPath'])
                if digest_file(plan_path)!=transfer['transferPlanSha256']:
                    raise ValueError('The apex transfer plan changed')
                plan=json.loads(plan_path.read_text());matrix=np.asarray(transfer['matrix'],float)
                if (matrix.shape!=(4,4)or not np.isfinite(matrix).all()or not np.array_equal(matrix,np.asarray(plan['matrix']))
                    or not np.array_equal(matrix[3],[0,0,0,1])or np.linalg.det(matrix[:3,:3])<=0
                    or not np.allclose(matrix[:3,:3].T@matrix[:3,:3],np.eye(3)*plan['uniform_scale']**2,rtol=1e-12,atol=1e-14)):
                    raise ValueError('The bladder apex must follow the established proper core similarity')
                source_record=dict(transfer['sourcePartRecord']);source_record['chunk']=0
                source=reader.read_part(dict(chunks=[transfer['sourceChunk']]),source_record,frame='model:local-ta98:assembled-metres')
                source_weights=np.asarray(transfer['sourceBarycentric'],float);source_point=np.asarray(transfer['sourcePoint'],float)
                source_face=transfer['sourceReferenceFace']
                if (source_record['id']!='BLADDER:male:body:adventitia'or transfer['sourcePartId']!=source_record['id']
                    or source.source['bank_sha256']!=transfer['sourceBankSha256']or type(source_face)is not int or not 0<=source_face<len(source.faces)
                    or source_weights.shape!=(3,)or not np.isfinite(source_weights).all()or np.any(source_weights<0)or np.any(source_weights>1)
                    or abs(source_weights.sum()-1)>64*np.finfo(float).eps
                    or np.linalg.norm(source_weights@source.triangles[source_face]-source_point)>1e-12):
                    raise ValueError('The original male apex must replay on its encoded source facet')
                mapped=matrix[:3,:3]@source_point+matrix[:3,3];distance=float(np.linalg.norm(mapped-point))
                valid=(atlas.get('sex')=='female'and transfer.get('schema')=='bladder-apex-similarity/v1'
                    and pid==transfer.get('targetPartId')=='BLADDER:female:body:adventitia'
                    and ref.get('sourcePlanSha256')==transfer.get('transferPlanSha256')
                    and 0<=transfer.get('maximumSourceOffsetM',float('inf'))<=1e-7
                    and distance<=transfer['maximumSourceOffsetM']+1e-12
                    and any(not g.get('neighbor')and g['start']//3<=face<(g['start']+g['count'])//3 for g in part.get('tissueDrawGroups',[])))
            else:valid=False
            if not valid:
                raise ValueError('The apex must use its declared wall-to-ligament attachment')
        checks.append(dict(term=group['term'],concept_id=group['conceptId'],host_part_id=pid,face=face,point_m=point.tolist(),
                           actual_source_point_m=actual.tolist(),roundtrip_error_m=error,numerical_bound_m=bound,
                           bank_path=mesh.source['bank_path'],bank_sha256=mesh.source['bank_sha256'],
                           passed=True,representation='point-landmark',new_tissue_solid=False,anatomical_acceptance=False))
    byterm={g['term']:g for g in marker.get('termGroups',[])}
    configuration_terms={'male':('A09.4.02.002','A08.3.01.027M'),
                         'female':('A09.2.03.002','A08.3.01.028F')}
    expected_configuration,expected_opening=configuration_terms.get(atlas.get('sex'),(None,None))
    for group in configurations:
        cfg=group['openingConfiguration'];parent=byterm.get(cfg.get('openingTerm'));concept=concepts.get(group['conceptId'])
        if any(not isinstance(cfg.get(k),str)or len(cfg[k])!=64 or any(c not in '0123456789abcdef' for c in cfg[k])for k in ('sourcePlanSha256','sourceBankSha256')):
            raise ValueError('Opening configuration source hashes must be explicit SHA256 identities')
        if (expected_configuration is None or group['term']!=expected_configuration or cfg.get('state')!='voiding'
            or cfg.get('openingTerm')!=expected_opening or not parent or not concept
            or concept.get('ta98Term')!=group['term']or concept.get('partOf')!=parent['conceptId']
            or group['partIds']!=parent['partIds']or concept['elements']!=group['partIds']
            or group.get('openingBoundaries')or cfg.get('sourcePlanSha256')!=marker.get('openingPlanSha256')):
            raise ValueError('An existing open funnel establishes only its explicit voiding configuration')
        boundaries=parent.get('openingBoundaries',[])
        if len(boundaries)!=1 or boundaries[0].get('name')!='neck':
            raise ValueError('The voiding configuration needs its clinical neck rim, not a histological join')
        boundary=boundaries[0]
        if boundary.get('sourcePlanSha256')!=cfg['sourcePlanSha256']or boundary.get('partIds')!=group['partIds']:
            raise ValueError('Clinical opening source ownership changed')
        n=len(boundary.get('points',[]))
        if not n or any(len(boundary.get(k,[]))!=n for k in ('samplePartIds','faceIds','barycentric')):
            raise ValueError('Clinical opening witnesses must be complete')
        cached={};errors=[]
        for pid,face,weights,point in zip(boundary['samplePartIds'],boundary['faceIds'],boundary['barycentric'],boundary['points']):
            if pid not in group['partIds']or pid not in parts:
                raise ValueError('Clinical opening host ownership changed')
            if pid not in cached:
                default_model='local-ta98'if atlas['sex']=='male'else'local-female-ta98'
                cached[pid]=reader.read_part(atlas,parts[pid],frame=frame or 'model:'+default_model+':assembled-metres')
            mesh=cached[pid]
            if mesh.source['bank_sha256']!=cfg.get('sourceBankSha256'):
                raise ValueError('Clinical opening source bank changed')
            w=np.asarray(weights,float);q=np.asarray(point,float)
            if (type(face)is not int or not 0<=face<len(mesh.faces)or w.shape!=(3,)or q.shape!=(3,)
                or not np.isfinite(np.r_[w,q]).all()or np.any(w< -1e-6)or np.any(w>1+1e-6)or abs(w.sum()-1)>1e-6
                or not any(not g.get('neighbor')and g['start']//3<=face<(g['start']+g['count'])//3 for g in parts[pid].get('tissueDrawGroups',[]))):
                raise ValueError('Clinical opening must bind actual free lining facets')
            error=float(np.linalg.norm(w@mesh.triangles[face]-q))
            bound=64*np.finfo(float).eps*max(1.,float(np.abs(mesh.triangles[face]).max()))
            if error>bound:
                raise ValueError('Clinical opening witness no longer replays on its literal source facet')
            errors.append(error)
        for pid,mesh in cached.items():
            checks.append(dict(term=group['term'],concept_id=group['conceptId'],host_part_id=pid,
                               representation='opening-configuration',configuration='voiding',witnesses=n,
                               maximum_roundtrip_error_m=max(errors),bank_path=mesh.source['bank_path'],
                               bank_sha256=mesh.source['bank_sha256'],passed=True,new_tissue_solid=False,
                               anatomical_acceptance=False))
    return checks


def validate_diaphragm_catalogue(atlas, *, reader=None, frame=None):
    """Replay opening/triangle references on the literal installed tissue.

    A semantic opening may select its actual rim tissue as context; this does
    not change the tissue's material identity or create a solid in the lumen.
    """
    marker=atlas.get('finiteDiaphragmIntegration')or{}
    if marker.get('schema')!='diaphragm-complete-local-v7':return []
    groups=marker.get('termGroups',[])
    expected={f'A04.4.02.{i:03d}'for i in range(1,18)}
    if len(groups)!=17 or{g.get('term')for g in groups}!=expected:
        raise ValueError('The complete diaphragm needs all seventeen authoritative term controls')
    is_sha=lambda s:isinstance(s,str)and len(s)==64 and all(c in '0123456789abcdef'for c in s)
    if not is_sha(marker.get('sourcePlanSha256'))or not is_sha(marker.get('bankSha256')):
        raise ValueError('Diaphragm controls need explicit source-plan and encoded-bank hashes')
    parts={p['id']:p for p in atlas['parts']};concepts={c['id']:c for c in atlas['concepts']}
    ids=marker.get('partIds',[])
    if not ids or len(ids)!=len(set(ids))or any(pid not in parts or parts[pid].get('suppressed')for pid in ids):
        raise ValueError('Diaphragm controls need unique active tissue hosts')
    reference_terms={'A04.4.02.010','A04.4.02.011','A04.4.02.014','A04.4.02.015','A04.4.02.016'}
    reader=reader or AtlasReader();meshes={};checks=[]
    for g in groups:
        term=g['term'];members=g.get('partIds',[]);concept=concepts.get(g.get('conceptId'))
        if not members or not concept or concept.get('ta98Kind')!='exact'or concept.get('ta98Term')!=term or concept['elements']!=members or any(pid not in ids for pid in members):
            raise ValueError('Diaphragm term identity or tissue membership changed')
        if term not in reference_terms:
            expected_representation='whole-organ'if term=='A04.4.02.001'else'tissue-assembly'if term=='A04.4.02.002'else'finite-tissue'
            if g.get('representation')!=expected_representation or g.get('boundaries'):
                raise ValueError('Material diaphragm controls cannot be boundary-only substitutes')
            if term=='A04.4.02.001':
                if set(members)!=set(ids):raise ValueError('Whole diaphragm membership is incomplete')
            elif term=='A04.4.02.002':
                lumbar={f'A04.4.02.{i:03d}'for i in range(2,8)}
                expected={pid for pid in ids if parts[pid].get('ta98Term')in lumbar or pid=='TA98PART:ZA:Diaphragm:Tendon:crural-attachment'}
                if set(members)!=expected:raise ValueError('The lumbar assembly must include both crura and arcuate ligaments')
            elif any(parts[pid].get('ta98Term')!=term for pid in members):
                raise ValueError('A finite tissue must retain its own anatomical term')
            continue
        boundaries=g.get('boundaries',[])
        if g.get('representation')!='boundary-reference'or not boundaries:
            raise ValueError('A diaphragm opening or triangle needs its actual boundary')
        if term in {'A04.4.02.015','A04.4.02.016'}and{b.get('name')for b in boundaries}!={'left','right'}:
            raise ValueError('Both sided diaphragm triangles must be represented')
        for boundary in boundaries:
            n=len(boundary.get('points',[]))
            if n<(3 if boundary.get('closed')else 2)or any(len(boundary.get(k,[]))!=n for k in ['samplePartIds','faceIds','barycentric']):
                raise ValueError('Diaphragm boundary witnesses are incomplete')
            if boundary.get('sourceBankSha256')!=marker['bankSha256']:
                raise ValueError('The diaphragm boundary source bank changed')
            errors=[]
            for pid,face,weights,point in zip(boundary['samplePartIds'],boundary['faceIds'],boundary['barycentric'],boundary['points']):
                if pid not in members:raise ValueError('Diaphragm boundary host ownership changed')
                if pid not in meshes:
                    meshes[pid]=reader.read_part(atlas,parts[pid],frame=frame or'model:local-ta98:assembled-metres')
                mesh=meshes[pid]
                if mesh.source['bank_sha256']!=marker['bankSha256']:
                    raise ValueError('The actual diaphragm reference bank is stale')
                w=np.asarray(weights,float);q=np.asarray(point,float)
                if(type(face)is not int or not 0<=face<len(mesh.faces)or w.shape!=(3,)or q.shape!=(3,)
                   or not np.isfinite(np.r_[w,q]).all()or np.any(w<0)or np.any(w>1)or abs(w.sum()-1)>64*np.finfo(float).eps):
                    raise ValueError('A diaphragm boundary must use finite barycentrics on an actual host facet')
                error=float(np.linalg.norm(w@mesh.triangles[face]-q))
                bound=64*np.finfo(float).eps*max(1.,float(np.abs(mesh.triangles[face]).max()))
                if error>bound:raise ValueError('A diaphragm boundary left its encoded host facet')
                errors.append(error)
            checks.append(dict(term=term,boundary=boundary.get('name'),representation='boundary-reference',
                witnesses=n,maximum_roundtrip_error_m=max(errors),bank_sha256=marker['bankSha256'],
                bank_path=mesh.source['bank_path'],
                passed=True,new_tissue_solid=False,anatomical_acceptance=False))
    return checks
