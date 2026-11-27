"""Independent replay of liver references and representation ownership."""
import numpy as np
from .model import AtlasReader
from .common import digest_file

def validate_liver_completion(atlas,*,reader=None,frame=None):
    m=atlas.get('liverCompletionIntegration')
    if m is None:return []
    reader=reader or AtlasReader();parts={p['id']:p for p in atlas['parts']};concepts={c['id']:c for c in atlas['concepts']}
    if atlas.get('sex')not in ['male','female']or m.get('schema')!=atlas['sex']+'-liver-completion/v1'or len(m.get('partIds',[]))!=16 or len(set(m['partIds']))!=16 or len(m.get('termGroups',[]))!=30:
        raise ValueError('Incomplete liver completion ownership')
    host=parts[m['hostPartId']];mesh=reader.read_part(atlas,host,frame=frame or'model:local-ta98:assembled-metres')
    if mesh.source['bank_sha256']!=m['hostBankSha256']:raise ValueError('Stale liver reference host')
    checks=[]
    all_groups=[*m['termGroups'],*atlas['liverFeatureIntegration']['termGroups']]
    for group in all_groups:
        c=concepts.get(group['conceptId']);ids=group['partIds']
        if not c or c.get('ta98Term')!=group['term']or c.get('ta98Kind')!='exact'or c['elements']!=ids or any(pid not in parts or parts[pid].get('suppressed')for pid in ids):
            raise ValueError('Liver representation has an invalid term or host')
        witnesses=([group['anchor']]if group.get('anchor')else[])+[w for curve in group.get('curves',[])for w in curve['witnesses']]
        errors=[]
        for w in witnesses:
            if w['partId']!=host['id']or type(w['face'])is not int or not 0<=w['face']<len(mesh.faces):raise ValueError('Liver reference left its actual host')
            weights=np.asarray(w['barycentric'],float);point=np.asarray(w['point'],float)
            if weights.shape!=(3,)or point.shape!=(3,)or not np.isfinite(np.r_[weights,point]).all()or weights.min()<0 or weights.max()>1 or abs(weights.sum()-1)>1e-12:raise ValueError('Invalid liver facet barycentrics')
            error=float(np.linalg.norm(weights@mesh.triangles[w['face']]-point));errors.append(error)
            if error>128*np.finfo(float).eps*max(1.,abs(mesh.triangles[w['face']]).max()):raise ValueError('Liver source-boundary point no longer replays')
        if 'surface'in group:
            faces=group['surface']['faces'];count=host.get('liverReferenceTriangleCount',host['tissueDrawGroups'][0]['count']//3)
            if not faces or len(set(faces))!=len(faces)or any(type(f)is not int or not 0<=f<count for f in faces):raise ValueError('Liver surface footprint left the actual outer skin')
        checks.append(dict(term=group['term'],representation=group['representation'],passed=True,witnesses=len(witnesses),maximum_replay_error_m=max(errors,default=0),bank_path=mesh.source['bank_path'],bank_sha256=mesh.source['bank_sha256'],anatomical_acceptance=False))
    sections=m.get('functionalSections',{})
    if sections.get('schema')!='functional-liver-sections/v1'or len(sections.get('exclusions',[]))!=5:raise ValueError('Functional territories need their conduit exclusions')
    path=reader.asset(sections['bankUrl']);bank=path.read_bytes()
    if digest_file(path)!=sections['bankSha256']:raise ValueError('Functional section bank changed')
    source=atlas['liverSectionIntegration'];sourcepath=reader.asset(source['bankUrl']);sourcebank=sourcepath.read_bytes()
    if sections.get('sourceBankSha256')!=digest_file(sourcepath):raise ValueError('Functional exclusion source changed')
    for a,b in zip(sections['exclusions'],source['exclusions']):
        if a['id']!=b['id']or a['geometrySha256']!=b['geometrySha256']or a['origin']!=b['origin']:raise ValueError('Functional section exclusion identity changed')
        for key,length in [('positions',a['vertexCount']*12),('indices',a['indexCount']*4)]:
            if bank[a[key]:a[key]+length]!=sourcebank[b[key]:b[key]+length]:raise ValueError('Functional section exclusions are not literal audited geometry')
    for pid in m['partIds']:
        p=parts[pid]
        if not p.get('pairedTissue')or not p.get('sectionSolid')or not p.get('defaultHidden'):raise ValueError('Optional liver structures need explicit finite geometry and visibility')
    for pid in m['functionalPartIds']:
        p=parts[pid]
        if atlas['chunks'][p['chunk']]['url']!=sections['bankUrl']:raise ValueError('Functional geometry lost its section classifier bank')
        actual=set(p.get('liverComputationalFaces',[]));declared=set()
        for g in p['tissueDrawGroups']:
            if g.get('computationalOnly'):declared.update(range(g['start']//3,(g['start']+g['count'])//3))
        if actual!=declared:raise ValueError('A computational hilar disk entered the displayed surface')
    represented={g['term']for marker in ['liverIntegration','hepaticConnectionIntegration','liverFeatureIntegration','liverCompletionIntegration']for g in atlas[marker]['termGroups']if g['term'].startswith('A05.8.01.')}
    expected=69 if atlas['liverIntegration']['schema'].endswith('/v2') else 63
    if len(represented)!=expected:raise ValueError('Liver coverage must count the actual represented authority terms: '+str(expected))
    return checks
