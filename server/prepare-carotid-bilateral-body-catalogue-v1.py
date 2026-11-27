"""Prepare one reviewed-bank catalogue candidate; never install live anatomy."""
from pathlib import Path
import copy,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.install import catalogue_diff

H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1'
EXPECTED='a9ae0fbaaf62fdf600a966e43641f54a4b20621944784e31ec43f43260428638'
NAME='Carotid sinus medial wall reference (left)'
DESCRIPTION='Qualified reference on the available medial outer wall of the proximal left internal carotid artery. Boundaries are provisional and coverage of the first centimetre is partial. Sinus bulb morphology, lumen extent and baroreceptor histology remain unverified.'

def main():
    began=time.monotonic();O=B/'bilateral-body-catalogue-v1';O.mkdir(exist_ok=False)
    master=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    female=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'
    assert H(master)==EXPECTED and H(female)=='8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'
    before=json.loads(master.read_bytes());old={p['id']:p for p in before['parts']}
    intermediate=B/'right-catalogue-adapter-v1/proposed.json'
    rp=json.loads((intermediate.parent/'adapter-proof.json').read_bytes())
    assert rp['passed']and rp['proposed_catalogue_sha256']==H(intermediate)and rp['source_catalogue_sha256']==EXPECTED
    after=json.loads(intermediate.read_bytes());assert len(after['parts'])==7010
    body=B/'body-at-installed-fork-v1/packed-body-v1'
    tissue=json.loads((body/'tissue.json').read_bytes());audit=json.loads((body/'independent-encoded-audit.json').read_bytes())
    assert tissue['passed']and audit['complete']and audit['passed']and H(body/'cells.bin')==audit['bank_sha256']==tissue['bank_sha256']
    assert audit['source_catalogue_sha256']==EXPECTED and audit['one_exterior_only']and not audit['added_capsule_geometry']and not audit['added_lumen_geometry']
    row=copy.deepcopy(tissue['parts']['body']);pid=row['id'];assert pid not in old
    chunk=len(after['chunks']);after['chunks'].append(dict(url='/local-models/'+str((body/'cells.bin').relative_to(ROOT/'.local-models')),bytes=(body/'cells.bin').stat().st_size))
    material='TA98BODY:carotid-body:illustrative-organ'
    assert material not in after.get('materials',{})
    after.setdefault('materials',{})[material]=dict(color=[.65,.42,.28],roughness=.65,metalness=0)
    description='A small chemoreceptor organ at the posteromedial aspect of the left carotid bifurcation. This original analytic macro-body is positioned from the actual installed fork. It is an organ, not an arterial opening or a capsule.'
    limits='Provisional macroscopic shape and size only. No distinct histological capsule, glomic arterial entry, capillary bed, carotid-sinus nerve or stromal attachment is constructed. Vessel and nerve interfaces remain unfinished.'
    row.update(conceptId=pid,name='Carotid body (left)',sourceId='Constructed left carotid chemoreceptor macro-body',
        chunk=chunk,material=material,system='arterial',tissue='organ',anatomicalEntityType='chemoreceptor-organ',
        displaySystemQualification='Arterial collection placement follows TA98 A12.2.04.006; this organ is not an artery.',
        ta98Term='A12.2.04.007',ta98Kind='exact',ta98Parent='A12.2.04.006',laterality='left',
        description=description,coverageLimitation=limits,pairedTissue=True,sectionSolid=True,
        sectionAssembly='carotid-body:male:left',license='CC BY-SA 4.0',
        provenance=dict(label='Human Atlas original analytic construction in the Z-Anatomy native frame',url='https://github.com/Z-Anatomy/Models-of-human-anatomy',detail='Original macro-body geometry; no source carotid-body mesh copied. Size and placement use cited primary anatomical evidence and the local installed arterial frame.'),
        groups=['5: Cardiovascular system','Aorta','Common carotid artery','Neck','Systemic arteries'],
        reviewStatus='qualified-local-macro-body; perfusion-and-attachment-unfinished')
    after['parts'].append(row)
    after['concepts'].append(dict(id=pid,name=row['name'],elements=[pid],partOf='ZA:Left common carotid artery',
        ta98Term=row['ta98Term'],ta98Kind='exact',description=description,coverageLimitation=limits))
    after.setdefault('ta98TermBindings',[]).append(dict(term=row['ta98Term'],conceptId=pid,partIds=[pid],kind='exact'))
    cca=next(c for c in after['concepts']if c['id']=='ZA:Left common carotid artery')
    cca['elements']=[*cca['elements'],pid]
    cca['coverageLimitation']='Stored finite common carotid trunk and bifurcation plus an adjacent provisional carotid body. Wall dimensions, regional boundaries, body attachment and microscopic supply remain unverified.'
    feature_path=B/'qualified-carotid-sinus-surface-reference-v1/feature.json';feature=json.loads(feature_path.read_bytes())
    assert feature['sourceCatalogueSha256']==EXPECTED and feature['referenceOnly']and feature['newMeshCount']==0
    g=copy.deepcopy(feature['applicationTermGroup']);cid=g['conceptId'];host=g['partIds'][0]
    assert not any(c['id']==cid for c in after['concepts'])and g['term']=='A12.2.04.008'
    assert host=='TA98PART:ZA:Internal carotid artery.l:A12.2.06.002'
    assert old[host]==next(p for p in after['parts']if p['id']==host)
    after['concepts'].append(dict(id=cid,name=NAME,elements=[host],partOf='ZA:Left common carotid artery',
        ta98Term=g['term'],ta98Kind='exact',description=DESCRIPTION,coverageLimitation=DESCRIPTION))
    after['carotidSurfaceFeatureIntegration']=dict(schema='carotid-surface-features/v1',side='left',
        hostPartId=host,bankSha256=feature['sourceBankSha256'],bankUrl=feature['bankUrl'],
        sourcePlanSha256=feature['selectionPlanSha256'],provisional=True,termGroups=[g])
    after['carotidBodyIntegration']=dict(schema='carotid-body-macro-v1',partIds=[pid],wholeParent='ZA:Left common carotid artery',
        sourceCatalogueSha256=EXPECTED,bankSha256=audit['bank_sha256'],auditSha256=H(body/'independent-encoded-audit.json'),
        anatomicalAcceptance=False,perfusionAndAttachmentAccepted=False)
    after['partsCount']=len(after['parts']);after['triangles']=sum(p.get('triangles',p['indexCount']//3)for p in after['parts'])
    expected=set(rp['actual_changed_parts'])|{pid};_,_,changed=catalogue_diff(before,after,expected)
    assert set(changed)==expected and len(after['parts'])==7011
    assert all(p==next(q for q in after['parts']if q['id']==p['id'])for p in before['parts']if p['id']not in expected)
    oldc={c['id']:c for c in before['concepts']};newc={c['id']:c for c in after['concepts']}
    changed_concepts={'ZA:Left common carotid artery','ZA:Right common carotid artery'}
    assert all(newc[id]==c for id,c in oldc.items()if id not in changed_concepts)
    reader=AtlasReader();m=reader.read_part(after,row,frame=tissue['frame'])
    assert m.geometry_sha==audit['geometry_sha256']and m.attributes_sha==audit['attributes_sha256']
    assert H(master)==EXPECTED and H(female)=='8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'
    out=O/'proposed.json';out.write_text(json.dumps(after,separators=(',',':'))+'\n')
    proof=dict(schema='bilateral-carotid-body-and-surface-reference-catalogue/v1',passed=True,
        source_catalogue_sha256=EXPECTED,proposed_catalogue_sha256=H(out),changed_parts=sorted(expected),
        unchanged_original_records=7000,new_physical_parts=2,new_surface_reference_meshes=0,
        body_geometry_sha256=m.geometry_sha,body_attributes_sha256=m.attributes_sha,
        qualified_surface_faces=len(g['surface']['faces']),source_catalogues_unchanged=True,
        anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    (O/'adapter-proof.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))

if __name__=='__main__':main()
