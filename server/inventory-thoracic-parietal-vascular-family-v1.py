"""Current-source authority inventory and literal thoracic vascular frame.

Inventory and source/interface measurements precede any authored branch.
No generic target centroid is treated as a branch takeoff or tissue entry.
"""
from pathlib import Path
import json,sys,time,sqlite3
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import patch_components
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.tube_station_topology import recover_station_rings

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
CENSUS=ROOT/'reports/ta98-blood-circulation-inventory-installed-v12.json'
AUTH=ROOT/'.local-models/ta98-authority.sqlite'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
FRAME='model:local-ta98:assembled-metres'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    a=json.loads(CAT.read_bytes());ps={p['id']:p for p in a['parts']};census=json.loads(CENSUS.read_bytes());assert census['source_catalogue_sha256']==EXPECTED
    terms={q['id']:q for q in census['terms']};db=sqlite3.connect('file:'+str(AUTH)+'?mode=ro',uri=True);db.row_factory=sqlite3.Row;authority={r['id']:dict(r)for r in db.execute('select * from ta98')}
    def under(t,parent):
        seen=set()
        while t and t not in seen:
            if t==parent:return True
            seen.add(t);t=authority[t].get('parent_id')if t in authority else None
        return False
    rows=[]
    for t,q in authority.items():
        if not(under(t,'A12.2.08.029')or under(t,'A12.2.11.001')):continue
        s=terms[t];rows.append(dict(id=t,parent_id=q['parent_id'],name_en=q['name_en'],name_la=q['name_la'],type_of_entity=q['type_of_entity'],
            classification='whole-organ',representationEntity='Named arterial conduit, branching network or collective named vessel group; an outer source cast does not establish a finite wall or microscopic vessel.',
            sourcePartIds=s['resolved_direct_part_ids'],activeExactSourcePartIds=s['active_exact_source_part_ids'],representationStatus=s['representation_status'],missingRequirements=s['missing_requirements']))
    rows.sort(key=lambda r:r['id']);assert len(rows)==35 and sum(bool(r['sourcePartIds'])for r in rows)==8
    # Save all applicable terms BEFORE any source derivation or construction.
    atomic_json(OUT/'authority-inventory.json',dict(sourceCatalogueSha256=EXPECTED,authoritySha256=digest_file(AUTH),censusSha256=digest_file(CENSUS),terms=rows,constructedGeometry=False))
    report=ROOT/'reports/mesh-audit/REVIEW.md';text=report.read_text();heading='## Thoracic parietal arterial family: complete inventory before construction v1';assert heading not in text
    lines=['',heading,'',f'Current male source `{EXPECTED}`. Authoritative roots A12.2.08.029 (internal thoracic artery) and A12.2.11.001 (thoracic aorta): 35 terms, eight current direct source/network references and 27 missing direct representations. All are named macroscopic arterial conduits/branches or collective networks; none is an opening, wall layer or microscopic tissue. Existing outer casts do not establish finite walls, shared annuli or patency. The 35-term inventory is source-bound below before new geometry.','',
        'The fixed skeleton, installed heart and original coronary courses cannot move. A branch must have an explicit native parent takeoff, a recipient tissue-interface patch, a justified anatomical route and complete current physical context; a structure name, bounding-box centroid or nearest unrelated surface is insufficient. Primary sources reviewed: PMID 32844814 (100-cadaver ITA study), PMID 26488008 / PMC4583587 (50-cadaver ITA branches), PMID 32435892 (43-cadaver posterior intercostal networks), PMID 15814947 / PMC7977093 (five cadaveric aortas and segmental origins).','',
        '| TA98 | Parent | Name | Current source/network controls | Current status |','|---|---|---|---|---|']
    for r in rows:lines.append('| '+' | '.join([r['id'],r['parent_id'],r['name_en'],', '.join(r['sourcePartIds'])or'—',r['representationStatus']])+' |')
    report.write_text(text+'\n'.join(lines)+'\n')
    pids=sorted(set(pid for r in rows for pid in r['sourcePartIds']))
    pids+=['ZA:Subclavian artery.l','ZA:Subclavian artery.r','ZA:First posterior intercostal artery.l','ZA:First posterior intercostal artery.r',
        'ZA:Second posterior intercostal artery.l','ZA:Second posterior intercostal artery.r','ZA:Supreme intercostal artery.l','ZA:Supreme intercostal artery.r',
        'ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6','ZA:Transversus thoracis muscle.l','ZA:Transversus thoracis muscle.r',
        'LOCAL:male:l_phrenic_nerve:1860','LOCAL:male:r_phrenic_nerve:1885','ZA:Internal intercostal muscles.l','ZA:Internal intercostal muscles.r',
        'ZA:Innermost intercostal muscles.l','ZA:Innermost intercostal muscles.r']
    sources=[];reader=AtlasReader()
    for pid in dict.fromkeys(pids):
        if pid not in ps:continue
        p=ps[pid];m=reader.read_part(a,p,frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED);path=OUT/(identity(pid)[:20]+'.npz');m.save(path)
        physical,first,iv=np.unique(m.vertices,axis=0,return_index=True,return_inverse=True);f=iv[m.faces]
        topo=topology(m);loop=None;course=None
        if p['system']=='arterial':
            try:ll=boundary_loops(f);loop=[dict(physicalVertexIds=q,compactVertexIds=first[q],points_m=physical[q],center_m=physical[q].mean(0),actualSourceBoundary=True,finiteAnnulus=False)for q in ll]
            except ValueError as exc:loop=dict(unresolved=str(exc))
            try:
                cc=recover_station_rings(physical,f);course=[]
                for g in cc:
                    centers=np.array([physical[q].mean(0)for q in g['rings']]);course.append(dict(sourceFaceIds=g['face_ids'],stationRingPhysicalVertexIds=g['rings'],centers_m=centers,
                        unbranchedSourceCylinder=True,sourceComponent=g['source_component']))
            except ValueError as exc:course=dict(unresolved=str(exc))
        sources.append(dict(id=pid,name=p['name'],system=p['system'],record=p,recordSha256=identity(p),mesh=str(path),meshSha256=digest_file(path),
            geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'],bounds_m=m.bounds,
            topology=topo,sourceMouths=loop,sourceCourses=course,sourceIdentityOnlyNotAnatomicalRegistration=True))
    result=dict(schema='thoracic-parietal-blood-family-source-inventory/v1',sourceCatalogueSha256=EXPECTED,terms=rows,sources=sources,
        authoritySha256=digest_file(AUTH),censusSha256=digest_file(CENSUS),allAuthorityParentsPreserved=True,noAuthoredBranchYet=True,
        sourceNameIsNotTakeoffOrRecipientEvidence=True,noPhysicalOrCatalogueChange=True,
        primaryReferences=[dict(url='https://pubmed.ncbi.nlm.nih.gov/32844814/',support='100 cadavers: ITA sternal and pericardiacophrenic branches vary; anterior intercostals can arise from musculophrenic artery.'),
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4583587/',support='50 cadavers: ITA descends behind costal cartilages; terminal branching and branch counts are variable.'),
            dict(url='https://pubmed.ncbi.nlm.nih.gov/32435892/',support='43 cadavers: intercostal collateral networks concentrated near costal angle; not a universal uniform ladder.'),
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC7977093/',support='Five cadaveric aortas: actual segmental origins vary by vertebral level and circumferential aortic location.')],
        elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'source-inventory.json',json.loads(json.dumps(result,default=clean)));assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(output=str(OUT/'source-inventory.json'),terms=len(rows),sources=len(sources),elapsedSeconds=result['elapsedSeconds'],
        arterialSources=[dict(id=r['id'],faces=r['topology']['triangles'],components=len(r['topology']['components']),openEdges=r['topology']['open_edges'],sourceCourses=len(r['sourceCourses'])if isinstance(r['sourceCourses'],list)else r['sourceCourses'])for r in sources if r['system']=='arterial'])))

if __name__=='__main__':main()
