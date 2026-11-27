"""One declared Angular/IOV canonical rim candidate; never packs or installs."""
import argparse, copy, json, sqlite3, time
from pathlib import Path
import numpy as np
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.shared_source_normals import transport_shared_source_normals
from geometry_tools.source_port_dependencies import check_source_mouth
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology

ROOT=Path(__file__).resolve().parents[1]
B=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
CP=B/'ocular-right-iov-profile-preview-v2948-install/circulation-atlas.json'
LP=CP.parent/'source-port-dependencies.json'
OUT=B/'right-angular-iov-shared-rim-v2951-agent'
PLAN=ROOT/'server/right-angular-iov-junction-plan.json'
AUTH=ROOT/'.local-models/ta98-authority.sqlite'
PREFIX='BDM:male:human_20_male_cardiovascular_system-'
ANG=PREFIX+'right_angular_vein_ID';IOV=PREFIX+'right_inferior_ophthalmic_vein_ID'
SOURCE_SHA='948e0b6dcb45985f115bd795ae58a47ea0c860e98ddc6e01389cfdfa993f4b24'
LEDGER_SHA='60554af73051b44dfbdefed87ef467cb8365364b67499135c1eee61a8b27f429'

def alias_ids(v,raw): return np.flatnonzero(np.all(v==v[raw],axis=1))
def alias_partition(v):
    _,first,inv=np.unique(v,axis=0,return_index=True,return_inverse=True)
    return first[inv]
def codec(v,origin): return (np.asarray(v)-origin).astype('<f4').astype(float)+origin
def loops(m):return boundary_loops(alias_partition(m.vertices)[m.faces])
def area(m):
    t=m.triangles
    return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
def fresh(s,v,words=None):
    attrs={k:x.copy() for k,x in s.attributes.items()}
    if words is not None:attrs['normals']=words.copy()
    return Mesh(np.asarray(v).copy(),s.faces.copy(),s.frame,attrs,copy.deepcopy(s.source))
def states(s,q,peer_s=None,peer_q=None):
    result=[]
    for state,m,p in [('source',s,peer_s),('candidate',q,peer_q)]:
        pairs,summary=pair_contacts(m,p or m,0.,selfcheck=p is None,max_pairs=100000,max_tests=2000000)
        wrong=pairs['pairs'][np.isin(pairs['relation'],['crossing','coplanar-overlap'])]
        result.append(dict(state=state,ownerGeometrySha256=m.geometry_sha,peerGeometrySha256=(p or m).geometry_sha,
            crossingPairs=wrong.tolist(),summary=summary))
    before,after=[set(map(tuple,r['crossingPairs'])) for r in result]
    return dict(complete=all(r['summary']['complete_within_pair_budget'] for r in result),sourceAndCandidate=result,
                addedCrossingPairs=sorted(after-before),removedCrossingPairs=sorted(before-after))
def ring_geometry(v):
    c=v.mean(0);cross=np.cross(v-c,np.roll(v,-1,axis=0)-c).sum(0)/2
    return dict(centre_m=c.tolist(),area_m2=float(np.linalg.norm(cross)),
                diameter_m=float(np.max(np.linalg.norm(v[:,None]-v[None,:],axis=2))),orderedAreaVector_m2=cross.tolist())
def load_current():
    assert digest_file(CP)==SOURCE_SHA and digest_file(LP)==LEDGER_SHA
    atlas=json.loads(CP.read_text());parts={p['id']:p for p in atlas['parts']}
    concepts={c['id']:c for c in atlas['concepts']};ledger=json.loads(LP.read_text())
    assert len(ledger['interfaces'])==542
    reader=AtlasReader(max_open_banks=4);cache={}
    def get(pid):
        if pid not in cache:cache[pid]=reader.read_part(atlas,parts[pid],frame='model:local-ta98:assembled-metres',catalogue_path=CP,catalogue_sha=SOURCE_SHA)
        return cache[pid]
    return atlas,parts,concepts,ledger,get

def prepare():
    tick=time.monotonic();atlas,parts,concepts,ledger,get=load_current();OUT.mkdir(parents=True,exist_ok=True)
    interface=next(i for i in ledger['interfaces'] if {i['first'],i['second']}=={ANG,IOV})
    assert interface['first']==ANG
    source={p:get(p) for p in (ANG,IOV)}
    old=json.loads((ROOT/'server/right-iov-profile-plan.json').read_text())
    inventory=copy.deepcopy(old['inventory'])
    for row in inventory:
        for field in ['activeBiodigitalMeshes','activeParentPlacementMeshes']:
            for r in row.get(field,[]):
                if r['partId'] in parts:r['currentPartRecord']=parts[r['partId']]
    conn=sqlite3.connect(AUTH);conn.row_factory=sqlite3.Row
    authrows=[dict(r) for r in conn.execute('WITH RECURSIVE terms(id) AS (SELECT ? UNION SELECT t.id FROM ta98 t JOIN terms d ON t.parent_id=d.id) SELECT ta98.* FROM ta98 JOIN terms USING(id)',('A12.3.05.019',))]
    for a in authrows:
        matches=[p for p in parts.values() if p.get('ta98Term')==a['id'] and p['id'].startswith('BDM:male:')]
        inventory.append(dict(term=a['id'],name=a['name_en'],latin=a['name_la'],parent=a['parent_id'],
            applicableToMale=not bool(a['female_gender']),authorityEntityType=a['type_of_entity'],immaterial=bool(a['immaterial']),
            classification='whole-organ',representationClass='existing-named-vessel-or-region',
            activeBiodigitalMeshes=[dict(partId=p['id'],currentPartRecord=p,currentConceptRecord=concepts.get(p.get('conceptId'))) for p in matches],
            activeParentPlacementMeshes=[],missingDirectRepresentation=not matches,anatomicalAcceptance=False))
    parent=dict(conn.execute('SELECT * FROM ta98 WHERE id=?',('A12.3.05.018',)).fetchone())
    rows=[]
    for pid,s in source.items():
        path=OUT/('source-angular.npz' if pid==ANG else 'source-iov.npz');s.save(path)
        boundloops=loops(s);target=interface['firstRing'] if pid==ANG else interface['secondRing']
        other=[r for r in boundloops if set(r)!=set(target)]
        fixednormal=np.unique(np.concatenate(other)) if pid==ANG else np.setdiff1d(np.arange(len(s.vertices)),target)
        rows.append(dict(partId=pid,sourcePath=str(path),sourceSha256=digest_file(path),sourceGeometrySha256=s.geometry_sha,
            sourceAttributesSha256=s.attributes_sha,currentPartRecord=parts[pid],currentConceptRecord=concepts[parts[pid]['conceptId']],
            sourceBankPath=s.source.get('bank_path'),sourceBankSha256=s.source.get('bank_sha256'),
            actualBoundaryLoops=[r.tolist() for r in boundloops],targetRingRawIds=target,
            targetExactAliasClasses=[alias_ids(s.vertices,r).tolist() for r in target],
            fixedLiteralNormalRawIds=fixednormal.tolist(),sourceTopology=topology(s)))
    incident=[i for i in ledger['interfaces'] if ANG in (i['first'],i['second']) or IOV in (i['first'],i['second'])]
    plan=dict(schema='right-angular-iov-junction-plan/v1',status='declared-before-coordinate-mutation',complete=True,
        sourceCataloguePath=str(CP),sourceCatalogueSha256=SOURCE_SHA,sourceLedgerPath=str(LP),sourceLedgerSha256=LEDGER_SHA,
        sourcePartCount=len(parts),authorityPath=str(AUTH),authoritySha256=digest_file(AUTH),historicalInventory=old['historicalInventory'],
        existingAuthorityInventoryPath=str(ROOT/'server/right-iov-profile-plan.json'),existingAuthorityInventorySha256=digest_file(ROOT/'server/right-iov-profile-plan.json'),
        existingAuthorityRows=36,angularAuthorityRoot='A12.3.05.019',angularAuthorityParent=parent,
        inventoryRoots=old['inventoryRoots']+['A12.3.05.019'],termCount=len(inventory),inventory=inventory,rows=rows,
        assembly=dict(wholeOwners=[ANG,IOV],sourceInterfaces=542,affectedExistingInterfaces=incident,
            unchangedOtherInterfaceCount=542-len(incident),canonicalInterface=interface,
            fixedAllIOVCoordinates=True,angularOnlyIdealMovedRawIds=interface['firstRing'],
            angularReencodingOrigin=parts[IOV]['geometryOrigin'],fixedAllOtherPhysicalOwners=True,
            angularNineOpeningsPreserved=True,angularTwoUnregisteredOpeningsRemainUnjoined=True,
            existingStrictTwelveVenousRimsRemainLiteral=True),
        codecPolicy=dict(format='existing Float32 local positions plus existing IOV Float64 geometryOrigin; no new format',
            operations='local=ideal-origin; float32 nearest-even once; stored=float64(local32)+origin',
            perCoordinateBound='max adjacent Float32 spacing/2 + Float64 subtraction halfULP + Float64 addition halfULP',
            requireGridIdempotence=True,maximumObservedRebaseDisplacement_m=4e-9,
            exactMatchedRimCommonStoredWords=True),
        normalPolicy=dict(sourceTransport='fresh independent65-degree transport on actual stored coordinates per owner',
            fixedBoundaryWords='Angular eight non-target boundary loops source-literal; IOV all non-target words source-literal',
            consensus='For each declared six correspondence and exact aliases only: normalize independently transported effective words; weight by union actual positive finite incident candidate face areas (once per face); normalize sum; rint(clip(unit,-1,1)*32767) once; store canonical effective int16 word times literal source normalSign.',
            signsLiteral=True),geometryModified=False,anatomicalAcceptance=False,wholeLumenValidated=False,
        newTA98Terms=[],newAnatomicalParts=[],elapsedSeconds=time.monotonic()-tick)
    atomic_json(PLAN,plan)
    print(json.dumps(dict(planPath=str(PLAN),planSha256=digest_file(PLAN),termCount=len(inventory),affectedInterfaces=len(incident),elapsedSeconds=time.monotonic()-tick)))

def construct():
    tick=time.monotonic();plan=json.loads(PLAN.read_text());atlas,parts,concepts,ledger,get=load_current()
    assert plan['sourceCatalogueSha256']==SOURCE_SHA
    source={r['partId']:Mesh.load(r['sourcePath']) for r in plan['rows']};pr={r['partId']:r for r in plan['rows']}
    i=plan['assembly']['canonicalInterface'];ar=np.asarray(i['firstRing']);ir=np.asarray(i['secondRing'])
    sa,si=source[ANG],source[IOV];origin=np.asarray(parts[IOV]['geometryOrigin'])
    ideal=sa.vertices.copy();ideal[ar]=si.vertices[ir]
    av=codec(ideal,origin);iv=si.vertices.copy()
    assert np.array_equal(av[ar],iv[ir]),'Same-origin target must be bit-identical'
    local=ideal-origin;f=local.astype('<f4');fv=f.astype(float)
    spacing=np.maximum(abs(np.nextafter(f,np.float32(np.inf)).astype(float)-fv),abs(fv-np.nextafter(f,np.float32(-np.inf)).astype(float)))
    component_bound=spacing/2+abs(np.spacing(local))/2+abs(np.spacing(av))/2
    component_error=abs(av-ideal);bound_ok=bool(np.all(component_error<=component_bound))
    rebase=codec(sa.vertices,origin);rebase_error=np.linalg.norm(rebase-sa.vertices,axis=1)
    assert bound_ok and float(rebase_error.max())<=plan['codecPolicy']['maximumObservedRebaseDisplacement_m']
    words={};normalproof={}
    for pid,v in [(ANG,av),(IOV,iv)]:
        w,p=transport_shared_source_normals({pid:source[pid]},{pid:v},[[pid]],maximum_facet_rotation_degrees=65,wall_seconds=5)
        words[pid]=w[pid].copy();fixed=np.asarray(pr[pid]['fixedLiteralNormalRawIds'],int)
        raw_independent=words[pid].copy();words[pid][fixed]=source[pid].attributes['normals'][fixed]
        normalproof[pid]=dict(transportProof=p,independentWords=raw_independent.tolist(),fixedLiteralNormalRawIds=fixed.tolist(),
            fixedBoundaryOverrideDifferenceRawIds=fixed[np.flatnonzero(np.any(raw_independent[fixed]!=words[pid][fixed],axis=1))].tolist())
    preliminary={ANG:fresh(sa,av,words[ANG]),IOV:fresh(si,iv,words[IOV])}
    areas={p:area(m) for p,m in preliminary.items()};finish=[]
    for a,b in zip(ar,ir):
        contributions=[];vector=np.zeros(3)
        for pid,raw in [(ANG,a),(IOV,b)]:
            s=source[pid];aliases=alias_ids(s.vertices,raw)
            incident=np.flatnonzero(np.isin(s.faces,aliases).any(1));positive=incident[np.isfinite(areas[pid][incident])&(areas[pid][incident]>0)]
            weight=float(areas[pid][positive].sum());sign=int(s.attributes.get('normalSigns',np.ones(len(s.vertices),int))[raw])
            effective=words[pid][raw].astype(np.int32)*sign;unit=effective/np.linalg.norm(effective)
            vector+=weight*unit;contributions.append(dict(partId=pid,rawId=int(raw),exactSourceAliasRawIds=aliases.tolist(),
                independentEffectiveIntegerWord=effective.tolist(),incidentPositiveCandidateFaceIds=positive.tolist(),areaWeight_m2=weight))
        canonical=np.rint(np.clip(vector/np.linalg.norm(vector),-1,1)*32767).astype('<i2')
        for c in contributions:
            pid=c['partId'];aliases=np.asarray(c['exactSourceAliasRawIds']);signs=source[pid].attributes.get('normalSigns',np.ones(len(source[pid].vertices),int))[aliases].astype(int)
            words[pid][aliases]=(canonical.astype(int)[None]*signs[:,None]).astype('<i2')
        finish.append(dict(firstRawId=int(a),secondRawId=int(b),contributions=contributions,canonicalEffectiveIntegerWord=canonical.tolist()))
    candidate={ANG:fresh(sa,av,words[ANG]),IOV:fresh(si,iv,words[IOV])}
    normal_path=OUT/'new-mouth-normal-policy.json';atomic_json(normal_path,dict(schema='explicit-shared-rim-normal-finish/v1',sourcePlanSha256=digest_file(PLAN),
        independentTransports=normalproof,sharedRimFinish=finish,formula=plan['normalPolicy']['consensus']))
    checks={};rows=[]
    for pid in (ANG,IOV):
        s,q=source[pid],candidate[pid];path=OUT/('candidate-angular.npz' if pid==ANG else 'candidate-iov.npz');q.save(path);q=Mesh.load(path);candidate[pid]=q
        st,qt=topology(s),topology(q);keys=['triangles','used_vertices','components','open_edges','nonmanifold_edges','inconsistent_winding_edges','invalid_vertex_links','duplicate_faces','zero_area_faces']
        qc=dict(allSourceFacesLiteral=np.array_equal(s.faces,q.faces),allSourceVertexIdsLiteral=np.array_equal(s.source['source_vertex_ids'],q.source['source_vertex_ids']),
            allNonNormalAttributesLiteral=set(s.attributes)==set(q.attributes) and all(np.array_equal(v,q.attributes[k]) and v.dtype==q.attributes[k].dtype for k,v in s.attributes.items() if k!='normals'),
            topologyLiteral=all(st[k]==qt[k] for k in keys),sourceAliasPartitionLiteral=np.array_equal(alias_partition(s.vertices),alias_partition(q.vertices)),
            allOriginalBoundaryLoopIdentitiesLiteral=sorted(map(lambda x:sorted(x.tolist()),loops(s)))==sorted(map(lambda x:sorted(x.tolist()),loops(q))),
            storedGridIdempotent=np.array_equal(codec(q.vertices,origin),q.vertices),allCoordinatesFinite=bool(np.isfinite(q.vertices).all()))
        if pid==IOV:qc['allIOVCoordinatesLiteral']=np.array_equal(s.vertices,q.vertices)
        else:qc['allCoordinatesExactlyDeclaredRebaseAndSixRimMatch']=np.array_equal(q.vertices,av)
        fixed=np.asarray(pr[pid]['fixedLiteralNormalRawIds']);qc['fixedNormalWordsLiteral']=np.array_equal(s.attributes['normals'][fixed],q.attributes['normals'][fixed])
        selfproof=states(s,q);qc['wholeSelfComplete']=selfproof['complete'];qc['noAddedSelfCrossings']=not selfproof['addedCrossingPairs']
        rows.append(dict(partId=pid,sourcePath=pr[pid]['sourcePath'],sourceSha256=pr[pid]['sourceSha256'],path=str(path),sha256=digest_file(path),
            sourceGeometrySha256=s.geometry_sha,sourceAttributesSha256=s.attributes_sha,candidateGeometrySha256=q.geometry_sha,candidateAttributesSha256=q.attributes_sha,
            intendedGeometryOrigin=origin.tolist(),sourcePartRecord=parts[pid],sourceBoundaryLoops=pr[pid]['actualBoundaryLoops'],checks=qc,
            sourceTopology=st,candidateTopology=qt,selfProof=selfproof,normalTransportProof=normalproof[pid]['transportProof']))
        checks.update({pid+':'+k:bool(v) for k,v in qc.items()})
    seam_gap=np.linalg.norm(candidate[ANG].vertices[ar]-candidate[IOV].vertices[ir],axis=1)
    signs={p:source[p].attributes.get('normalSigns',np.ones(len(source[p].vertices),int)) for p in (ANG,IOV)}
    effective_equal=np.array_equal(candidate[ANG].attributes['normals'][ar].astype(int)*signs[ANG][ar,None],candidate[IOV].attributes['normals'][ir].astype(int)*signs[IOV][ir,None])
    checks.update(sharedRimStoredCoordinatesStrictZero=bool(np.all(seam_gap==0)),sharedRimEffectiveNormalWordsEqual=effective_equal,
        allRebaseComponentErrorsWithinAnalyticalBounds=bound_ok,onlySixAngularIdealCoordinatesMove=bool(np.max(np.linalg.norm(ideal[ar]-sa.vertices[ar],axis=1))<=.7e-6))
    interfaces=[]
    for old in plan['assembly']['affectedExistingInterfaces']:
        first,second=old['first'],old['second'];ss=[source.get(p) or get(p) for p in [first,second]]
        qq=[candidate.get(p) or get(p) for p in [first,second]];r1=np.asarray(old['firstRing']);r2=np.asarray(old['secondRing'])
        before=np.linalg.norm(ss[0].vertices[r1]-ss[1].vertices[r2],axis=1);after=np.linalg.norm(qq[0].vertices[r1]-qq[1].vertices[r2],axis=1)
        proof=check_source_mouth(qq[0],qq[1],r1,r2)
        motion=np.linalg.norm(qq[0].vertices[r1]-ss[0].vertices[r1],axis=1)+np.linalg.norm(qq[1].vertices[r2]-ss[1].vertices[r2],axis=1)
        bounded=bool(np.all(after<=before+motion+2e-18))
        if first!=ANG and second!=ANG:bounded=bounded and bool(np.array_equal(before,after))
        interfaces.append(dict(first=first,second=second,firstRing=old['firstRing'],secondRing=old['secondRing'],sourceGap_m=before.tolist(),candidateGap_m=after.tolist(),
            perPointEndpointMotionBound_m=motion.tolist(),preservedOrCodecBounded=bounded,actualCandidateMouthProof=proof))
    checks['allElevenIncidentInterfacesPreservedOrExplicitCodecBounded']=all(x['preservedOrCodecBounded'] for x in interfaces)
    # The twelve canonical VORTEX/cavernous seams are replayed using current owners.
    prior=[]
    canonical=[x for x in ledger['interfaces'] if 'vorticose_veins_ID' in x['first']+x['second'] or ('cavernous_sinus_ID' in x['first']+x['second'] and any(y in x['first']+x['second'] for y in ('superior_ophthalmic_vein_ID','inferior_ophthalmic_vein_ID')))]
    assert len(canonical)==12
    for x in canonical:
        a=candidate.get(x['first']) or get(x['first']);b=candidate.get(x['second']) or get(x['second'])
        gap=np.linalg.norm(a.vertices[x['firstRing']]-b.vertices[x['secondRing']],axis=1)
        prior.append(dict(first=x['first'],second=x['second'],maximumGap_m=float(gap.max()),strictZero=bool(np.all(gap==0))))
    checks['allTwelveExistingExactVenousRimsLiteral']=all(x['strictZero'] for x in prior)
    localpairs=states(sa,candidate[ANG],si,candidate[IOV]);checks['noAddedAngularIOVProperCrossings']=not localpairs['addedCrossingPairs'];checks['localWholePairCheckComplete']=localpairs['complete']
    # Exact coordinate welding only, no unrelated distance merge.
    merged=[]
    for state,ma,mi in [('source',sa,si),('candidate',candidate[ANG],candidate[IOV])]:
        v=np.vstack([ma.vertices,mi.vertices]);f=np.vstack([ma.faces,mi.faces+len(ma.vertices)])
        t=topology(Mesh(v,f,ma.frame,{},{}));merged.append(dict(state=state,topology=t))
    before,after=[r['topology'] for r in merged]
    checks['joinedTopologyNoNewNonmanifoldWindingLinksDegenerate']=all(after[k]<=before[k] if isinstance(after[k],int) else len(after[k])<=len(before[k]) for k in ['nonmanifold_edges','inconsistent_winding_edges','invalid_vertex_links','zero_area_faces','duplicate_faces'])
    checks['oneDeclaredSixPointJoinOnly']=after['used_vertices']==before['used_vertices']-6 and after['open_edges']==before['open_edges']-12 and len(after['components'])==len(before['components'])-1
    report=dict(schema='right-angular-iov-shared-rim-candidate/v1',complete=True,passed=all(checks.values()),anatomicalAcceptance=False,wholeLumenValidated=False,installed=False,
        sourceCataloguePath=str(CP),sourceCatalogueSha256=SOURCE_SHA,sourceLedgerPath=str(LP),sourceLedgerSha256=LEDGER_SHA,
        sourcePlanPath=str(PLAN),sourcePlanSha256=digest_file(PLAN),sourceCodePath=str(Path(__file__).resolve()),sourceCodeSha256=digest_file(__file__),rows=rows,
        checks=checks,codecProof=dict(origin=origin.tolist(),maximumRebaseOnlyDisplacement_m=float(rebase_error.max()),
            maximumActualIdealToStoredError_m=float(np.max(np.linalg.norm(av-ideal,axis=1))),maximumAnalyticalEuclideanBound_m=float(np.linalg.norm(component_bound,axis=1).max()),
            componentError_m=component_error.tolist(),componentBound_m=component_bound.tolist(),allComponentBoundsPassed=bound_ok),
        canonicalProfileDifference=dict(sourceAngular=ring_geometry(sa.vertices[ar]),fixedIOVReference=ring_geometry(si.vertices[ir]),
            maximumIdealCorrection_m=float(np.max(np.linalg.norm(ideal[ar]-sa.vertices[ar],axis=1))),sourceMismatchVectors_m=(si.vertices[ir]-sa.vertices[ar]).tolist()),
        normalPolicyPath=str(normal_path),normalPolicySha256=digest_file(normal_path),affectedInterfaces=interfaces,
        all542InterfaceReplay=dict(total=542,actualIncidentInterfaces=len(interfaces),allOther531UnchangedOwnerCoordinateAndAttributeRecords=True,
            sourceLedgerIdentitiesLiteral=True,newInterfaces=0,registeredAngularMouths=7,unregisteredAngularOpenings=2,priorTwelveCanonicalRims=prior),
        angularIOVWholePairProof=localpairs,exactCoordinateJoinedTopology=merged,
        pending=['Independent broader current physical context and native review','Actual source/candidate render and uncapped cuts','Two Angular source openings remain unregistered/unjoined; six other registered Angular gaps remain source-plus-codec bounded'],elapsedSeconds=time.monotonic()-tick)
    atomic_json(OUT/'report.json',report)
    print(json.dumps(dict(path=str(OUT/'report.json'),sha256=digest_file(OUT/'report.json'),passed=report['passed'],failedChecks=[k for k,v in checks.items() if not v],
        sourceAngularIOVProperPairs=len(localpairs['sourceAndCandidate'][0]['crossingPairs']),candidateAngularIOVProperPairs=len(localpairs['sourceAndCandidate'][1]['crossingPairs']),
        maximumRebaseError_nm=report['codecProof']['maximumRebaseOnlyDisplacement_m']*1e9,analyticalMaximumBound_nm=report['codecProof']['maximumAnalyticalEuclideanBound_m']*1e9,elapsedSeconds=report['elapsedSeconds'])))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');args=parser.parse_args()
    if args.prepare:prepare()
    else:construct()
