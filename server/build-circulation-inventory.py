"""Authoritative non-cardiac blood-circulation census; inventory, not acceptance.

Requires the saved actual frontend resolver export. Does not construct geometry,
infer patent lumens from labels, change authority parents or edit the atlas.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,re,sqlite3,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.ta98_authority import corrected_terms,term_sex

AUTH=ROOT/'.local-models/ta98-authority.sqlite'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
RESOLVED=ROOT/'.local-models/geometry-workbench/circulation-census-v1/resolved-male.json'
REPORT=ROOT/'reports/ta98-circulation-inventory.json'
DOC=ROOT/'TA98-CIRCULATION.md'
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def validate_effective_visibility(data,snapshot):
    """Require an ordinary-view snapshot, not raw hidden-donor availability."""
    expected=['withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)']
    view=snapshot.get('view',{})
    if (snapshot.get('schema')!='circulation-effective-frontend-source-census/v2'
        or view.get('mode')!='ordinary-local-atlas' or view.get('version')!=1
        or view.get('inspectAlternatives') is not False or view.get('pipeline')!=expected):
        raise ValueError('Ordinary-view census requires resolver snapshot v2 with Reviewed -> Preferred -> Hidden(false); regenerate the old snapshot.')
    parts={p['id']:p for p in data['parts']}
    if len(parts)!=len(data['parts']):raise ValueError('Duplicate source part IDs')
    rows=snapshot.get('effectivePartVisibility')
    if not isinstance(rows,list):raise ValueError('Missing complete effectivePartVisibility')
    visibility={r.get('id'):r for r in rows if isinstance(r,dict)}
    if len(visibility)!=len(rows) or set(visibility)!=set(parts):
        raise ValueError('effectivePartVisibility must contain every raw source ID exactly once')
    hidden=set(data.get('defaultHidden',[]))
    for pid,row in visibility.items():
        for flag in ['available','suppressed','defaultHidden','partDefaultHidden','rawSuppressed','rawDefaultHidden']:
            if type(row.get(flag))is not bool:raise ValueError('Visibility flags must be booleans: '+pid+'/'+flag)
        if row['available']==row['suppressed']:raise ValueError('Availability must be the inverse of effective suppression: '+pid)
        source=parts[pid]
        if row['rawSuppressed']!=bool(source.get('suppressed')) or row['rawDefaultHidden']!=(pid in hidden or bool(source.get('defaultHidden'))):
            raise ValueError('Visibility source flags do not match the frozen catalogue: '+pid)
        if row['rawSuppressed'] and not row['suppressed']:raise ValueError('An originally suppressed source became available: '+pid)
    for node in snapshot.get('resolved',[]):
        for key in ['directPartIds','descendantPartIds']:
            for pid in node.get(key,[]):
                if pid not in visibility or not visibility[pid]['available']:
                    raise ValueError('Resolved hierarchy contains unavailable source '+pid+' at '+str(node.get('id')))
    for key in ['rawConceptBindings','reviewedConceptBindings','conceptBindings','sourceTopLevelBindings','topLevelBindings']:
        if not isinstance(snapshot.get(key),list):raise ValueError('Snapshot is missing source/effective bindings: '+key)
    return visibility


def normal(s):
    s=re.sub(r'\s*\((?:left|right)\)\s*$','',s,flags=re.I).lower().strip()
    return re.sub(r'\s+',' ',s).strip('(). ')


def classification(row,children):
    name=row['name_en'].lower();code=row['id']
    if code in {'A12.0.00.000','A12.2.00.001','A12.3.00.001','A12.4.01.001','A13.0.00.000','A13.1.00.001','A13.2.00.001','A12.4.02.001'}:
        return 'whole-system-or-anatomical-container'
    if name in {'blood','lymph'}:return 'fluid-compartment'
    if any(x in name for x in ('tunica ','internal elastic','external elastic')):return 'wall-layer'
    if any(x in name for x in ('opening','orifice','ostium','hiatus','aperture')):return 'opening-not-a-solid'
    if row.get('immaterial'):return 'space-or-anatomical-boundary'
    # Recipient names do not change a vessel's representation: a tonsillar
    # branch is an artery, and an artery of a sulcus is not the sulcus itself.
    # Keep explicit openings, immaterial sites and wall coats above.
    vascular_name = re.search(r'\b(?:artery|arteries|vein|veins|branch|branches|anastomosis|plexus)\b', name)
    vascular_latin = re.match(r'^(?:a|aa|v|vv|r|rr)\.\s', row.get('name_la','').lower())
    if code.startswith(('A12.2.','A12.3.')) and (vascular_name or vascular_latin or (code.startswith('A12.3.') and 'sinus' in name)):
        return 'vessel-or-vascular-network'
    if any(x in name for x in ('capillar','arteriole','venule','sinusoid','vasa vasorum','vessels of nerves','lymphoid nodule','white pulp','red pulp','germinal','cords','glandula')):return 'microscopic-tissue-or-network'
    if any(x in name for x in ('surface','impression','area','bare area')):return 'surface-region'
    if any(x in name for x in ('margin','border','sulcus','fissure','groove')):return 'boundary-curve-or-groove'
    if any(x in name for x in ('valve','valvule','leaflet','cusp')):return 'valve-or-leaflet'
    if code.startswith('A13.')or code=='A12.4.02.001'or any(x in name for x in ('tonsil','lymphoid','lymph node','thym','spleen','splenic pulp','bone marrow')):
        return 'lymphoid-organ-tissue-or-node-group'
    return 'vessel-or-vascular-network'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--resolved',type=Path,default=RESOLVED)
    parser.add_argument('--report',type=Path,default=REPORT)
    parser.add_argument('--document',type=Path,default=DOC)
    parser.add_argument('--json-only',action='store_true',help='Save a fresh machine-readable snapshot without another document')
    parser.add_argument('--expected-source-sha256',default='1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39')
    args=parser.parse_args()
    source=args.source.resolve();resolved=args.resolved.resolve();report=args.report.resolve();doc=args.document.resolve()
    began=time.monotonic();assert not report.exists()and(args.json_only or not doc.exists()),'Never overwrite an existing inventory'
    source_sha=SHA(source);assert source_sha==args.expected_source_sha256,'Explicit expected catalogue revision changed'
    with sqlite3.connect('file:'+str(AUTH)+'?mode=ro',uri=True)as conn:
        conn.row_factory=sqlite3.Row;terms=corrected_terms({r['id']:dict(r)for r in conn.execute('select * from ta98')})
        synonyms=[dict(r)for r in conn.execute('select id,synonym,lang from synonyms')]
    children=collections.defaultdict(list)
    for code,row in terms.items():children[row.get('parent_id')].append(code)
    def closure(root):
        found={root};stack=[root]
        while stack:
            for code in children[stack.pop()]:
                if code not in found:found.add(code);stack.append(code)
        return found
    heart=closure('A12.1.00.001')
    # User scope is blood circulation only. Chapter 12 also contains lymphatic
    # terminology; a chapter number alone does not authorize lymphatic work.
    lymphatic_roots={'A12.0.00.038','A12.0.00.043','A12.0.00.044','A12.0.00.045','A12.4.01.001'}
    lymphatic=set().union(*(closure(code)for code in lymphatic_roots))
    lymphoid=closure('A13.0.00.000')
    scope=closure('A12.0.00.000')-heart-lymphatic-lymphoid
    data=json.loads(source.read_bytes());effective=json.loads(resolved.read_bytes());assert effective['sourceSha256']==source_sha
    visibility=validate_effective_visibility(data,effective)
    assert effective['hierarchySha256']==SHA(ROOT/'public/assets/ta98.json'),'Resolver snapshot uses a stale hierarchy'
    assert effective['metadataSha256']==SHA(ROOT/'app/data/ta98-metadata.json'),'Resolver snapshot uses stale source identities'
    nodes={n['id']:n for n in effective['resolved']};parts={p['id']:p for p in data['parts']};hidden=set(data.get('defaultHidden',[]))
    bindings=collections.defaultdict(list)
    effective_bindings=collections.defaultdict(list);raw_bindings=collections.defaultdict(list)
    for key,stage,target in [('rawConceptBindings','raw-source-concept-and-authority-metadata',raw_bindings),
                             ('reviewedConceptBindings','reviewed-source-concept-and-authority-metadata',bindings),
                             ('conceptBindings','effective-concept-and-authority-metadata',effective_bindings)]:
        for b in effective[key]:
            if b['term']in scope:target[b['term']].append(dict(b,binding_source=stage))
    for key,stage,target in [('sourceTopLevelBindings','source-ta98TermBindings',bindings),('topLevelBindings','effective-ta98TermBindings',effective_bindings)]:
        for b in effective[key]:
            if b['term']in scope:target[b['term']].append(dict(term=b['term'],conceptId=b.get('conceptId'),partIds=b['partIds'],kind=b.get('kind','exact'),method='explicit-top-level-binding',binding_source=stage))
    for code,records in effective_bindings.items():bindings[code].extend(records)
    typed=collections.defaultdict(list)
    def scan(value,path):
        if isinstance(value,dict):
            term=value.get('term')or value.get('ta98Term')
            if term in scope and value.get('representation'):
                surface=value.get('surface',{});anchor=value.get('anchor')
                typed[term].append(dict(path=path,representation=value['representation'],partIds=value.get('partIds',[]),
                    surfaceFaceCount=len(surface.get('faces',[])),anchor=anchor,curveCount=len(value.get('curves',[])),anatomicalAcceptance=False))
            for k,v in value.items():scan(v,path+'/'+str(k))
        elif isinstance(value,list):
            for i,v in enumerate(value):scan(v,path+'/'+str(i))
    scan({k:v for k,v in data.items()if k not in {'parts','concepts','materials'}},'catalogue')
    names=collections.defaultdict(set)
    for p in parts.values():names[normal(p['name'])].add(p['id'])
    authority_names={code:{normal(terms[code]['name_en']),normal(terms[code]['name_la'])}for code in scope}
    for row in synonyms:
        if row['id']in scope and row['lang'].lower()in {'en','eng','english'}:authority_names[row['id']].add(normal(row['synonym']))
    rows=[];part_terms=collections.defaultdict(set);set_terms=collections.defaultdict(list)
    for code in sorted(scope):
        r=terms[code];node=nodes.get(code,{});direct=sorted(set(node.get('directPartIds',[])))
        exact={pid for b in bindings[code]if b['kind']=='exact'for pid in b['partIds']if pid in parts}
        parent={pid for b in bindings[code]if b['kind']=='parent'for pid in b['partIds']if pid in parts}
        effective_exact={pid for b in effective_bindings[code]if b['kind']=='exact'for pid in b['partIds']if pid in parts}
        effective_parent={pid for b in effective_bindings[code]if b['kind']=='parent'for pid in b['partIds']if pid in parts}
        active=lambda ids:sorted(pid for pid in ids if visibility[pid]['available'])
        inactive=lambda ids:sorted(pid for pid in ids if not visibility[pid]['available'])
        candidates=set().union(*(names.get(name,set())for name in authority_names[code]))
        # Exact name equivalence is a retrieval lead only, never an assignment.
        unbound=sorted(candidates-set(direct)-exact-parent)
        sex=term_sex(terms,code)
        row=dict(id=code,name_en=r['name_en'],name_la=r['name_la'],parent_id=r['parent_id'],parent_name=terms.get(r['parent_id'],{}).get('name_en'),
            source_parent_id=r.get('source_parent_id',r['parent_id']),authority_type=r['type_of_entity'],authority_immaterial=bool(r['immaterial']),
            variant=bool(r['variant']),bilaterality=bool(r['bilaterality']),sex_qualification=sex,male_applicable=sex!='female',female_applicable=sex!='male',
            scope='blood-circulation',
            planned_representation_class=classification(r,children),classification_is_planning_only=True,
            active_exact_source_part_ids=active(effective_exact),active_parent_mapping_part_ids=active(effective_parent),suppressed_exact_source_part_ids=inactive(exact),
            suppressed_parent_source_part_ids=inactive(parent),resolved_direct_part_ids=direct,
            resolved_descendant_only_part_ids=sorted(set(node.get('descendantPartIds',[]))-set(direct)),
            named_but_unbound_source_candidates=unbound,bindings=bindings[code],raw_source_bindings=raw_bindings[code],effective_bindings=effective_bindings[code],typed_feature_controls=typed[code],
            anatomy_accepted=False,wall_lumen_patency_and_organ_interfaces_verified=False)
        if direct:
            for pid in direct:part_terms[pid].add(code)
            set_terms[tuple(direct)].append(code)
        rows.append(row)
    referenced=set()
    for row in rows:
        code=row['id'];direct=row['resolved_direct_part_ids'];aliases=set_terms.get(tuple(direct),[])if direct else[]
        row['identical_carrier_set_other_terms']=[x for x in aliases if x!=code]
        row['carrier_parts_shared_with_other_terms']={p:sorted(part_terms[p]-{code})for p in direct if part_terms[p]-{code}}
        bounded=any(x['representation']=='point-landmark'and x.get('anchor')or x['representation']=='surface-region'and x['surfaceFaceCount']or x['representation']=='boundary-curve'and x['curveCount']for x in row['typed_feature_controls'])
        if direct:status='typed-feature-carrier-reference'if bounded else'shared-carrier-alias'if len(aliases)>1 else'source-mesh-or-network-reference'
        elif row['active_exact_source_part_ids']:status='exact-source-not-resolved-in-current-hierarchy'
        elif row['active_parent_mapping_part_ids']:status='parent-placement-only'
        elif row['suppressed_exact_source_part_ids']:status='suppressed-source-only'
        elif row['named_but_unbound_source_candidates']:status='named-unbound-candidate-only'
        else:status='missing-direct-representation'
        row['representation_status']=status
        row['missing_requirements']=['Actual geometry/reference review, neighboring packing and vessel/organ interfaces remain unaccepted.']
        if not direct:row['missing_requirements'].insert(0,'No directly selectable active representation resolves for this term in the current TA98 hierarchy.')
        if row['identical_carrier_set_other_terms']and not bounded:row['missing_requirements'].append('This control shares the entire same carrier set with another term; distinct term geometry is not established.')
        for key in ['active_exact_source_part_ids','active_parent_mapping_part_ids','suppressed_exact_source_part_ids','suppressed_parent_source_part_ids','resolved_direct_part_ids','named_but_unbound_source_candidates']:referenced.update(row[key])
        for binding in row['raw_source_bindings']:referenced.update(pid for pid in binding['partIds']if pid in parts)
    bank_hashes={};source_records=[]
    for pid in sorted(referenced):
        p=parts[pid];url=data['chunks'][p['chunk']]['url'];assert url.startswith('/local-models/')
        bank=ROOT/'.local-models'/url.removeprefix('/local-models/')
        if url not in bank_hashes:bank_hashes[url]=SHA(bank)
        fields={k:p[k]for k in ['id','name','conceptId','system','positions','normals','normalSigns','indices','vertexCount','indexCount','geometryOrigin','sourceOffset','geometryScale','geometryInstance','bounds','license','provenance','coverageLimitation']if k in p}
        source_records.append(dict(**fields,chunk_url=url,bank_sha256=bank_hashes[url],suppressed=visibility[pid]['suppressed'],default_hidden=visibility[pid]['defaultHidden'],
            ordinary_view_available=visibility[pid]['available'],raw_source_suppressed=bool(p.get('suppressed')),raw_source_default_hidden=pid in hidden or bool(p.get('defaultHidden'))))
    byid={r['id']:r for r in rows}
    def counts(codes):
        chosen=[byid[x]for x in codes];male=[r for r in chosen if r['male_applicable']]
        count=lambda rs:sum(bool(r['resolved_direct_part_ids'])for r in rs)
        return dict(all_terms=len(chosen),male_applicable_terms=len(male),active_direct_reference_terms=count(chosen),male_active_direct_reference_terms=count(male),
            male_no_direct_reference_terms=len(male)-count(male),male_reference_availability_percent=round(100*count(male)/len(male),3)if male else 0,
            male_ninety_percent_target=math.ceil(.9*len(male)),male_status_counts=dict(collections.Counter(r['representation_status']for r in male)),accepted_anatomy_terms=0)
    branches=[]
    for code in sorted(scope):
        sub=closure(code)&scope
        if len(sub)>=10:branches.append(dict(root=code,name=terms[code]['name_en'],**counts(sub)))
    summary={'blood_circulation_excluding_heart':counts(scope),
        'arteries':counts(closure('A12.2.00.001')),'veins':counts(closure('A12.3.00.001'))}
    excluded=[dict(id=k,name_en=terms[k]['name_en'],parent_id=terms[k]['parent_id'])for k in sorted(heart|{k for k in terms if k.startswith('A12.1.')})]
    result=dict(schema='ta98-blood-circulation-authoritative-inventory/v3',view=effective['view'],source_catalogue=str(source.relative_to(ROOT)),source_catalogue_sha256=source_sha,
        authority_cache=str(AUTH.relative_to(ROOT)),authority_cache_sha256=SHA(AUTH),frontend_resolver_export=str(resolved.relative_to(ROOT)),frontend_resolver_sha256=SHA(resolved),
        authority_roots=dict(cardiovascular='A12.0.00.000',excluded_heart='A12.1.00.001',excluded_lymphatic=sorted(lymphatic_roots),excluded_lymphoid='A13.0.00.000'),
        scope_policy='User-requested blood circulation only, excluding heart anatomy, lymphatic vessels/ducts and lymphoid organs/nodes. Coronary arteries and cardiac veins remain blood-vessel scope. All sex-qualified blood terms remain in rows; male applicability only defines the male subtotal.',
        denominator_not_reduced_for_missing_microstructures_variants_containers_or_aliases=True,
        primary_sources=['https://ifaa.unifr.ch/Public/EntryPage/PDF/TA98%20Chapter%2012.pdf'],
        summary=summary,branches=branches,terms=rows,source_mesh_records=source_records,bank_hashes=bank_hashes,excluded_heart_associated_terms=excluded,
        excluded_non_blood_terms=[dict(id=k,name_en=terms[k]['name_en'])for k in sorted(lymphatic|lymphoid)],
        literal_carrier_aliases_do_not_establish_distinct_representation=True,all_geometry_and_interface_acceptance_pending=True,
        no_geometry_catalogue_ui_or_authority_modified=True,elapsed_seconds=time.monotonic()-began)
    assert SHA(source)==source_sha and not(scope&(lymphatic|lymphoid|heart))
    assert closure('A12.2.00.001')|closure('A12.3.00.001')<=scope
    report.parent.mkdir(parents=True,exist_ok=True)
    # This census deliberately runs on stdlib Python. Serialize completely
    # before opening the new output, reject NaN, and retain exclusive creation.
    report_bytes=(json.dumps(result,indent=2,allow_nan=False)+'\n').encode('utf-8')
    with report.open('xb')as stream:stream.write(report_bytes)
    if args.json_only:
        print(json.dumps(dict(report=str(report),terms=len(rows),summary=summary,source_mesh_records=len(source_records),bank_count=len(bank_hashes),elapsed_seconds=time.monotonic()-began)))
        return
    lines=['# TA98 blood-circulation inventory — male source snapshot','',
        'This complete source and control census establishes availability and missing work; it accepts no vessel geometry, wall, lumen, blood-flow continuity or organ interface.',
        '',f'Male source SHA256: `{source_sha}`. Full machine-readable term, binding, suppressed alternative and source-bank evidence: [{report.relative_to(ROOT)}]({report.relative_to(ROOT)}).',
        '', '## Scope and denominator','',
        'The full canonical cardiovascular root `A12.0.00.000` has 1,301 descendants including itself. The canonical heart subtree `A12.1.00.001` contributes 101 excluded terms. Twelve further pericardial chapter terms have canonical mediastinal/pericardial parents and are also recorded as excluded heart-associated context. Coronary arteries and cardiac veins remain in the vascular inventory.',
        '', 'The requested scope is blood circulation only. Lymphatic vessels, ducts, valves, capillaries and the separate lymphoid organ/node system are excluded from this work and its denominator. No missing blood-circulation microscopic feature, variant, container, shared-carrier alias or female-specific term is silently deleted.',
        '', '| Scope | All sex terms | Male applicable | Active direct references | Male without direct reference | 90% male target |', '|---|---:|---:|---:|---:|---:|']
    for name,c in summary.items():lines.append(f'| {name} | {c["all_terms"]} | {c["male_applicable_terms"]} | {c["male_active_direct_reference_terms"]} | {c["male_no_direct_reference_terms"]} | {c["male_ninety_percent_target"]} |')
    lines.extend(['','The direct-reference numerator is an **availability upper bound**, not a count of correct distinct anatomy. Whole-carrier aliases are marked; parent-only placement is not exact identity. Named source candidates without a binding are listed separately rather than called absent meshes. The frontend resolver includes reviewed identities, source metadata/name mappings, top-level `ta98TermBindings`, installed typed controls and the actual TA98 hierarchy.',
        '', 'Representation classes are explicit construction-planning categories based on authoritative entity type and anatomical names. They are not proof that the required layer, microscopic tissue or opening has been constructed. A point requires a source-bound point, a surface requires a bounded surface, a boundary requires a curve, and an opening must not become a solid plug.',
        '', '## Largest missing branches','', '| Root | Branch | Male terms | Direct references | No direct reference |','|---|---|---:|---:|---:|'])
    for b in sorted(branches,key=lambda b:b['male_no_direct_reference_terms'],reverse=True)[:20]:lines.append(f'| {b["root"]} | {b["name"]} | {b["male_applicable_terms"]} | {b["male_active_direct_reference_terms"]} | {b["male_no_direct_reference_terms"]} |')
    lines.extend(['','## Every applicable authority term','', 'The versionable table lists every canonical parent, exact/parent/direct source identity, suppressed alternative, named unbound candidate and whole-carrier alias. All source-bank hashes and full binding evidence are retained in the local JSON. Female-only terms required for the later transfer remain listed. A dash means no source of that type.','', '| Term | Name | Canonical parent | Planned representation | Male | Active exact sources | Parent placement sources | Resolved direct sources | Suppressed sources | Named unbound candidates | Carrier-set aliases | Status |','|---|---|---|---|---|---|---|---|---|---|---|---|'])
    compact=lambda values:'<br>'.join('`'+x.replace('|','/')+'`'for x in values)or'—'
    for r in rows:lines.append('| '+' | '.join([r['id'],r['name_en'].replace('|','/'),r['parent_id'],r['planned_representation_class'],'yes'if r['male_applicable']else'female-only',compact(r['active_exact_source_part_ids']),compact(r['active_parent_mapping_part_ids']),compact(r['resolved_direct_part_ids']),compact(sorted(set(r['suppressed_exact_source_part_ids']+r['suppressed_parent_source_part_ids']))),compact(r['named_but_unbound_source_candidates']),compact(r['identical_carrier_set_other_terms']),r['representation_status']])+' |')
    lines.extend(['','## Primary terminology references','', '[IFAA TA98 chapter 12](https://ifaa.unifr.ch/Public/EntryPage/PDF/TA98%20Chapter%2012.pdf) identifies the cardiovascular terminology. The local authority cache and explicit reviewed parent corrections supply the enumerated blood-circulation closure; lymphatic entries are excluded by the user-requested scope.',''])
    with doc.open('x')as f:f.write('\n'.join(lines))
    print(json.dumps(dict(report=str(report),document=str(doc),terms=len(rows),summary=summary,source_mesh_records=len(source_records),bank_count=len(bank_hashes),elapsed_seconds=time.monotonic()-began)))


if __name__=='__main__':main()
