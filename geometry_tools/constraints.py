"""Explicit mechanical relationships, forest closure and geometric checks."""
from pathlib import Path
import json
import numpy as np
from psycopg.types.json import Jsonb
from .common import digest_file,identity
from .registry import load_registered
from .contact import pair_contacts,require_common_frame
from .spatial import BVH
from .audit import topology,winding_contains

KINDS={'fixed','shared-interface','disjoint','nested','adjacent-unjoined','sliding','open-conduit','clearance'}


def import_forest(db,path,model,*,fixed_systems=('skeletal',)):
    path=Path(path).resolve();forest=json.loads(path.read_bytes());forest_sha=digest_file(path);groups=[]
    def walk(node):
        if 'partIds'in node:groups.append(node)
        else:
            for child in node.get('children',[]):walk(child)
    for root in forest['forest']:walk(root)
    with db.connect()as c:
        cat=c.execute('SELECT a.catalogue_sha256,c.metadata FROM atlas_geom.current_catalogues a JOIN atlas_geom.catalogues c ON c.sha256=a.catalogue_sha256 WHERE a.model_id=%s',(model,)).fetchone()
        if not cat or cat['metadata']['json_sha256']!=forest['sourceSha256']:raise ValueError('Dependency forest is stale relative to the imported literal catalogue')
        rows=c.execute('SELECT * FROM atlas_geom.meshes WHERE catalogue_sha256=%s',(cat['catalogue_sha256'],)).fetchall();byid={r['part_id']:r for r in rows}
        members=[pid for g in groups for pid in g['partIds']]
        if len(members)!=len(set(members))or set(members)!=set(byid):raise ValueError('Forest must own every catalogue item exactly once')
        for group in groups:
            gid=identity(cat['catalogue_sha256'],'assembly',group['id'])
            c.execute('INSERT INTO atlas_geom.assemblies VALUES(%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',
                      (gid,cat['catalogue_sha256'],group['name'],group.get('basis','explicit-forest-ownership'),Jsonb(dict(forest_group_id=group['id'],source_sha256=forest_sha,move_scale_together=group.get('moveScaleTogether',True),source_frame=group.get('sourceFrame')))))
            active=set(group.get('activePartIds',[]))
            for pid in group['partIds']:c.execute('INSERT INTO atlas_geom.assembly_members VALUES(%s,%s,%s)ON CONFLICT DO NOTHING',(gid,byid[pid]['id'],pid in active))
        fixed=0
        for row in rows:
            if row['system']in fixed_systems and not row['suppressed']:
                cid=identity('fixed-forest',row['id']);evidence=dict(source=str(path),source_sha256=forest_sha,reason='Explicit fixed-system policy supplied to forest import')
                c.execute('INSERT INTO atlas_geom.constraints(id,kind,first_mesh_id,evidence,parameters)VALUES(%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',(cid,'fixed',row['id'],Jsonb(evidence),Jsonb({'geometry_sha256':row['geometry_sha256']})));fixed+=1
        imported=0;declarations=[];unresolved_endpoints=[]
        mapping={'shared-tissue-interface':'shared-interface','attachment-refit-obligation':'adjacent-unjoined','disjoint-tissue-refit-obligation':'disjoint','anatomical-open-conduit':'open-conduit','open-conduit':'open-conduit','adjacent-unjoined':'adjacent-unjoined','nested':'nested','sliding':'sliding','clearance':'clearance','disjoint':'disjoint'}
        for interface in forest.get('interfaces',[]):
            a,b=interface.get('from'),interface.get('to');kind=mapping.get(interface['kind'])
            if a not in byid or b not in byid:
                unresolved_endpoints.append(interface);continue
            evidence=dict(source=str(path),source_sha256=forest_sha,forest_interface=interface,verification='Imported declaration; geometric verification remains pending')
            cid=identity('forest-interface',byid[a]['id'],byid[b]['id'],interface)
            declarations.append((cid,cat['catalogue_sha256'],byid[a]['id'],byid[b]['id'],interface['kind'],bool(interface.get('coMove',False)),Jsonb(evidence)))
            if not kind:continue
            c.execute('INSERT INTO atlas_geom.constraints(id,kind,first_mesh_id,second_mesh_id,evidence,parameters)VALUES(%s,%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',(cid,kind,byid[a]['id'],byid[b]['id'],Jsonb(evidence),Jsonb(dict(co_move=interface.get('coMove',False),verification_pending=True))));imported+=1
        if declarations:
            with c.cursor()as cursor:cursor.executemany('INSERT INTO atlas_geom.relation_declarations(id,catalogue_sha256,first_mesh_id,second_mesh_id,kind,co_move,evidence)VALUES(%s,%s,%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',declarations)
    return dict(model=model,assembly_count=len(groups),owned_items=len(members),fixed_items=fixed,imported_interface_declarations=imported,total_forest_relations_preserved=len(declarations),unresolved_relation_endpoints=unresolved_endpoints,geometric_interfaces_verified=False,anatomical_acceptance=False)


def add_constraint(db,spec,model=None):
    if spec.get('kind')not in KINDS:raise ValueError('Unknown constraint kind')
    if not isinstance(spec.get('evidence'),dict)or not spec['evidence'].get('source'):raise ValueError('A constraint requires an explicit evidence source')
    a=db.mesh(spec['first'],model);b=db.mesh(spec['second'],model)if spec.get('second')else None
    if spec['kind']!='fixed'and not b:raise ValueError('Two mesh identities are required')
    if b and a['frame_id']!=b['frame_id']:raise ValueError('Constraint endpoints occupy different declared frames')
    params=spec.get('parameters',{});cid=identity('declared-constraint/v1',spec['kind'],a['id'],b['id']if b else None,spec['evidence'],params)
    with db.connect()as c:
        c.execute('INSERT INTO atlas_geom.constraints(id,kind,first_mesh_id,second_mesh_id,evidence,parameters)VALUES(%s,%s,%s,%s,%s,%s)ON CONFLICT DO NOTHING',(cid,spec['kind'],a['id'],b['id']if b else None,Jsonb(spec['evidence']),Jsonb(params)))
    return dict(id=cid,kind=spec['kind'],first=a['id'],second=b['id']if b else None)


def refine(db,cid,spec):
    if not isinstance(spec.get('evidence'),dict)or not spec['evidence'].get('source'):raise ValueError('Refinement requires an explicit evidence source')
    parameters=spec['parameters']
    if parameters.get('verification_pending'):raise ValueError('A refinement must provide actual reviewed geometric parameters')
    with db.connect()as c:
        old=c.execute('SELECT * FROM atlas_geom.constraints WHERE id=%s FOR UPDATE',(cid,)).fetchone()
        if not old:raise ValueError('Unknown constraint')
        evidence=dict(spec['evidence'],previous_declarations=[*old['evidence'].get('previous_declarations',[]),dict(evidence=old['evidence'],parameters=old['parameters'])])
        c.execute('UPDATE atlas_geom.constraints SET evidence=%s,parameters=%s WHERE id=%s',(Jsonb(evidence),Jsonb(parameters),cid))
    return dict(id=cid,kind=old['kind'],refined=True,geometric_check_required=True)


def closure(db,references,model=None):
    initial={db.mesh(x,model)['id']for x in references};selected=set(initial)
    with db.connect()as c:
        # Only explicit ownership and co_move declarations propagate movement.
        # AABB overlap or a near-contact result never becomes an attachment.
        while True:
            before=len(selected)
            members=c.execute('SELECT DISTINCT b.mesh_id FROM atlas_geom.assembly_members a JOIN atlas_geom.assembly_members b USING(assembly_id) WHERE a.mesh_id=ANY(%s)',(list(selected),)).fetchall()
            selected.update(x['mesh_id']for x in members)
            links=c.execute("SELECT first_mesh_id,second_mesh_id FROM atlas_geom.constraints WHERE parameters->>'co_move'='true' AND (first_mesh_id=ANY(%s)OR second_mesh_id=ANY(%s))",(list(selected),list(selected))).fetchall()
            for link in links:selected.update(x for x in link.values()if x)
            links=c.execute('SELECT first_mesh_id,second_mesh_id FROM atlas_geom.relation_declarations WHERE co_move AND (first_mesh_id=ANY(%s)OR second_mesh_id=ANY(%s))',(list(selected),list(selected))).fetchall()
            for link in links:selected.update(link.values())
            if len(selected)==before:break
        rows=c.execute('SELECT id,part_id,name,geometry_sha256,hidden,suppressed FROM atlas_geom.meshes WHERE id=ANY(%s) ORDER BY part_id',(list(selected),)).fetchall()
        fixed=c.execute("SELECT first_mesh_id AS id,evidence FROM atlas_geom.constraints WHERE kind='fixed'AND first_mesh_id=ANY(%s)",(list(selected),)).fetchall()
    return dict(initial_ids=sorted(initial),members=rows,fixed=fixed,movable=not fixed,dependency_expansion=len(selected)-len(initial),proximity_inferred=False)


def list_constraints(db,model=None,mesh_ids=None):
    with db.connect()as c:
        return c.execute('''SELECT r.* FROM atlas_geom.constraints r JOIN atlas_geom.meshes m ON m.id=r.first_mesh_id
          LEFT JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=m.catalogue_sha256
          WHERE (%s::text IS NULL OR a.model_id=%s)AND (%s::text[]IS NULL OR first_mesh_id=ANY(%s)OR second_mesh_id=ANY(%s)) ORDER BY r.id''',(model,model,mesh_ids,mesh_ids,mesh_ids)).fetchall()


def _literal_cyclic_equal(first_fields,second_fields):
    """Exact oriented facet multisets, including repeated positions/words.

    Sort numeric columns without casting mixed integer/float corner fields.
    This has the same literal semantics as the tuple Counter, while avoiding
    Python objects and whole-mesh triangle copies for every small seam.
    """
    def canonical(fields):
        columns=[np.stack([field[:,(corner+r)%3]for r in range(3)],axis=1)
                 for corner in range(3)for field in fields]
        phase=np.lexsort(tuple(columns[::-1]),axis=-1)[:,0]
        rows=[np.take_along_axis(c,phase[:,None],axis=1)[:,0]for c in columns]
        order=np.lexsort(tuple(rows[::-1]));return [c[order]for c in rows]
    first=canonical(first_fields);second=canonical(second_fields)
    return all(np.array_equal(a,b)for a,b in zip(first,second))


def check_constraint(spec,a,b=None,original=None):
    kind=spec['kind'];p=spec.get('parameters',{})
    if kind=='fixed':
        expected=p.get('geometry_sha256')or(original.geometry_sha if original is not None else None)
        return dict(kind=kind,passed=expected is not None and expected==a.geometry_sha,expected=expected,actual=a.geometry_sha)
    require_common_frame(a,b)
    if p.get('verification_pending'):return dict(kind=kind,passed=False,pending=True,reason='Declaration lacks reviewed geometric parameters')
    if kind=='shared-interface':
        ai=a.faces[np.asarray(p['first_faces'],int)];bi=b.faces[np.asarray(p['second_faces'],int)][:,[0,2,1]]
        af=a.vertices[ai];bf=b.vertices[bi]
        first=[af[:,:,k]for k in range(3)];second=[bf[:,:,k]for k in range(3)]
        matched=bool(len(af))and _literal_cyclic_equal(first,second);normal_match=False
        if 'normals'in a.attributes and 'normals'in b.attributes:
            def words(mesh,indices):
                n=mesh.attributes['normals'][indices].astype('int64')
                if 'normalSigns'in mesh.attributes:n=n*mesh.attributes['normalSigns'][indices][:,:,None]
                return n
            aw=words(a,ai);bw=-words(b,bi)
            normal_match=_literal_cyclic_equal(first+[aw[:,:,k]for k in range(3)],second+[bw[:,:,k]for k in range(3)])
        return dict(kind=kind,passed=matched and normal_match,literal_opposed_triangles=matched,literal_opposed_normal_words=normal_match)
    if kind=='open-conduit':
        from geometry_tools.algorithms.tube_continuity import match_mouths
        aa=a.vertices[np.asarray(p['first_ring'],int)][None];bb=b.vertices[np.asarray(p['second_ring'],int)][None]
        order,r=match_mouths(aa,bb,position_tolerance_m=p.get('position_tolerance_m',0.),area_tolerance_relative=p.get('area_tolerance_relative',.001))
        from .model import Mesh
        ring=aa[0];center=ring.mean(0);normal=np.cross(ring-center,np.roll(ring,-1,axis=0)-center).sum(0);normal=normal/np.linalg.norm(normal)
        if np.max(abs((ring-center)@normal))>p.get('planarity_tolerance_m',1e-9):return dict(kind=kind,passed=False,reason='Nonplanar mouth requires an explicit reviewed gate surface')
        turns=np.einsum('ij,j->i',np.cross(np.roll(ring,-1,axis=0)-ring,np.roll(ring,-2,axis=0)-np.roll(ring,-1,axis=0)),normal)
        if np.any(turns<=0):return dict(kind=kind,passed=False,reason='Nonconvex/collinear mouth requires an explicit reviewed triangulation')
        gate=Mesh(np.r_[ring,center[None]],np.array([[len(ring),k,(k+1)%len(ring)]for k in range(len(ring))]),a.frame);blocks=[]
        for mesh in [a,b]:
            arr,summary=pair_contacts(mesh,gate,0.)
            blocks.append(dict(obstructing_face_pairs=int(np.isin(arr['relation'],['crossing','coplanar-overlap']).sum()),complete=not summary['truncated']))
        return dict(kind=kind,passed=all(x['obstructing_face_pairs']==0 and x['complete']for x in blocks),mouth=r,mouth_patency=blocks,whole_lumen_obstruction_checked=False)
    ba,bb=BVH(a),BVH(b);nearest=ba.closest_mesh(bb);arr,contact=pair_contacts(a,b,0.,first_bvh=ba,second_bvh=bb,max_pairs=p.get('max_pairs',250000))
    crossings=contact['crossing_face_pairs'];overlaps=contact['coplanar_overlap_face_pairs'];passed=not contact['truncated']and crossings==0
    out=dict(kind=kind,minimum_distance_m=float(nearest['distance_m']),crossings=crossings,coplanar_overlaps=overlaps,complete=not contact['truncated'])
    if kind=='nested':
        if not topology(b)['closed_oriented_manifold']:return {**out,'passed':False,'reason':'Container is not a closed oriented volume'}
        # With exhaustive no-crossing evidence, each connected surface remains
        # in one component of the complement. Use one seed per component.
        used=np.unique(a.faces);inside=winding_contains(b,a.vertices[used[:1]])if len(topology(a)['components'])==1 else winding_contains(b,a.vertices[used])
        passed=passed and all(inside['inside'])and not any(inside['boundary_ambiguous']);out['containment']=inside
    elif kind=='disjoint'and p.get('occupancy')=='solid':
        if not topology(a)['closed_oriented_manifold']or not topology(b)['closed_oriented_manifold']:return {**out,'passed':False,'reason':'Solid separation requires two closed volumes'}
        for host,other in [(a,b),(b,a)]:
            probes=other.vertices[np.unique(other.faces)]
            inside=winding_contains(host,probes[:1]if len(topology(other)['components'])==1 else probes)
            passed=passed and not any(inside['inside'])and not any(inside['boundary_ambiguous'])
    if overlaps and not p.get('allow_coplanar_contact',False):passed=False
    if kind in ['clearance','sliding','adjacent-unjoined']:
        low=p.get('minimum_m',0.);high=p.get('maximum_m',float('inf'));passed=passed and low<=out['minimum_distance_m']<=high
    out['passed']=bool(passed);return out


def check(db,model=None,candidates=None,*,parameter_overrides=None,endpoint_overrides=None):
    candidates=candidates or {};rows=list_constraints(db,model,list(candidates)or None);loaded={}
    def get(mid):
        if mid in candidates:return candidates[mid]
        if mid not in loaded:loaded[mid]=load_registered(db,mid)[1]
        return loaded[mid]
    results=[]
    for row in rows:
        actual=dict(row)
        if row['id']in(parameter_overrides or{}):actual['parameters']=parameter_overrides[row['id']]
        try:
            if row['id']in(endpoint_overrides or{}):
                if row['kind']!='shared-interface':raise ValueError('Only shared tissue interface ownership may be transferred')
                first,second=endpoint_overrides[row['id']]
                r=check_constraint(actual,first,second)
            else:r=check_constraint(actual,get(row['first_mesh_id']),get(row['second_mesh_id'])if row['second_mesh_id']else None)
        except (ValueError,KeyError)as e:r=dict(passed=False,reason=str(e))
        results.append(dict(id=row['id'],**r))
    return dict(passed=all(r['passed']for r in results),constraints=results,anatomical_acceptance=False)


def add_parser(sub):
    p=sub.add_parser('constraints');s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('import-forest');q.add_argument('forest');q.add_argument('--model',required=True);q.add_argument('--fixed-system',action='append',default=['skeletal'])
    q=s.add_parser('add');q.add_argument('spec');q.add_argument('--model')
    q=s.add_parser('refine');q.add_argument('id');q.add_argument('spec')
    q=s.add_parser('closure');q.add_argument('meshes',nargs='+');q.add_argument('--model')
    q=s.add_parser('list');q.add_argument('--model')
    q=s.add_parser('relations');q.add_argument('--model',required=True);q.add_argument('--kind');q.add_argument('--mesh')
    q=s.add_parser('check');q.add_argument('--model')
    p.set_defaults(handler=run)


def run(args,db):
    if args.action=='import-forest':return import_forest(db,args.forest,args.model,fixed_systems=args.fixed_system)
    if args.action=='add':return add_constraint(db,json.loads(Path(args.spec).read_text()),args.model)
    if args.action=='refine':return refine(db,args.id,json.loads(Path(args.spec).read_text()))
    if args.action=='closure':return closure(db,args.meshes,args.model)
    if args.action=='list':return dict(constraints=list_constraints(db,args.model))
    if args.action=='relations':
        mid=db.mesh(args.mesh,args.model)['id']if args.mesh else None
        with db.connect()as c:
            rows=c.execute('''SELECT r.* FROM atlas_geom.relation_declarations r JOIN atlas_geom.current_catalogues a ON a.catalogue_sha256=r.catalogue_sha256
            WHERE a.model_id=%s AND (%s::text IS NULL OR r.kind=%s)AND (%s::text IS NULL OR r.first_mesh_id=%s OR r.second_mesh_id=%s) ORDER BY r.id''',(args.model,args.kind,args.kind,mid,mid,mid)).fetchall()
        return dict(relations=rows,geometric_constraint_inferred=False)
    return check(db,args.model)
