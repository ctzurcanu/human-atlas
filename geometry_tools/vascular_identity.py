"""Apply explicitly reviewed source identities without constructing geometry."""
import json
import sqlite3
from contextlib import closing
from pathlib import Path

from .course_io import sha
from .model import AtlasReader


def apply_reviewed_identities(atlas, source, catalogue, plan_path):
    plan_path = Path(plan_path)
    plan = json.loads(plan_path.read_text())
    assert plan['complete'] and plan['sourceAtlasSha256'] == sha(catalogue)
    assert sha(plan['authorityCachePath']) == plan['authorityCacheSha256']
    assert sha(plan['auditPath']) == plan['auditSha256']
    assert sha(plan['reviewPath']) == plan['reviewSha256']
    audit = json.loads(Path(plan['auditPath']).read_text())
    review = json.loads(Path(plan['reviewPath']).read_text())
    assert audit['complete'] and review['complete']
    assert review['reviewedBy'] == 'Codex literal rendered surface/section inspection'
    assert audit['sourceAtlasSha256'] == review['sourceAtlasSha256'] == sha(catalogue)
    assert all(sha(image['path']) == image['sha256'] for image in review['images'])
    assert 1 <= len(plan['bindings']) <= 16
    parts = {part['id']: part for part in atlas['parts']}
    concepts = {concept['id']: concept for concept in atlas['concepts']}
    old_parts = {part['id']: part for part in source['parts']}
    reader = AtlasReader()
    changed = []
    terms = []
    for binding in plan['bindings']:
        pid, term = binding['partId'], binding['term']
        assert audit['sourcePartId'] == review['sourcePartId'] == pid
        part, old = parts[pid], old_parts[pid]
        assert pid.startswith('BDM:male:') and not old.get('suppressed')
        assert old.get('ta98Term') in (None, term)
        concept = concepts[binding['sourceConceptId']]
        recipient = concepts[binding['recipientConceptId']]
        with closing(sqlite3.connect(plan['authorityCachePath'])) as db:
            db.row_factory = sqlite3.Row
            authority = dict(db.execute('select id,name_en,parent_id,immaterial from ta98 where id=?', (term,)).fetchone())
        assert authority['parent_id'] == recipient['ta98Term'] and not authority['immaterial']
        assert part['conceptId'] == concept['id'] and concept['elements'] == [pid]
        assert not any(any(member.get('partId') == pid for member in collection.get('members', []))
                       for collection in atlas.get('sourceCollectionIdentities', []))
        mesh = reader.read_part(source, old, frame='model:local-ta98:assembled-metres', catalogue_path=Path(catalogue))
        assert mesh.geometry_sha == binding['sourceGeometrySha256'] == audit['sourceGeometrySha256']
        roots = review['sourceRecipientRoots']
        assert len(roots) == audit['sourceComponents'] == binding['sourceComponents']
        assert 0 < binding['maximumRecipientGap_m'] <= 2e-6
        assert {root['component'] for root in roots} == set(range(audit['sourceComponents']))
        for root in roots:
            assert root in audit['mouths']
            matches = [target for target in root['targets'] if target['partId'] in recipient['elements']]
            assert len(matches) == 1 and matches[0]['nearestTriangleCertified']
            assert matches[0]['maximumSurfaceDistance_m'] <= binding['maximumRecipientGap_m']
        tags = dict(ta98Term=term, ta98Kind='exact', coverageLimitation=binding['coverageLimitation'])
        part.update(tags)
        concept.update(tags)
        if concept.get('partOf') != recipient['id']:
            concept.update(sourceClinicalParentConceptId=concept.get('partOf'), partOf=recipient['id'])
        previous = [row for row in atlas['ta98TermBindings'] if row['conceptId'] == concept['id']]
        assert all(row['term'] == term for row in previous)
        atlas['ta98TermBindings'] = [row for row in atlas['ta98TermBindings'] if row['conceptId'] != concept['id']]
        atlas['ta98TermBindings'].append(dict(term=term, conceptId=concept['id'], partIds=[pid],
                                            kind='exact', representation='source-identity'))
        changed.append(pid)
        terms.append(term)
    record = dict(planPath=str(plan_path), planSha256=sha(plan_path), partIds=changed, terms=terms,
                  geometryModified=False, wholeCirculationAccepted=False)
    atlas['bloodSourceIdentityAdditions'] = source.get('bloodSourceIdentityAdditions', []) + [record]
    return record
