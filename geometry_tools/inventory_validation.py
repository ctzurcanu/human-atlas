"""Authoritative inventory integrity, independent of anatomical completion.

This bounded, read-only gate accepts a full inventory with explicitly unfinished
representations. It does not turn donor meshes, mappings, or caller booleans into
verified geometry or ordinary-view coverage.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import sqlite3

from .ta98_authority import CONTEXT_REFERENCES, corrected_terms, term_sex

ROOT = Path(__file__).resolve().parents[1]
LEGAL_CLASSES = frozenset({'tissue-region', 'surface-region', 'wall-layer', 'opening',
                          'boundary-curve', 'point-landmark', 'microscopic-structure', 'whole-organ',
                          'luminal-content'})


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _roots(inventory):
    roots = []
    for key in ['root_term', 'additional_context_root']:
        value = inventory.get(key)
        if value is not None:
            if not isinstance(value, str) or not value:
                raise ValueError('Inventory root must be a nonempty authority code: ' + key)
            roots.append(value)
    for key in ['root_terms', 'roots', 'additional_context_roots', 'whole_source_roots']:
        values = inventory.get(key, [])
        if not isinstance(values, list) or any(not isinstance(v, str) or not v for v in values):
            raise ValueError('Inventory roots must be authority-code arrays: ' + key)
        roots.extend(values)
    if not roots:
        raise ValueError('Declare authoritative inventory root_term/root_terms before installation')
    return sorted(set(roots))


def validate_authoritative_inventory(inventory, *, model, current_catalogue,
                                     authority_path=None, corrections_path=None):
    """Return exact scope/gap evidence; never infer anatomical coverage.

    All declared roots include their full corrected descendants. Additional
    explicit context terms preserve canonical parents. The historical heart
    cross-reference flag means the full TA98 heart chapter occurrence set,
    including its pericardial entries under canonical parents outside the heart.
    """
    sex = {'local-ta98': 'male', 'local-female-ta98': 'female'}.get(model)
    if not sex or inventory.get('schema') != 'ta98-geometry-inventory/v1' or inventory.get('model') != model:
        raise ValueError('A source-bound TA98 geometry inventory for this local model is required')
    if inventory.get('model_sex', sex) != sex:
        raise ValueError('Inventory model sex disagrees with the actual model')
    authority_path = Path(authority_path or ROOT / '.local-models/ta98-authority.sqlite').resolve()
    corrections_path = Path(corrections_path or Path(__file__).with_name('ta98_authority.py')).resolve()
    authority_sha, corrections_sha = _sha(authority_path), _sha(corrections_path)
    declared_authority = Path(inventory['authority_path']).resolve()
    if inventory.get('authority_sha256') != authority_sha or _sha(declared_authority) != authority_sha:
        raise ValueError('Inventory authority must match the actual current authoritative cache')
    bindings = {str(authority_path): authority_sha, str(declared_authority): authority_sha,
                str(corrections_path): corrections_sha}
    if inventory.get('authority_corrections_sha256') not in [None, corrections_sha]:
        raise ValueError('Inventory uses stale authority corrections')
    if inventory.get('authority_corrections_path'):
        declared_corrections = Path(inventory['authority_corrections_path']).resolve()
        if _sha(declared_corrections) != corrections_sha:
            raise ValueError('Inventory corrections do not match the current reviewed corrections')
        bindings[str(declared_corrections)] = corrections_sha
    current_catalogue = Path(current_catalogue).resolve()
    current_sha = _sha(current_catalogue)
    source_hashes = [inventory[k] for k in ['source_catalogue_sha256', 'current_catalogue_sha256'] if inventory.get(k)]
    source = inventory.get('source_catalogue')
    if isinstance(source, dict):
        if source.get('sha256'):
            source_hashes.append(source['sha256'])
        source = source.get('path')
    if not source_hashes or any(value != current_sha for value in source_hashes):
        raise ValueError('Inventory must bind the exact current catalogue revision; refresh the source inventory')
    if source is not None:
        if not isinstance(source, str) or _sha(source) != current_sha:
            raise ValueError('Inventory source catalogue path does not match its current source binding')
        bindings[str(Path(source).resolve())] = current_sha
    bindings[str(current_catalogue)] = current_sha
    with sqlite3.connect('file:' + str(authority_path) + '?mode=ro', uri=True) as conn:
        conn.row_factory = sqlite3.Row
        authority = corrected_terms({row['id']: dict(row) for row in conn.execute('SELECT * FROM ta98')})
    roots = _roots(inventory)
    if any(root not in authority for root in roots):
        raise ValueError('Inventory contains an unknown authoritative root')
    children = {}
    for code, row in authority.items():
        children.setdefault(row.get('parent_id'), []).append(code)
    required = set()
    for root in roots:
        pending = [root]
        while pending:
            code = pending.pop()
            if code in required:
                continue
            required.add(code)
            pending.extend(children.get(code, []))
        required.update(CONTEXT_REFERENCES.get(root, {}))
    context = inventory.get('additional_context_terms', [])
    if not isinstance(context, list) or any(not isinstance(code, str) or code not in authority for code in context):
        raise ValueError('Explicit context terms must name actual authoritative rows')
    required.update(context)
    if inventory.get('includes_pericardial_cross_references'):
        if 'A12.1.00.001' not in roots:
            raise ValueError('Pericardial chapter cross references require the declared heart root')
        required.update(code for code in authority if code.startswith('A12.1.'))
    if not required <= authority.keys():
        raise ValueError('An explicit authoritative context reference is missing from the cache')
    rows = inventory.get('terms')
    if not isinstance(rows, list) or not rows:
        raise ValueError('A complete nonempty authoritative inventory is required')
    codes = [row.get('authority', {}).get('id') if isinstance(row, dict) and isinstance(row.get('authority'), dict) else None for row in rows]
    if any(not isinstance(code, str) or not code for code in codes):
        raise ValueError('Every inventory row must contain its authoritative code')
    if len(set(codes)) != len(codes):
        raise ValueError('Inventory contains duplicate authoritative term rows')
    if type(inventory.get('term_count')) is not int or inventory['term_count'] != len(rows):
        raise ValueError('Inventory term_count must equal its actual unique row count')
    missing, extra = sorted(required - set(codes)), sorted(set(codes) - required, key=str)
    if missing or extra:
        raise ValueError('Inventory is not the full declared authority closure: missing=' + str(missing[:12]) + '; extra=' + str(extra[:12]))
    applicable = []
    missing_representations, provisional, unverified, declared_verified = [], [], [], []
    unspecified = []
    for row in rows:
        code = row['authority']['id']
        actual = authority[code]
        # Require the entire authoritative row, including canonical parents,
        # sex/type metadata and explicit reviewed parent corrections. Extra
        # inventory annotations may live alongside it, never replace it.
        if any(key not in row['authority'] or row['authority'][key] != value for key, value in actual.items()):
            raise ValueError('Inventory authoritative row/parent metadata differs from the current cache: ' + code)
        if row.get('representation_class') not in LEGAL_CLASSES:
            raise ValueError('Inventory representation_class must be a legal classified value: ' + code)
        qualified = term_sex(authority, code)
        is_applicable = not qualified or qualified == sex
        if 'applicable_to_model' in row and (type(row['applicable_to_model']) is not bool or row['applicable_to_model'] != is_applicable):
            raise ValueError('Inventory applicability disagrees with authoritative model sex: ' + code)
        if 'authority_sex' in row and row['authority_sex'] != qualified:
            raise ValueError('Inventory authority_sex differs from its canonical ancestry: ' + code)
        if not is_applicable:
            continue
        applicable.append(code)
        if row.get('missing_distinct_representation') is True:
            missing_representations.append(code)
        explicit_provisional = any(row.get(key) is True for key in ['provisional', 'representation_provisional', 'boundary_provisional'])
        planned = row.get('planned_representation')
        if isinstance(planned, dict):
            explicit_provisional |= any(planned.get(key) is True for key in ['provisional', 'boundaryProvisional'])
        if explicit_provisional:
            provisional.append(code)
        verified = row.get('distinct_representation_verified') is True
        if verified and row.get('missing_distinct_representation') is True:
            raise ValueError('Inventory representation is simultaneously declared missing and verified: ' + code)
        if not verified or explicit_provisional:
            unverified.append(code)
        if verified and row.get('missing_distinct_representation') is False and not explicit_provisional:
            declared_verified.append(code)
        if type(row.get('missing_distinct_representation')) is not bool or type(row.get('distinct_representation_verified')) is not bool:
            unspecified.append(code)
    if 'applicable_term_count' in inventory and inventory['applicable_term_count'] != len(applicable):
        raise ValueError('Inventory applicable_term_count differs from authoritative applicability')
    return dict(schema='authoritative-inventory-integrity/v1', passed=True, authority_inventory_complete=True,
                source_catalogue_sha256=current_sha, authority_sha256=authority_sha, authority_corrections_sha256=corrections_sha,
                model=model, model_sex=sex, roots=roots, term_count=len(rows), applicable_term_count=len(applicable),
                declared_missing_representation_terms=sorted(missing_representations), declared_missing_representation_count=len(missing_representations),
                explicitly_provisional_terms=sorted(provisional), explicitly_provisional_count=len(provisional),
                unverified_representation_terms=sorted(unverified), unverified_representation_count=len(unverified),
                declared_verified_representation_terms=sorted(declared_verified), declared_verified_representation_count=len(declared_verified),
                unspecified_representation_status_terms=sorted(unspecified),
                anatomical_coverage_complete=False, anatomical_coverage_evaluated=False,
                ordinary_view_representation_coverage_evaluated=False,
                source_mesh_or_candidate_mapping_does_not_establish_representation=True,
                coverage_qualification='Authority inventory integrity only. Missing/provisional/unverified rows remain permitted for local previews; anatomy and effective-view coverage are not certified.',
                input_bindings=bindings)
