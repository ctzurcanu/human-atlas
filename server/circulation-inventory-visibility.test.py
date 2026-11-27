"""Ordinary-view availability cannot be inferred from raw source flags."""
import copy
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('circulation_inventory',ROOT/'server/build-circulation-inventory.py')
inventory=importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory)


def fixture():
    source={'parts':[{'id':'donor'},{'id':'owned-view'},{'id':'finite-wall'}],
            'defaultHidden':['donor','owned-view']}
    visibility=[dict(id=pid,available=pid!='donor',suppressed=pid=='donor',
                     defaultHidden=pid!='finite-wall',partDefaultHidden=False,
                     rawSuppressed=False,rawDefaultHidden=pid!='finite-wall')
                for pid in ['donor','owned-view','finite-wall']]
    snapshot=dict(schema='circulation-effective-frontend-source-census/v2',
                  view=dict(mode='ordinary-local-atlas',version=1,inspectAlternatives=False,
                            pipeline=['withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)']),
                  effectivePartVisibility=visibility,
                  resolved=[dict(id='term',directPartIds=['owned-view'],descendantPartIds=['owned-view','finite-wall'])],
                  rawConceptBindings=[],reviewedConceptBindings=[],conceptBindings=[],sourceTopLevelBindings=[],topLevelBindings=[])
    return source,snapshot


class OrdinaryVisibility(unittest.TestCase):
    def test_effective_suppression_retires_donor_but_not_verified_hidden_control(self):
        source,snapshot=fixture()
        self.assertNotIn('suppressed',source['parts'][0])
        actual=inventory.validate_effective_visibility(source,snapshot)
        self.assertFalse(actual['donor']['available'])
        self.assertTrue(actual['owned-view']['defaultHidden'])
        self.assertTrue(actual['owned-view']['available'])

    def test_old_or_repair_view_snapshots_are_explicitly_rejected(self):
        source,snapshot=fixture()
        for mutate in [lambda s:s.update(schema='circulation-effective-frontend-source-census/v1'),
                       lambda s:s.pop('view'),lambda s:s['view'].update(inspectAlternatives=True),
                       lambda s:s['view'].update(pipeline=['withReviewedCatalogueIdentities'])]:
            with self.subTest(mutate=mutate):
                changed=copy.deepcopy(snapshot);mutate(changed)
                with self.assertRaisesRegex(ValueError,'Ordinary-view census requires'):
                    inventory.validate_effective_visibility(source,changed)

    def test_missing_duplicate_or_inconsistent_visibility_cannot_enter_counts(self):
        source,snapshot=fixture()
        for mutate in [lambda s:s['effectivePartVisibility'].pop(),
                       lambda s:s['effectivePartVisibility'].append(s['effectivePartVisibility'][0]),
                       lambda s:s['effectivePartVisibility'][0].update(available=True),
                       lambda s:s['effectivePartVisibility'][0].update(rawDefaultHidden=False)]:
            with self.subTest(mutate=mutate):
                changed=copy.deepcopy(snapshot);mutate(changed)
                with self.assertRaises(ValueError):inventory.validate_effective_visibility(source,changed)

    def test_unavailable_direct_and_descendant_ids_fail_instead_of_inflating_counts(self):
        source,snapshot=fixture()
        for key in ['directPartIds','descendantPartIds']:
            changed=copy.deepcopy(snapshot);changed['resolved'][0][key].append('donor')
            with self.assertRaisesRegex(ValueError,'Resolved hierarchy contains unavailable source donor'):
                inventory.validate_effective_visibility(source,changed)


if __name__=='__main__':unittest.main()
