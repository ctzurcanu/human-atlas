import sqlite3
import unittest

from geometry_tools.common import ROOT
from geometry_tools.ta98_authority import corrected_terms,CONTEXT_REFERENCES,term_sex


class Ta98UrinaryAuthorityTests(unittest.TestCase):
    def terms(self):
        path=ROOT/'.local-models/ta98-authority.sqlite'
        if not path.is_file():self.skipTest('Local primary extract unavailable')
        with sqlite3.connect('file:'+str(path)+'?mode=ro',uri=True)as db:
            db.row_factory=sqlite3.Row
            return {r['id']:dict(r)for r in db.execute('SELECT * FROM ta98')}

    def test_female_states_follow_female_opening_without_changing_raw_authority(self):
        source=self.terms();before=source['A09.2.03.001'].copy();terms=corrected_terms(source)
        self.assertEqual(source['A09.2.03.001'],before)
        children=lambda parent:{code for code,row in terms.items()if row['parent_id']==parent}
        self.assertTrue({'A09.2.03.001','A09.2.03.002'}<=children('A08.3.01.028F'))
        self.assertFalse({'A09.2.03.001','A09.2.03.002'}&children('A08.3.01.027M'))
        self.assertTrue({'A09.4.02.001','A09.4.02.002'}<=children('A08.3.01.027M'))
        self.assertEqual(terms['A09.2.03.001']['female_gender'],1)
        self.assertEqual(terms['A09.4.02.001']['male_gender'],1)
        self.assertEqual(term_sex(terms,'A09.2.03.001'),'female')
        self.assertEqual(term_sex(terms,'A09.4.02.001'),'male')
        self.assertIsNone(term_sex(terms,'A08.3.01.002'))
        if before['parent_id']!='A08.3.01.028F':self.assertEqual(terms['A09.2.03.001']['source_parent_id'],before['parent_id'])

    def test_rectovesicalis_context_retains_canonical_muscular_parent(self):
        source=self.terms();terms=corrected_terms(source)
        self.assertIn('A04.5.03.017',CONTEXT_REFERENCES['A08.3.01.001'])
        self.assertEqual(terms['A04.5.03.017']['parent_id'],source['A04.5.03.017']['parent_id'])
        self.assertNotEqual(terms['A04.5.03.017']['parent_id'],'A08.3.01.001')

    def test_unexpected_extract_drift_is_rejected(self):
        terms=self.terms();terms['A09.2.03.001']['parent_id']='A05.5.01.001'
        with self.assertRaisesRegex(ValueError,'fresh review'):corrected_terms(terms)


if __name__=='__main__':unittest.main()
