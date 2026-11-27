"""Actual input/output artifact identity must survive a guarded install job."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from geometry_tools.common import digest_file
from geometry_tools.jobs import file_bindings, result_file_bindings, verify_bindings


class ResultBindingsTests(unittest.TestCase):
    def test_repeated_witness_banks_are_hashed_once_and_still_reverified(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve()
            bank=root/'bank.bin';bank.write_bytes(b'actual geometry')
            alias=root/'alias.bin';alias.symlink_to(bank.name)
            witnesses=[dict(source=str(bank),peer=str(alias))for _ in range(1000)]
            with patch('geometry_tools.jobs.digest_file',wraps=digest_file)as measured:
                bindings=file_bindings(witnesses)
            self.assertEqual(bindings,{str(bank):digest_file(bank)})
            self.assertEqual(measured.call_count,1)
            bank.write_bytes(b'changed geometry')
            with self.assertRaisesRegex(ValueError,'changed'):
                verify_bindings(bindings)

    def test_audited_layout_retained_when_also_an_input(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve()
            layout=root/'tissue.json'
            bank=root/'cells.bin'
            source=root/'source.npz'
            layout.write_bytes(b'{"origin":[0,1.25,0]}')
            bank.write_bytes(b'actual emitted bytes')
            source.write_bytes(b'ordinary input, not audited output')
            expected={str(p):digest_file(p)for p in [layout,bank]}
            inputs={str(p):digest_file(p)for p in [layout,source]}
            result=dict(audited_banks=list(expected),proposal=dict(bank_bindings=expected),source=str(source))
            self.assertEqual(result_file_bindings(result,inputs),expected)

    def test_audited_artifact_cannot_change_after_the_independent_audit(self):
        with tempfile.TemporaryDirectory() as folder:
            bank=Path(folder)/'cells.bin';bank.write_bytes(b'audited geometry')
            expected={str(bank):digest_file(bank)}
            bank.write_bytes(b'later geometry')
            with self.assertRaisesRegex(ValueError,'missing or changed'):
                result_file_bindings(dict(audited_banks=[str(bank)],proposal=dict(bank_bindings=expected)),expected)

    def test_an_input_alias_is_not_a_new_output_artifact(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve()
            source=root/'source.npz';source.write_bytes(b'unchanged source')
            alias=root/'alias.npz';alias.symlink_to(source.name)
            self.assertEqual(result_file_bindings(dict(source=str(alias)),{str(alias):digest_file(alias)}),{})

    def test_missing_or_undeclared_artifact_cannot_be_promoted(self):
        with tempfile.TemporaryDirectory() as folder:
            bank=Path(folder)/'cells.bin'
            for exists in [False,True]:
                if exists:bank.write_bytes(b'unaudited geometry')
                with self.assertRaisesRegex(ValueError,'missing or changed'):
                    result_file_bindings(dict(audited_banks=[str(bank)],proposal=dict(bank_bindings={})),{})


if __name__=='__main__':
    unittest.main()
