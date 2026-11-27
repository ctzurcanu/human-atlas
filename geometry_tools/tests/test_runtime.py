import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

from geometry_tools.runtime import DEADLINE_ENV, run


class RuntimeTests(unittest.TestCase):
    def test_completion_and_failure_are_metered_without_retry(self):
        with tempfile.TemporaryDirectory() as folder:
            for code, state in [(0, 'completed'), (7, 'failed')]:
                path = Path(folder) / f'{code}.json'
                result = run([sys.executable, '-c', f'raise SystemExit({code})'], seconds=2, report=path)
                self.assertEqual((result['state'], result['returncode']), (state, code))
                self.assertEqual(json.loads(path.read_text()), result)
                self.assertEqual(result['attempts'], 1)
                self.assertGreater(result['elapsed_seconds'], 0)

    def test_native_style_unresponsive_worker_is_killed_at_deadline(self):
        with tempfile.TemporaryDirectory() as folder:
            result = run([sys.executable, '-c', 'import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(20)'],
                         seconds=.5, report=Path(folder)/'timeout.json')
            self.assertEqual(result['state'], 'timed_out')
            self.assertEqual(result['returncode'], 124)
            self.assertLess(result['elapsed_seconds'], 2)

    def test_nested_budget_cannot_restart_the_clock(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {DEADLINE_ENV: str(time.time()+.25)}):
            result = run([sys.executable, '-c', 'import time; time.sleep(20)'], seconds=3, report=Path(folder)/'nested.json')
            self.assertEqual(result['state'], 'timed_out')
            self.assertLess(result['limit_seconds'], .3)
            self.assertLess(result['elapsed_seconds'], 2)

    def test_root_exit_does_not_leave_its_background_worker_running(self):
        with tempfile.TemporaryDirectory() as folder:
            marker = Path(folder)/'escaped.txt'
            code = "import subprocess,sys; subprocess.Popen([sys.executable,'-c',sys.argv[1]])"
            child = f"import time,pathlib; time.sleep(.6); pathlib.Path({str(marker)!r}).write_text('still running')"
            result = run([sys.executable, '-c', code, child], seconds=2, report=Path(folder)/'parent.json')
            time.sleep(.7)
            self.assertEqual(result['state'], 'completed')
            self.assertFalse(marker.exists())

    def test_invalid_budget_or_reused_report_cannot_launch(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'prior.json'
            path.write_text('prior evidence')
            for seconds in (0, -1, 301, float('nan'), float('inf')):
                with self.assertRaises(ValueError):
                    run([sys.executable], seconds=seconds, report=path)
            with self.assertRaises(FileExistsError):
                run([sys.executable], seconds=1, report=path)
            self.assertEqual(path.read_text(), 'prior evidence')

    def test_cancelled_supervisor_stops_owned_worker(self):
        with tempfile.TemporaryDirectory() as folder:
            report = Path(folder)/'cancel.json'
            process = subprocess.Popen([sys.executable, '-B', '-m', 'geometry_tools.runtime', '--seconds', '5', '--report', str(report),
                                        '--', sys.executable, '-c', 'import time; time.sleep(20)'])
            deadline = time.monotonic()+3
            while time.monotonic() < deadline:
                if report.exists():
                    try:
                        if json.loads(report.read_text()).get('state') == 'running':
                            break
                    except json.JSONDecodeError:
                        pass
                time.sleep(.02)
            process.send_signal(signal.SIGTERM)
            self.assertEqual(process.wait(timeout=2), 130)
            self.assertEqual(json.loads(report.read_text())['state'], 'cancelled')

    def test_timed_out_job_status_requires_terminal_evidence_for_its_worker(self):
        from geometry_tools.jobs import stopped_runtime
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'runtime.json'
            parameters=dict(owner_pid=12345,runtime_report=str(path))
            for state in ('running','timed_out','cancelled','completed'):
                path.write_text(json.dumps(dict(state=state,pid=12345,elapsed_seconds=12.5)))
                result=stopped_runtime(parameters)
                self.assertEqual(result is not None,state!='running')
                if result:self.assertEqual(result['elapsed_seconds'],12.5)
            path.write_text(json.dumps(dict(state='timed_out',pid=999)))
            self.assertIsNone(stopped_runtime(parameters))


if __name__ == '__main__':
    unittest.main()
