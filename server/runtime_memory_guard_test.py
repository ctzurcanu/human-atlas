"""Resource-limit behavior, including an actual small owned-group allocation."""
from pathlib import Path
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
import shutil
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools import runtime


class MemoryGuardTest(unittest.TestCase):
    def setUp(self):
        from unittest.mock import Mock
        self.global_backend = Mock()
        self.global_backend.sample.return_value = ([], [])
        self.global_patch = patch.object(runtime, 'process_observer', return_value=self.global_backend)
        self.global_patch.start()
        self.addCleanup(self.global_patch.stop)

    def test_nested_limit_cannot_be_relaxed(self):
        self.assertEqual(runtime._memory_limit(4096, {runtime.MEMORY_ENV: '512'}), 512)
        self.assertEqual(runtime._memory_limit(128, {runtime.MEMORY_ENV: '512'}), 128)

    def test_invalid_limit_is_rejected(self):
        for value in [0, -1, math.inf, math.nan]:
            with self.assertRaises(ValueError):
                runtime._memory_limit(value, {})
        with self.assertRaises(ValueError):
            runtime._memory_limit(128, {runtime.MEMORY_ENV: 'nan'})

    def test_deadline_wins_if_both_limits_are_reached(self):
        self.assertEqual(runtime._termination_reason(10, 10, 100, 100), 'timed_out')
        self.assertEqual(runtime._termination_reason(9, 10, 100, 100), 'memory_limit')
        self.assertIsNone(runtime._termination_reason(9, 10, 99, 100))

    def test_monitor_failure_never_starts_unbounded_worker(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(runtime, 'ProcessGroupMemory', side_effect=runtime.MemoryMonitorError('denied')):
                with patch.object(runtime.subprocess, 'Popen') as launch:
                    result = runtime.run([sys.executable, '-c', 'print(1)'], seconds=2,
                                         report=Path(directory) / 'failure.json')
            launch.assert_not_called()
            self.assertEqual(result['state'], 'memory_monitor_failed')
            self.assertEqual(result['returncode'], 125)
            self.assertFalse(result['worker_started'])

    def test_native_sampler_sees_current_group(self):
        value = runtime.ProcessGroupMemory().sample(os.getpgrp())
        self.assertGreater(value['rss_bytes'], 0)
        self.assertIn(os.getpid(), [row['pid'] for row in value['processes']])
        self.assertEqual(value['rss_bytes'], sum(row['rss_bytes'] for row in value['processes']))
        self.assertGreater(value['physical_footprint_bytes'], 0)
        self.assertGreaterEqual(value['accounted_memory_bytes'], value['rss_bytes'])
        self.assertGreaterEqual(value['accounted_memory_bytes'], value['physical_footprint_bytes'])

    def test_compressed_footprint_is_enforced_when_rss_is_small(self):
        value = runtime.ProcessGroupMemory._totals([
            dict(pid=1, rss_bytes=16, physical_footprint_bytes=900, accounted_memory_bytes=900),
            dict(pid=2, rss_bytes=20, physical_footprint_bytes=800, accounted_memory_bytes=800)])
        self.assertEqual(value['rss_bytes'], 36)
        self.assertEqual(value['accounted_memory_bytes'], 1700)
        self.assertEqual(runtime._termination_reason(1, 10, value['accounted_memory_bytes'], 1000), 'memory_limit')

    def test_cleanup_permission_error_requires_dead_root_and_empty_group(self):
        from unittest.mock import Mock
        process = Mock(pid=1234)
        process.poll.return_value = 0
        with patch.object(runtime.os, 'killpg', side_effect=PermissionError('denied')):
            with patch.object(runtime, 'ProcessGroupMemory') as sampler:
                sampler.return_value.sample.return_value = {'processes': []}
                runtime._signal_group(process, 9)
                sampler.return_value.sample.return_value = {'processes': [{'pid': 1235}]}
                with self.assertRaises(PermissionError):
                    runtime._signal_group(process, 9)
                process.poll.return_value = None
                sampler.return_value.sample.return_value = {'processes': []}
                with self.assertRaises(PermissionError):
                    runtime._signal_group(process, 9)

    def test_global_alert_captures_identity_and_stops_only_owned_work(self):
        observed = dict(pid=999999, process_start_abstime=123, executable='/Example/ExistingApp',
                        rss_bytes=20*1024*1024, physical_footprint_bytes=16_000_000_000)
        self.global_backend.sample.side_effect = [([], [dict(pid=1, errno=1)]),
                                                  ([observed], [dict(pid=1, errno=1)])]
        self.global_backend.command.return_value = ['/Example/ExistingApp', '--original']
        self.global_backend.ancestry.return_value = [dict(pid=999999, parent_pid=42,
                                                         argv=['/Example/ExistingApp', '--original'])]
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(runtime, 'GLOBAL_MEMORY_INTERVAL_SECONDS', .01):
                with patch.object(runtime, '_signal_group', wraps=runtime._signal_group) as signals:
                    result = runtime.run([sys.executable, '-c', 'import time;time.sleep(2)'],
                                         seconds=1, memory_mib=256,
                                         observe_all_processes=True,
                                         report=Path(directory)/'global-alert.json')
            self.assertEqual(result['state'], 'observed_memory_alert')
            self.assertEqual(result['returncode'], 138)
            self.assertGreaterEqual(result['global_memory_sample_count'], 2)
            self.assertEqual(result['global_memory_unreadable_process_count'], 1)
            event = result['observed_memory_alerts'][0]
            self.assertEqual(event['argv'], ['/Example/ExistingApp', '--original'])
            self.assertEqual(event['parent_chain'][0]['parent_pid'], 42)
            self.assertFalse(event['observed_process_signalled'])
            self.assertTrue(signals.call_args_list)
            for call in signals.call_args_list:
                self.assertEqual(call.args[0].pid, result['pid'])
                self.assertNotEqual(call.args[0].pid, observed['pid'])

    def test_real_global_inventory_is_readable_without_signalling(self):
        from geometry_tools.process_memory_observer import process_observer
        observer = process_observer()
        rows, unreadable = observer.sample()
        self.assertIn(os.getpid(), [row['pid'] for row in rows])
        self.assertTrue(all(isinstance(row['pid'], int) for row in unreadable))
    @unittest.skipUnless(shutil.which('node'), 'Node is required for the small actual allocation check')
    def test_actual_node_footprint_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            evidence_dir = Path(os.environ.get('ATLAS_RUNTIME_TEST_EVIDENCE_DIR', directory))
            evidence_dir.mkdir(parents=True, exist_ok=True)
            path = evidence_dir / 'node-footprint-limit.json'
            result = runtime.run([shutil.which('node'), '-e',
                                  'const x=Buffer.alloc(96*1024*1024,1);setTimeout(()=>console.log(x.length),8000)'],
                                 seconds=4, memory_mib=64, report=path)
            self.assertEqual(result['state'], 'memory_limit')
            self.assertEqual(result['returncode'], 137)
            self.assertGreaterEqual(result['peak_sampled_accounted_memory_bytes'], 64 * 1024 * 1024)
            self.assertGreater(result['peak_sampled_physical_footprint_bytes'], 0)
            self.assertLess(result['elapsed_seconds'], 4)

    def test_actual_group_allocation_is_stopped_and_other_group_survives(self):
        with tempfile.TemporaryDirectory() as directory:
            # A separate, test-owned group represents an existing application.
            sentinel = subprocess.Popen([sys.executable, '-c', 'import time;time.sleep(8)'],
                                        start_new_session=True)
            try:
                path = Path(directory) / 'allocation.json'
                worker = ('import subprocess,sys,time;'
                          'p=subprocess.Popen([sys.executable,"-c",'
                          '"import time;x=bytearray(96*1024*1024);time.sleep(8)"]);'
                          'time.sleep(8)')
                result = runtime.run([sys.executable, '-c', worker], seconds=4,
                                     memory_mib=64, report=path)
                self.assertEqual(result['state'], 'memory_limit')
                self.assertEqual(result['returncode'], 137)
                self.assertGreaterEqual(result['peak_sampled_rss_bytes'], 64 * 1024 * 1024)
                self.assertLess(result['elapsed_seconds'], 4)
                self.assertGreater(result['memory_sample_count'], 0)
                self.assertGreaterEqual(len(result['peak_sampled_processes']), 2)
                self.assertIsNone(sentinel.poll())
                self.assertEqual(json.loads(path.read_bytes())['state'], 'memory_limit')
            finally:
                sentinel.terminate()
                sentinel.wait(timeout=2)

    def test_small_worker_completes_with_recorded_peak(self):
        with tempfile.TemporaryDirectory() as directory:
            result = runtime.run([sys.executable, '-c', 'import time;time.sleep(.12)'],
                                 seconds=2, memory_mib=256,
                                 report=Path(directory) / 'small.json')
            self.assertEqual(result['state'], 'completed')
            self.assertEqual(result['returncode'], 0)
            self.assertGreater(result['peak_sampled_rss_bytes'], 0)
            self.assertLess(result['peak_sampled_rss_bytes'], 256 * 1024 * 1024)

    def test_default_worker_does_not_start_system_wide_observation(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(runtime, 'GLOBAL_MEMORY_INTERVAL_SECONDS', .01):
                result = runtime.run([sys.executable, '-c', 'import time;time.sleep(.12)'],
                                     seconds=2, memory_mib=256,
                                     report=Path(directory) / 'no-global-observer.json')
            runtime.process_observer.assert_not_called()
            self.assertEqual(result['state'], 'completed')
            self.assertFalse(result['global_memory_observation_enabled'])
            self.assertEqual(result['global_memory_sample_count'], 0)


if __name__ == '__main__':
    unittest.main()
