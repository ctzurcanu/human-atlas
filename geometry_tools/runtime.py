"""Meter and bound local commands, including native geometry libraries.

The supervisor lives outside the worker, so a stuck C/C++ call cannot prevent
the deadline from being enforced. Only the newly created process group is
stopped; an existing Blender session is never part of that group.
"""
from pathlib import Path
import argparse
import ctypes
import datetime
import errno
import json
import math
import os
import signal
import subprocess
import sys
import time
import uuid
from .process_memory_observer import process_observer

DEFAULT_SECONDS = 120.0
MAX_SECONDS = 300.0
DEFAULT_MEMORY_MIB = 5120.0
GLOBAL_MEMORY_ALERT_BYTES = 15_000_000_000
GLOBAL_MEMORY_INTERVAL_SECONDS = 1.0
DEADLINE_ENV = 'ATLAS_RUNTIME_DEADLINE_UTC'
MEMORY_ENV = 'ATLAS_RUNTIME_MEMORY_LIMIT_MIB'
REPORT_ENV = 'ATLAS_RUNTIME_REPORT'
ROOT = Path(__file__).resolve().parent.parent


class MemoryMonitorError(RuntimeError):
    """The owned process group's RSS could not be measured reliably."""


class ProcessGroupMemory:
    """Read only the owned POSIX group, without shelling out to ps.

    Summed RSS is conservative: shared pages may be counted for each process.
    On macOS, also account for physical footprint, including compressed memory
    that is absent from RSS. Linux uses the explicit RSS + VmSwap estimate.
    This is sampled enforcement, not a kernel allocation ceiling. Brief peaks
    between samples cannot be recovered after a worker exits.
    """
    def __init__(self):
        self.pending = {}
        if sys.platform == 'darwin':
            self.lib = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True)
            self.lib.proc_listpids.argtypes = [ctypes.c_uint32, ctypes.c_uint32,
                                               ctypes.c_void_p, ctypes.c_int]
            self.lib.proc_listpids.restype = ctypes.c_int
            self.lib.proc_pidinfo.argtypes = [ctypes.c_int, ctypes.c_int,
                                              ctypes.c_uint64, ctypes.c_void_p, ctypes.c_int]
            self.lib.proc_pidinfo.restype = ctypes.c_int
            self.lib.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
            self.lib.proc_pid_rusage.restype = ctypes.c_int
            self.footprint_kind = 'macOS ri_phys_footprint (includes compressed memory)'
        elif sys.platform.startswith('linux'):
            self.page_size = os.sysconf('SC_PAGE_SIZE')
            self.footprint_kind = 'Linux RSS + VmSwap conservative estimate'
        else:
            raise MemoryMonitorError('RSS monitoring supports macOS and Linux only')

    def sample(self, group):
        if sys.platform == 'darwin':
            # PROC_PGRP_ONLY=2. The sizing call may return an all-PID estimate;
            # the filled call filters to the supplied owned group.
            ctypes.set_errno(0)
            estimate = self.lib.proc_listpids(2, group, None, 0)
            if estimate < 0 or (estimate == 0 and ctypes.get_errno() not in (0, errno.ESRCH)):
                raise MemoryMonitorError('Cannot enumerate the owned process group')
            capacity = max(4096, estimate + 64)
            for _ in range(3):
                buffer = ctypes.create_string_buffer(capacity)
                ctypes.set_errno(0)
                count = self.lib.proc_listpids(2, group, buffer, capacity)
                if count < 0 or (count == 0 and ctypes.get_errno() not in (0, errno.ESRCH)):
                    raise MemoryMonitorError('Cannot enumerate the owned process group')
                if count < capacity:
                    break
                capacity *= 2
            else:
                raise MemoryMonitorError('Owned process group changed faster than it could be sampled')
            pids = list((ctypes.c_int * (count // ctypes.sizeof(ctypes.c_int))).from_buffer(buffer))
            rows = []
            for pid in pids:
                if pid <= 0:
                    continue
                # PROC_PIDTASKINFO=4; its first words are virtual_size,
                # resident_size. Use a larger buffer for ABI-compatible tails.
                task = ctypes.create_string_buffer(256)
                ctypes.set_errno(0)
                size = self.lib.proc_pidinfo(pid, 4, 0, task, len(task))
                if size >= 16:
                    resident = int((ctypes.c_uint64 * 2).from_buffer(task)[1])
                    # SDK rusage_info_v0: UUID[16], then ten uint64 fields;
                    # resident is field 6 and physical footprint is field 7.
                    usage = ctypes.create_string_buffer(96)
                    ctypes.set_errno(0)
                    result = self.lib.proc_pid_rusage(pid, 0, usage)
                    if result != 0:
                        if ctypes.get_errno() in (errno.ESRCH, errno.ENOENT):
                            continue
                        raise MemoryMonitorError('Cannot read physical footprint for owned PID ' + str(pid))
                    footprint = int((ctypes.c_uint64 * 10).from_buffer(usage, 16)[7])
                    rows.append({'pid': pid, 'rss_bytes': resident,
                                 'physical_footprint_bytes': footprint,
                                 'accounted_memory_bytes': max(resident, footprint)})
                    self.pending.pop(pid, None)
                elif ctypes.get_errno() not in (0, errno.ESRCH, errno.ENOENT):
                    # A short-lived child can be listed during exec/exit while
                    # task info temporarily returns EPERM. Retain explicit
                    # pending telemetry for at most three sampling intervals;
                    # persistent unreadability still fails closed.
                    since = self.pending.setdefault(pid, time.monotonic())
                    if time.monotonic() - since >= .15:
                        raise MemoryMonitorError('Cannot read RSS for owned PID ' + str(pid))
                    rows.append({'pid': pid, 'rss_bytes': 0,
                                 'physical_footprint_bytes': 0,
                                 'accounted_memory_bytes': 0,
                                 'measurement_pending': True})
                else:
                    self.pending.pop(pid, None)
            observed = set(pids)
            self.pending = {pid: since for pid, since in self.pending.items() if pid in observed}
            return self._totals(rows)
        rows = []
        for path in Path('/proc').iterdir():
            if not path.name.isdigit():
                continue
            try:
                text = (path / 'stat').read_text()
                fields = text[text.rfind(')') + 2:].split()
                if int(fields[2]) != group:
                    continue
                resident = int(fields[21]) * self.page_size
                swap_lines = [line for line in (path / 'status').read_text().splitlines() if line.startswith('VmSwap:')]
                swap = int(swap_lines[0].split()[1]) * 1024 if swap_lines else 0
                footprint = resident + swap
                rows.append({'pid': int(path.name), 'rss_bytes': resident,
                             'physical_footprint_bytes': footprint,
                             'accounted_memory_bytes': max(resident, footprint)})
            except (FileNotFoundError, ProcessLookupError):
                continue
            except PermissionError as error:
                raise MemoryMonitorError('Cannot inspect /proc process-group membership') from error
        return self._totals(rows)

    @staticmethod
    def _totals(rows):
        return {'rss_bytes': sum(r['rss_bytes'] for r in rows),
                'physical_footprint_bytes': sum(r['physical_footprint_bytes'] for r in rows),
                'accounted_memory_bytes': sum(r['accounted_memory_bytes'] for r in rows),
                'processes': rows}


def _memory_limit(requested, environment):
    if not math.isfinite(requested) or requested <= 0:
        raise ValueError('Memory limit must be finite and greater than zero')
    inherited = float(environment.get(MEMORY_ENV, requested))
    if not math.isfinite(inherited) or inherited <= 0:
        raise ValueError('Invalid inherited memory limit')
    return min(requested, inherited)


def _termination_reason(now, deadline, rss_bytes, memory_bytes):
    if now >= deadline:
        return 'timed_out'
    if rss_bytes >= memory_bytes:
        return 'memory_limit'
    return None


def _observe_global_memory(observer, record):
    """Capture high-memory identities while alive; never signal them here."""
    rows, unreadable = observer.sample()
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    record['global_memory_sample_count'] += 1
    record['global_memory_latest_at'] = stamp
    record['global_memory_readable_process_count'] = len(rows)
    record['global_memory_unreadable_process_count'] = len(unreadable)
    record['global_memory_unreadable_processes'] = unreadable
    ordered = sorted(rows, key=lambda row:max(row['rss_bytes'], row['physical_footprint_bytes']), reverse=True)
    record['global_memory_current_top'] = ordered[:10]
    if ordered:
        peak = record.get('global_memory_peak_process')
        top = ordered[0]
        if peak is None or max(top['rss_bytes'], top['physical_footprint_bytes']) > max(peak['rss_bytes'], peak['physical_footprint_bytes']):
            record['global_memory_peak_process'] = dict(top, sampled_at=stamp)
    alerts = []
    for row in ordered:
        if max(row['rss_bytes'], row['physical_footprint_bytes']) < GLOBAL_MEMORY_ALERT_BYTES:
            break
        alerts.append(dict(row, sampled_at=stamp, argv=observer.command(row['pid']),
                           parent_chain=observer.ancestry(row['pid']),
                           threshold_bytes=GLOBAL_MEMORY_ALERT_BYTES,
                           observed_process_signalled=False))
    if alerts:
        record['observed_memory_alerts'] = alerts
    return alerts


def _save(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def _signal_group(process, number):
    try:
        os.killpg(process.pid, number)
    except ProcessLookupError:
        pass
    except PermissionError:
        # On macOS an already-reaped group can return EPERM rather than ESRCH.
        # Ignore only after both root termination and an actual empty-group
        # sample are established. A live/unreadable/nonempty group still raises.
        if process.poll() is not None and not ProcessGroupMemory().sample(process.pid)['processes']:
            return
        raise


def run(command, *, seconds=DEFAULT_SECONDS, report=None, cwd=None,
        heartbeat_seconds=15.0, label=None, memory_mib=DEFAULT_MEMORY_MIB,
        observe_all_processes=False):
    """Run once; retain live/final timing evidence and stop the owned group.

    ``seconds`` includes the termination grace period. Nested commands inherit
    the remaining deadline, never a fresh full budget. Exit 124 means timeout.
    Exit 137 means sampled owned-group memory exceeded its cap; opt-in global
    observation uses exit 138 when another accessible process reaches 15 GB.
    Exit 125 means memory monitoring failed. Only the owned group's resource
    cap is enabled by default; all-process observation requires an explicit opt-in.
    No retry, deadline extension, installation or acceptance is performed.
    """
    if not math.isfinite(seconds) or not 0 < seconds <= MAX_SECONDS:
        raise ValueError('Runtime limit must be greater than zero and at most 300 seconds')
    if not command or not math.isfinite(heartbeat_seconds) or heartbeat_seconds <= 0:
        raise ValueError('A command and a positive heartbeat interval are required')
    if os.name != 'posix':
        raise RuntimeError('This local supervisor requires POSIX process groups')
    memory_limit_mib = _memory_limit(memory_mib, os.environ)
    memory_limit_bytes = int(memory_limit_mib * 1024 * 1024)
    started = time.monotonic()
    utc_started = time.time()
    deadline = utc_started + seconds
    if DEADLINE_ENV in os.environ:
        inherited = float(os.environ[DEADLINE_ENV])
        if not math.isfinite(inherited):
            raise ValueError('Invalid inherited runtime deadline')
        deadline = min(deadline, inherited)
    budget = max(0.0, deadline - utc_started)
    monotonic_deadline = started + budget
    grace = min(0.5, budget * 0.1)
    if report is None:
        report = ROOT / '.local-models/geometry-workbench/runtime' / (uuid.uuid4().hex + '.json')
    report = Path(report).resolve()
    report.parent.mkdir(parents=True, exist_ok=True)
    # Never replace evidence belonging to an earlier invocation.
    with report.open('x'):
        pass
    record = dict(schema='bounded-command/v1', label=label or Path(command[0]).name,
                  command=list(command), cwd=str(Path(cwd or os.getcwd()).resolve()),
                  started_at=datetime.datetime.fromtimestamp(utc_started, datetime.timezone.utc).isoformat(),
                  limit_seconds=budget, requested_limit_seconds=seconds,
                  deadline_utc=deadline, elapsed_seconds=0.0, state='starting',
                  requested_memory_limit_mib=memory_mib, memory_limit_mib=memory_limit_mib,
                  memory_scope='sum of per-process max(RSS, physical footprint) in the newly owned process group; shared pages may be counted repeatedly',
                  memory_sample_interval_seconds=0.05, memory_sample_count=0,
                  current_rss_bytes=0, peak_sampled_rss_bytes=0,
                  current_physical_footprint_bytes=0, peak_sampled_physical_footprint_bytes=0,
                  current_accounted_memory_bytes=0, peak_sampled_accounted_memory_bytes=0,
                  global_memory_alert_bytes=GLOBAL_MEMORY_ALERT_BYTES,
                  global_memory_interval_seconds=GLOBAL_MEMORY_INTERVAL_SECONDS,
                  global_memory_sample_count=0,
                  global_memory_observation_enabled=observe_all_processes,
                  global_memory_observation=('all accessible processes; macOS physical footprint/RSS, Linux RSS+VmSwap estimate'
                                             if observe_all_processes else 'disabled'),
                  observed_unrelated_processes_signalled=False,
                  attempts=1, automatic_retry=False, report=str(report))
    process = None
    old_handlers = {}

    def interrupted(number, frame):
        raise KeyboardInterrupt('Supervisor interrupted by signal ' + str(number))

    try:
        _save(report, record)
        if budget <= 0:
            record.update(state='timed_out', returncode=124, worker_started=False)
            return record
        try:
            monitor = ProcessGroupMemory()
            preflight = monitor.sample(os.getpgrp())
            if preflight['rss_bytes'] <= 0 or os.getpid() not in [row['pid'] for row in preflight['processes']]:
                raise MemoryMonitorError('RSS preflight did not observe this supervisor in its process group')
            record['physical_footprint_measurement'] = monitor.footprint_kind
            global_monitor = process_observer() if observe_all_processes else None
            if global_monitor is not None and _observe_global_memory(global_monitor, record):
                record.update(state='observed_memory_alert', returncode=138,
                              worker_started=False)
                return record
        except (OSError, RuntimeError) as error:
            record.update(state='memory_monitor_failed', returncode=125,
                          worker_started=False, error=str(error))
            return record
        environment = dict(os.environ, **{DEADLINE_ENV: str(deadline), REPORT_ENV: str(report),
                                          MEMORY_ENV: str(memory_limit_mib)})
        for number in (signal.SIGTERM, signal.SIGINT):
            old_handlers[number] = signal.signal(number, interrupted)
        process = subprocess.Popen(command, cwd=cwd, env=environment, start_new_session=True)
        record.update(state='running', pid=process.pid, worker_started=True)
        _save(report, record)
        print(f"[runtime] {record['label']}: limit {budget:.1f}s, owned-group memory {memory_limit_mib:g} MiB; {report}", file=sys.stderr, flush=True)
        next_heartbeat = started + heartbeat_seconds
        next_global_sample = time.monotonic() + GLOBAL_MEMORY_INTERVAL_SECONDS
        terminating = False
        stop_reason = None
        while process.poll() is None:
            now = time.monotonic()
            record['elapsed_seconds'] = now - started
            if now >= monotonic_deadline:
                terminating = True
                stop_reason = 'timed_out'
                _signal_group(process, signal.SIGKILL)
                break
            try:
                memory = monitor.sample(process.pid)
            except (OSError, MemoryMonitorError) as error:
                if process.poll() is not None:
                    break
                stop_reason = 'memory_monitor_failed'
                record['error'] = str(error)
                _signal_group(process, signal.SIGKILL)
                break
            record['memory_sample_count'] += 1
            record['current_rss_bytes'] = memory['rss_bytes']
            record['current_physical_footprint_bytes'] = memory['physical_footprint_bytes']
            record['peak_sampled_physical_footprint_bytes'] = max(record['peak_sampled_physical_footprint_bytes'], memory['physical_footprint_bytes'])
            record['current_accounted_memory_bytes'] = memory['accounted_memory_bytes']
            if memory['accounted_memory_bytes'] > record['peak_sampled_accounted_memory_bytes']:
                record['peak_sampled_accounted_memory_bytes'] = memory['accounted_memory_bytes']
                record['peak_sampled_accounted_processes'] = memory['processes']
            if memory['rss_bytes'] > record['peak_sampled_rss_bytes']:
                record['peak_sampled_rss_bytes'] = memory['rss_bytes']
                record['peak_sampled_processes'] = memory['processes']
            reason = _termination_reason(now, monotonic_deadline,
                                         memory['accounted_memory_bytes'], memory_limit_bytes)
            if reason:
                stop_reason = reason
                record['state'] = reason
                _save(report, record)
                _signal_group(process, signal.SIGKILL)
                break
            if global_monitor is not None and now >= next_global_sample:
                try:
                    alerts = _observe_global_memory(global_monitor, record)
                except (OSError, RuntimeError) as error:
                    stop_reason = 'memory_monitor_failed'
                    record['error'] = 'Global memory observation: ' + str(error)
                    _save(report, record)
                    _signal_group(process, signal.SIGKILL)
                    break
                next_global_sample = time.monotonic() + GLOBAL_MEMORY_INTERVAL_SECONDS
                _save(report, record)
                if alerts:
                    stop_reason = 'observed_memory_alert'
                    record['state'] = stop_reason
                    _save(report, record)
                    print('[runtime] Accessible process exceeded 15 GB; stopped only this owned work group. Identity retained in ' + str(report), file=sys.stderr, flush=True)
                    _signal_group(process, signal.SIGKILL)
                    break
            if now >= monotonic_deadline - grace and not terminating:
                terminating = True
                stop_reason = 'timed_out'
                record['state'] = 'stopping_at_deadline'
                _signal_group(process, signal.SIGTERM)
                _save(report, record)
            if now >= next_heartbeat:
                _save(report, record)
                print(f"[runtime] {record['label']}: {record['elapsed_seconds']:.1f}/{budget:.1f}s; accounted memory {record['current_accounted_memory_bytes']/1048576:.1f} MiB, sampled peak {record['peak_sampled_accounted_memory_bytes']/1048576:.1f} MiB", file=sys.stderr, flush=True)
                next_heartbeat = now + heartbeat_seconds
            time.sleep(min(0.05, max(0.0, monotonic_deadline - now)))
        code = process.wait(timeout=1)
        record.update(state=stop_reason or ('completed' if code == 0 else 'failed'),
                      returncode={'timed_out':124, 'memory_limit':137, 'memory_monitor_failed':125, 'observed_memory_alert':138}.get(stop_reason, code), worker_returncode=code)
    except KeyboardInterrupt as error:
        record.update(state='cancelled', returncode=130, error=str(error))
    except BaseException as error:
        record.update(state='failed', returncode=1, error=repr(error))
        raise
    finally:
        if process is not None:
            # Also stop children left behind after the root command exits.
            _signal_group(process, signal.SIGKILL)
            process.wait(timeout=1)
        for number, handler in old_handlers.items():
            signal.signal(number, handler)
        record['elapsed_seconds'] = time.monotonic() - started
        _save(report, record)
        print(f"[runtime] {record['label']}: {record['state']} after {record['elapsed_seconds']:.2f}s", file=sys.stderr, flush=True)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=DEFAULT_SECONDS)
    parser.add_argument('--report')
    parser.add_argument('--label')
    parser.add_argument('--cwd')
    parser.add_argument('--memory-mib', type=float, default=DEFAULT_MEMORY_MIB)
    parser.add_argument('--observe-all-processes', action='store_true',
                        help='Explicitly enable the otherwise-disabled system-wide 15 GB observer')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('Supply a command after --')
    return run(command, seconds=args.seconds, report=args.report, cwd=args.cwd, label=args.label,
               memory_mib=args.memory_mib, observe_all_processes=args.observe_all_processes)['returncode']


def geometry_main():
    """Supervise the ordinary geometry CLI before importing native libraries."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--wall-seconds', type=float, default=DEFAULT_SECONDS)
    parser.add_argument('--runtime-report')
    parser.add_argument('--memory-mib', type=float, default=DEFAULT_MEMORY_MIB)
    parser.add_argument('--observe-all-processes', action='store_true')
    args, rest = parser.parse_known_args()
    command = [sys.executable, '-B', '-c',
               'from geometry_tools.cli import main; raise SystemExit(main())', *rest]
    print('Geometry runtime: --wall-seconds N (default 120, maximum 300); --runtime-report FILE', file=sys.stderr, flush=True)
    return run(command, seconds=args.wall_seconds, report=args.runtime_report,
               label='geometry ' + (' '.join(rest[:2]) or 'CLI'), memory_mib=args.memory_mib,
               observe_all_processes=args.observe_all_processes)['returncode']


if __name__ == '__main__':
    raise SystemExit(main())
