"""Bounded macOS observer: capture every readable process above a memory threshold.

Read-only: this observer never signals or restarts the processes it observes.
Physical footprint includes compressed memory; RSS is retained separately.
"""
import argparse
import ctypes
import ctypes.util
import datetime
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0, str(ROOT))


def atomic(path, value):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, separators=(',', ':'), allow_nan=False) + '\n')
    os.replace(temp, path)


from geometry_tools.process_memory_observer import MacProcesses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=300)
    parser.add_argument('--interval', type=float, default=1)
    parser.add_argument('--threshold-gb', type=float, default=15)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not (0 < args.seconds <= 300 and .2 <= args.interval <= 10 and 0 < args.threshold_gb <= 100):
        parser.error('Require 0–300 seconds, .2–10 second interval and positive finite threshold')
    args.output.mkdir(parents=True, exist_ok=False)
    observer = MacProcesses(); start = time.monotonic(); end = start + args.seconds
    threshold = int(args.threshold_gb * 1_000_000_000); seen = set(); peaks = {}; alerts = 0; samples = 0
    record = dict(schema='process-physical-memory-observer/v1', observer_pid=os.getpid(),
                  started_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  wall_seconds=args.seconds, threshold_bytes=threshold, interval_seconds=args.interval,
                  measurement='macOS physical footprint and RSS; threshold is max of both; decimal GB',
                  signals_or_restarts_observed_processes=False, state='running')
    atomic(args.output/'status.json', record)
    try:
        while time.monotonic() < end:
            rows, unreadable = observer.sample(); samples += 1
            stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
            for row in rows:
                key = (row['pid'], row['process_start_abstime']); usage = max(row['rss_bytes'], row['physical_footprint_bytes'])
                if key not in peaks or usage > max(peaks[key]['rss_bytes'], peaks[key]['physical_footprint_bytes']):
                    peaks[key] = dict(row, sampled_at=stamp)
                if usage > threshold and key not in seen:
                    seen.add(key); alerts += 1
                    event = dict(row, sampled_at=stamp, argv=observer.command(row['pid']),
                                 parent_chain=observer.ancestry(row['pid']),
                                 threshold_bytes=threshold, event='memory-threshold-exceeded')
                    atomic(args.output/f'alert-{alerts:04}-{row["pid"]}.json', event)
                    print(json.dumps(event), flush=True)
            ordered = sorted(rows, key=lambda r:max(r['rss_bytes'], r['physical_footprint_bytes']), reverse=True)
            record.update(samples=samples, alerts=alerts, elapsed_seconds=time.monotonic()-start,
                          latest_at=stamp, readable_processes=len(rows), unreadable_processes=unreadable,
                          current_top=ordered[:15], current_node=[r for r in ordered if 'node' in r['executable'].lower()])
            atomic(args.output/'status.json', record)
            remaining = end-time.monotonic()
            if remaining > 0:
                time.sleep(min(args.interval, remaining))
        record['state'] = 'completed'
    except BaseException as error:
        record.update(state='stopped', error=repr(error))
        raise
    finally:
        record.update(elapsed_seconds=time.monotonic()-start, samples=samples, alerts=alerts)
        atomic(args.output/'status.json', record)
        atomic(args.output/'peaks.json', sorted(peaks.values(), key=lambda r:max(r['rss_bytes'], r['physical_footprint_bytes']), reverse=True))
    print(json.dumps(dict(state=record['state'], samples=samples, alerts=alerts, output=str(args.output))), flush=True)


if __name__ == '__main__':
    main()
