#!/usr/bin/env python3
"""Run one local benchmark, recording RSS and stopping before memory/disk exhaustion."""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--log', type=Path, required=True)
p.add_argument('--resources', type=Path, required=True)
p.add_argument('--max-rss-gib', type=float, default=26)
p.add_argument('--reserve-disk-gib', type=float, default=4)
p.add_argument('command', nargs=argparse.REMAINDER)
a = p.parse_args()
command = a.command[1:] if a.command[:1] == ['--'] else a.command
if not command:
    p.error('missing command')
start = time.monotonic()
with a.log.open('x') as log, a.resources.open('x', newline='') as stats:
    out = csv.writer(stats)
    out.writerow(['elapsed_seconds', 'rss_gib', 'disk_free_gib'])
    child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
    reason = ''
    while child.poll() is None:
        result = subprocess.run(['ps', '-o', 'rss=', '-p', str(child.pid)], capture_output=True, text=True)
        rss = int(result.stdout.strip() or '0') / 1024**2  # macOS ps reports KiB
        free = shutil.disk_usage(os.getcwd()).free / 1024**3
        out.writerow([time.monotonic()-start, rss, free]); stats.flush()
        if rss > a.max_rss_gib or free < a.reserve_disk_gib:
            reason = f'Guard stopped PID {child.pid}: rss={rss:.3f} GiB, free disk={free:.3f} GiB'
            print(reason, flush=True); log.write('\n'+reason+'\n'); log.flush()
            child.terminate()
            try:
                child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                child.kill(); child.wait()
            break
        time.sleep(5)
    code = child.wait()
    print(f'Benchmark exit={code}; elapsed={time.monotonic()-start:.3f}s', flush=True)
    raise SystemExit(3 if reason else code)
