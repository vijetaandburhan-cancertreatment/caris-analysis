"""Run a task-local command with conservative RSS/disk monitoring; no network access."""
import datetime, json, os, pathlib, resource, shutil, signal, subprocess, sys, time

ROOT = pathlib.Path(__file__).resolve().parent
name, *cmd = sys.argv[1:]
if not name or not cmd or '/' in name:
    raise SystemExit('Usage: run_guarded.py job-name executable [args...]')
logs = ROOT / 'logs'
logs.mkdir(exist_ok=True)
start = time.monotonic()
cap_bytes = 6_000_000_000
minimum_free_bytes = 6 * 1024**3
rss_limit_kb = 8 * 1024 * 1024
rec = {'name': name, 'command': cmd, 'cwd': str(pathlib.Path.cwd()),
       'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
       'rss_limit_kb': rss_limit_kb, 'new_disk_cap_bytes': cap_bytes, 'minimum_free_disk_bytes': minimum_free_bytes,
       'max_observed_rss_kb': 0, 'max_observed_dir_bytes': 0}
runtime_env = os.environ.copy()
if '/kallisto-v0.51.1/' in cmd[0]:
    runtime_env['DYLD_LIBRARY_PATH'] = str(ROOT / 'public/h5py-runtime/h5py/.dylibs')
rec['task_local_dynamic_library_path'] = runtime_env.get('DYLD_LIBRARY_PATH')

def dir_bytes():
    return sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file())

with open(logs / f'{name}.stdout.log', 'w') as out, open(logs / f'{name}.stderr.log', 'w') as err:
    proc = subprocess.Popen(cmd, stdout=out, stderr=err, start_new_session=True, env=runtime_env)
    last_disk = 0
    while proc.poll() is None:
        ps = subprocess.run(['ps', '-axo', 'pgid=,rss='], capture_output=True, text=True, check=True)
        rss = sum(int(line.split()[1]) for line in ps.stdout.splitlines()
                  if len(line.split()) == 2 and int(line.split()[0]) == proc.pid)
        rec['max_observed_rss_kb'] = max(rec['max_observed_rss_kb'], rss)
        reason = None
        if rss > rss_limit_kb:
            reason = 'RSS exceeded 8 GiB preferred cap'
        if time.monotonic() - last_disk > 5:
            n = dir_bytes()
            last_disk = time.monotonic()
            rec['max_observed_dir_bytes'] = max(rec['max_observed_dir_bytes'], n)
            if n > cap_bytes:
                reason = 'Additional local directory exceeded 6 GB cap'
            if shutil.disk_usage(ROOT).free < minimum_free_bytes:
                reason = 'Free local disk fell below 6 GiB floor'
        if reason:
            rec['stopped_reason'] = reason
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
            break
        time.sleep(0.5)
    rec['exit_code'] = proc.wait()
rec['elapsed_seconds'] = round(time.monotonic() - start, 2)
rec['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
rec['dir_bytes_at_finish'] = dir_bytes()
rec['child_maxrss_native_bytes_on_macos'] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
(logs / f'{name}.job.json').write_text(json.dumps(rec, indent=2) + '\n')
print(json.dumps(rec, indent=2))
raise SystemExit(0 if rec['exit_code'] == 0 else 1)
