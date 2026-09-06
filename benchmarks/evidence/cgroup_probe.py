#!/usr/bin/env python3
"""Check native evidence failure contracts inside real, bounded Docker cgroups."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import uuid


CONTROLLER = r'''set -u
mkdir -p /output
for key in memory.max memory.swap.max memory.events memory.peak; do
  cat "/sys/fs/cgroup/$key" > "/output/$key.before" 2>/dev/null || true
done
if [ ! -s /output/memory.max.before ]; then
  echo 'cgroup-v2 memory.max is unavailable' >&2
  exit 79
fi
"$@" > /output/analysis.stdout 2> /output/analysis.stderr
code=$?
printf '%s\n' "$code" > /output/analysis.exit_code
for key in memory.max memory.swap.max memory.events memory.peak; do
  cat "/sys/fs/cgroup/$key" > "/output/$key.after" 2>/dev/null || true
done
for receipt in /output/result.tsv.manifest.json /output/result.tsv.manifest.json.partial; do
  if [ -f "$receipt" ]; then
    /inputs/rosalind verify --manifest "$receipt" --json > "$receipt.verify.json" 2> "$receipt.verify.stderr"
    printf '%s\n' "$?" > "$receipt.verify.exit_code"
  fi
done
exit "$code"
'''


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path, limit=32 << 20):
    with path.open('rb') as stream:
        raw = stream.read(limit + 1)
    if len(raw) > limit:
        raise ValueError(f'metadata exceeds {limit} bytes: {path}')
    return json.loads(raw)


def events(path):
    if not path.is_file():
        return None
    result = {}
    for line in path.read_text().splitlines():
        key, value = line.split()
        result[key] = int(value)
    return result


def number(path):
    if not path.is_file():
        return None
    value = path.read_text().strip()
    return None if value == 'max' else int(value)


def check_case(case, result):
    """Do not infer an OOM from exit137 alone or a partial from its filename."""
    errors = []
    code = result['analysis_exit_code']
    before, after = result['events_before'], result['events_after']
    event_oom = (before is not None and after is not None
                 and after.get('oom_kill', 0) > before.get('oom_kill', 0))
    oom = result['docker_oom_killed'] or event_oom
    success = result['success_artifact']
    receipt = result['success_receipt']
    partial = result['partial_artifact']
    partial_receipt = result['partial_receipt']
    if result['timed_out']:
        errors.append('container exceeded the bounded timeout')
    if result['memory_max_bytes'] != case['hard_limit_bytes']:
        errors.append('observed memory.max does not match the requested hard limit')
    if result['swap_max_bytes'] != 0:
        errors.append('memory.swap.max was not zero')
    if code != 0 and (success or receipt):
        errors.append('failed native invocation left a successful destination')
    if case['expected'] == 'completed':
        if code != 0 or oom or not success or not receipt:
            errors.append('admitted native run did not complete with an artifact and receipt')
        if result.get('receipt_status') != 'completed' or not result.get('verified'):
            errors.append('completed artifact did not pass independent native verification')
        if case.get('require_os_limit') and result.get('assurance') != 'cgroup-v2':
            errors.append('required OS assurance was not recorded')
    elif case['expected'] == 'refused':
        if code != 3 or oom or partial or partial_receipt:
            errors.append('expected a preflight refusal without artifacts, partials, or OOM')
    elif case['expected'] == 'startup-breach':
        if code != 4 or oom or partial or partial_receipt:
            errors.append('expected an observed startup-budget breach before output creation')
    elif case['expected'] == 'capacity':
        if code != 4 or oom or not partial or not partial_receipt:
            errors.append('expected declared-capacity failure with explicit partial evidence')
        if result.get('partial_status') in (None, 'completed'):
            errors.append('partial receipt does not identify a failed run')
        if not result.get('partial_integrity'):
            errors.append('partial receipt integrity was not independently verified')
    elif case['expected'] == 'oom':
        if not oom or code in (None, 0):
            errors.append('kernel OOM was not established by memory.events or Docker OOMKilled')
    else:
        raise ValueError('unknown expected outcome')
    return errors


class Docker:
    def __init__(self, context, raw):
        self.prefix = ['docker', '--context', context]
        self.raw = raw
        self.sequence = 0

    def call(self, args, *, check=True, timeout=120):
        self.sequence += 1
        stem = self.raw / f'docker-{self.sequence:03}'
        argv = self.prefix + list(map(str, args))
        stem.with_suffix('.argv.json').write_text(json.dumps(argv, indent=2) + '\n')
        with stem.with_suffix('.stdout').open('wb') as out, stem.with_suffix('.stderr').open('wb') as err:
            proc = subprocess.run(argv, stdout=out, stderr=err, timeout=timeout, check=False)
        stdout = stem.with_suffix('.stdout').read_text(errors='replace')
        if check and proc.returncode:
            raise RuntimeError(f'Docker command failed ({proc.returncode}); see {stem}.stderr')
        return proc.returncode, stdout


def collect(root, case, inspect, timed_out):
    state = inspect['State']
    success = root / 'result.tsv'
    manifest = root / 'result.tsv.manifest.json'
    partial = root / 'result.tsv.partial'
    partial_manifest = root / 'result.tsv.manifest.json.partial'
    result = {
        'case': case,
        'analysis_exit_code': number(root / 'analysis.exit_code'),
        'container_exit_code': state['ExitCode'],
        'docker_oom_killed': state['OOMKilled'],
        'timed_out': timed_out,
        'memory_max_bytes': number(root / 'memory.max.before'),
        'swap_max_bytes': number(root / 'memory.swap.max.before'),
        'memory_peak_bytes': number(root / 'memory.peak.after'),
        'events_before': events(root / 'memory.events.before'),
        'events_after': events(root / 'memory.events.after'),
        'success_artifact': success.is_file(),
        'success_receipt': False,
        'partial_artifact': partial.is_file(),
        'partial_receipt': False,
        'files': {str(p.relative_to(root)): {'bytes': p.stat().st_size, 'sha256': sha256(p)}
                  for p in sorted(root.rglob('*')) if p.is_file()},
    }
    # If the kernel killed the controller too, retain that absence explicitly;
    # Docker's container exit still establishes a failed invocation.
    if result['analysis_exit_code'] is None:
        result['analysis_exit_code'] = state['ExitCode']
        result['analysis_exit_source'] = 'container-only; controller did not record child exit'
    else:
        result['analysis_exit_source'] = 'controller'
    for path in [manifest, partial_manifest]:
        if path.is_file():
            m = read_json(path)
            params = m.get('params', {})
            # The first-party CLI keeps the requested receipt path for identified
            # resource failures; the managed SDK suffixes that path with .partial.
            # Neither convention makes a failed receipt a completed artifact.
            completed = params.get('run_status') == 'completed'
            prefix = '' if path == manifest and completed else 'partial_'
            result['success_receipt' if prefix == '' else 'partial_receipt'] = True
            result[prefix + 'receipt_path'] = path.name
            result[prefix + 'claim_hash'] = params.get('manifest_blake3')
            result[prefix + 'producer'] = {
                key: value for key, value in params.items()
                if key.startswith(('producer.', 'code_', 'rustc_', 'target_', 'deps_'))}
            result[prefix + 'measurements'] = m.get('measurements', {})
            result['partial_status' if prefix else 'receipt_status'] = params.get('run_status')
            result[prefix + 'assurance'] = params.get('contract.assurance')
            verified_path = Path(str(path) + '.verify.json')
            if verified_path.is_file():
                v = read_json(verified_path)
                result[prefix + 'verified'] = v.get('ok') is True
                result[prefix + 'integrity'] = (
                    v.get('trust', {}).get('receipt_integrity', {}).get('state') == 'satisfied')
    result['errors'] = check_case(case, result)
    result['passed'] = not result['errors']
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--docker-context', required=True)
    parser.add_argument('--image', required=True, help='existing Linux/amd64 image with sh, awk, and the binary runtime libraries')
    parser.add_argument('--binary', required=True, type=Path, help='Linux/amd64 ELF Rosalind binary')
    parser.add_argument('--reference', required=True, type=Path, help='uncompressed indexed FASTA')
    parser.add_argument('--alignments', required=True, type=Path)
    parser.add_argument('--alignment-index', type=Path)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--sites', type=Path)
    selection.add_argument('--regions', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--memory-mib', type=int, default=128)
    parser.add_argument('--native-oom-mib', type=int, default=0,
                        help='optionally require a native OOM at this intentionally tiny hard limit (minimum6MiB)')
    parser.add_argument('--timeout', type=int, default=180)
    parser.add_argument('--fields', default='depths,alleles')
    parser.add_argument('--sample', help='select a declared read-group sample')
    args = parser.parse_args(argv)
    if not 64 <= args.memory_mib <= 512 or args.timeout <= 0:
        parser.error('admitted memory must be64–512MiB and timeout must be positive')
    if args.native_oom_mib and not 6 <= args.native_oom_mib < 32:
        parser.error('native OOM hard limit must be6–31MiB')
    if args.reference.suffix == '.gz':
        parser.error('use an uncompressed indexed FASTA for this bounded probe')
    with args.binary.open('rb') as stream:
        header = stream.read(20)
    if header[:4] != b'\x7fELF' or header[18:20] != b'>\x00':
        parser.error('--binary must be a Linux/amd64 ELF executable')
    cram = args.alignments.suffix == '.cram'
    suffixes = ['.crai'] if cram else ['.bai', '.csi']
    index = args.alignment_index
    if index is None:
        index = next((p for s in suffixes for p in [Path(str(args.alignments) + s), args.alignments.with_suffix(s)] if p.is_file()), None)
    if index is None:
        parser.error('alignment index missing; supply --alignment-index')
    selected = args.sites or args.regions
    selection_name = ('selection.bcf' if selected.suffix == '.bcf' else
                      'selection.vcf.gz' if str(selected).endswith('.gz') else
                      'selection.vcf') if args.sites else 'selection.bed'
    files = {'rosalind': args.binary, 'reference.fa': args.reference,
             'reference.fa.fai': Path(str(args.reference) + '.fai'),
             'alignments.cram' if cram else 'alignments.bam': args.alignments,
             'alignment.index': index, selection_name: selected}
    for path in files.values():
        if not path.is_file():
            parser.error(f'missing input: {path}')
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    (output / 'harness.py').write_bytes(Path(__file__).read_bytes())
    raw = output / 'raw'
    raw.mkdir()
    docker = Docker(args.docker_context, raw)
    report = {'schema': 1, 'kind': 'rosalind-cgroup-evidence-probe', 'status': 'running',
              'source_hashes': {str(p.resolve()): {'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in files.values()},
              'harness_sha256': sha256(output / 'harness.py'), 'docker_context': args.docker_context,
              'cases': [], 'limitations': ['The allocation OOM control is not a Rosalind workload.',
              'Kernel SIGKILL cannot run cooperative cleanup or guarantee a partial receipt.',
              'Docker cgroup memory charges include page cache and the small controller process; they are not process RSS.',
              'This is a failure-contract probe, not a performance benchmark.']}
    save = lambda: (output / 'report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    save()
    identity = uuid.uuid4().hex[:12]
    volume, carrier = f'rosalind-cgroup-{identity}', f'rosalind-cgroup-stage-{identity}'
    created = []
    volume_created = False
    try:
        _, daemon_json = docker.call(['info', '--format',
            '{"kernel":{{json .KernelVersion}},"architecture":{{json .Architecture}},"os":{{json .OperatingSystem}},"cgroup_version":{{json .CgroupVersion}},"memory_bytes":{{.MemTotal}},"cpus":{{.NCPU}}}'])
        report['docker_daemon'] = json.loads(daemon_json)
        _, image_json = docker.call(['image', 'inspect', '--platform', 'linux/amd64', args.image])
        image = json.loads(image_json)[0]
        report['image'] = {key: image.get(key) for key in ['Id', 'RepoDigests', 'Architecture', 'Os']}
        if image.get('Architecture') != 'amd64' or image.get('Os') != 'linux':
            raise ValueError('the supplied image must resolve to Linux/amd64')
        run_image = (image.get('RepoDigests') or [image['Id']])[0]
        report['image']['execution_reference'] = run_image
        docker.call(['volume', 'create', volume]); volume_created = True
        docker.call(['create', '--name', carrier, '--platform', 'linux/amd64', '--network', 'none',
                     '--memory', '128m', '--memory-swap', '128m', '--cpus', '2', '--pids-limit', '32',
                     '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges',
                     '--mount', f'type=volume,src={volume},dst=/inputs', run_image,
                     '/bin/sh', '-c', 'sha256sum /inputs/*'])
        created.append(carrier)
        for name, path in files.items():
            docker.call(['cp', path.resolve(), f'{carrier}:/inputs/{name}'])
        _, staged = docker.call(['start', '--attach', carrier], timeout=args.timeout)
        staged_hashes = {line.split(None, 1)[1].strip(): line.split(None, 1)[0]
                         for line in staged.splitlines()}
        report['staged_source_sha256'] = staged_hashes
        for name, path in files.items():
            expected = report['source_hashes'][str(path.resolve())]['sha256']
            if staged_hashes.get('/inputs/' + name) != expected:
                raise ValueError(f'staged input changed during transfer: {name}')
        common = ['/inputs/rosalind', 'analyze', 'evidence', '--reference', '/inputs/reference.fa',
                  '--alignments', '/inputs/alignments.cram' if cram else '/inputs/alignments.bam',
                  '--alignment-index', '/inputs/alignment.index', '--sites' if args.sites else '--regions',
                  '/inputs/' + selection_name, '--fields', args.fields, '--format', 'tsv',
                  '--output', '/output/result.tsv', '--manifest', '/output/result.tsv.manifest.json', '--enforce']
        if args.sample:
            common += ['--sample', args.sample]
        cases = [
            {'name': 'admitted-cooperative', 'expected': 'completed', 'budget_mib': args.memory_mib},
            {'name': 'admitted-os-limit', 'expected': 'completed', 'budget_mib': args.memory_mib, 'require_os_limit': True},
            {'name': 'preflight-refusal', 'expected': 'refused', 'budget_mib': args.memory_mib,
             'extra': ['--max-record-bytes', str(args.memory_mib << 21)]},
            {'name': 'startup-budget-breach', 'expected': 'startup-breach', 'budget_mib': 1},
            {'name': 'declared-capacity', 'expected': 'capacity', 'budget_mib': args.memory_mib, 'extra': ['--max-read-len', '1']},
            {'name': 'allocation-oom-control', 'expected': 'oom', 'hard_mib': 32, 'control': True},
        ]
        if args.native_oom_mib:
            cases.append({'name': 'native-kernel-oom', 'expected': 'oom', 'hard_mib': args.native_oom_mib, 'budget_mib': args.memory_mib})
        for case in cases:
            hard = case.get('hard_mib', args.memory_mib)
            case['hard_limit_bytes'] = hard << 20
            label = case['name']
            name = f'rosalind-cgroup-{identity}-{label}'
            destination = output / label
            destination.mkdir()
            command = (['awk', 'BEGIN { s=sprintf("%1024s", ""); for(i=0;i<262144;i++) a[i]=s i; print length(a) }']
                       if case.get('control') else common + ['--memory-budget-mb', str(case['budget_mib'])]
                       + (['--require-os-limit'] if case.get('require_os_limit') else []) + case.get('extra', []))
            case['argv'] = command
            docker.call(['create', '--name', name, '--platform', 'linux/amd64', '--network', 'none',
                         '--memory', f'{hard}m', '--memory-swap', f'{hard}m', '--cpus', '2', '--pids-limit', '32',
                         '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges',
                         '--mount', f'type=volume,src={volume},dst=/inputs,readonly',
                         run_image, '/bin/sh', '-c', CONTROLLER, 'rosalind-cgroup-controller', *command])
            created.append(name)
            timed_out = False
            started = time.monotonic()
            try:
                docker.call(['start', '--attach', name], check=False, timeout=args.timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                docker.call(['kill', name], check=False)
            _, inspection = docker.call(['inspect', name])
            (destination / 'docker-inspect.json').write_text(inspection)
            docker.call(['cp', f'{name}:/output/.', destination], check=False)
            item = collect(destination, case, json.loads(inspection)[0], timed_out)
            item['wall_seconds'] = time.monotonic() - started
            report['cases'].append(item)
            save()
            docker.call(['rm', name]); created.remove(name)
        completed = [x for x in report['cases'] if x['case']['expected'] == 'completed']
        report['completed_output_equality'] = (len(completed) == 2 and all(x['success_artifact'] for x in completed)
            and len({x['files']['result.tsv']['sha256'] for x in completed}) == 1)
        report['status'] = 'passed' if all(x['passed'] for x in report['cases']) and report['completed_output_equality'] else 'failed'
    except Exception as error:
        report['status'] = 'failed'
        report['harness_error'] = f'{type(error).__name__}: {error}'
    finally:
        cleanup_errors = []
        for name in reversed(created):
            try:
                code, _ = docker.call(['rm', '--force', name], check=False)
                if code:
                    cleanup_errors.append('container ' + name)
            except Exception as error:
                cleanup_errors.append(f'container {name}: {error}')
        if volume_created:
            try:
                code, _ = docker.call(['volume', 'rm', volume], check=False)
                if code:
                    cleanup_errors.append('volume ' + volume)
            except Exception as error:
                cleanup_errors.append(f'volume {volume}: {error}')
        if cleanup_errors:
            report['cleanup_errors'] = cleanup_errors
            report['status'] = 'failed'
        save()
    print(json.dumps({'status': report['status'], 'report': str(output / 'report.json')}))
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    sys.exit(main())
