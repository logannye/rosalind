#!/usr/bin/env python3
"""Allow publication retries only when existing wheel filenames have exact bytes."""
import argparse
import email.parser
import hashlib
import json
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


def fetch_json(url):
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            return json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise


def verify(directory: Path, index: str, fetch=fetch_json):
    host = {'pypi': 'https://pypi.org', 'testpypi': 'https://test.pypi.org'}[index]
    wheels = sorted(directory.glob('*.whl'))
    if not wheels:
        raise ValueError('no wheels to publish')
    reports = []
    for wheel in wheels:
        with zipfile.ZipFile(wheel) as archive:
            metadata_paths = [name for name in archive.namelist() if name.endswith('.dist-info/METADATA')]
            if len(metadata_paths) != 1:
                raise ValueError(f'{wheel.name}: ambiguous wheel metadata')
            metadata = email.parser.Parser().parsestr(archive.read(metadata_paths[0]).decode())
        if metadata['Name'] != 'rosalind-bio':
            raise ValueError(f'{wheel.name}: unexpected package name')
        remote = fetch(f"{host}/pypi/rosalind-bio/{metadata['Version']}/json")
        digest = hashlib.sha256()
        with wheel.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        matching = [entry for entry in (remote or {}).get('urls', []) if entry.get('filename') == wheel.name]
        if matching and (len(matching) != 1 or matching[0].get('digests', {}).get('sha256') != digest.hexdigest()):
            raise ValueError(f'{wheel.name}: existing index artifact has different bytes; refuse overwrite/skip')
        reports.append({'filename': wheel.name, 'sha256': digest.hexdigest(), 'existing_exact': bool(matching)})
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('dist'))
    parser.add_argument('--index', choices=['pypi', 'testpypi'], required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.directory, args.index), indent=2))


if __name__ == '__main__':
    main()
