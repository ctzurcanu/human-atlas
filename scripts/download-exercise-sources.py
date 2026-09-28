#!/usr/bin/env python3
"""Download missing public exercise snapshots and verify the pinned content hashes."""
import hashlib
import json
from pathlib import Path
import time
import urllib.request

manifest = json.loads(Path('scripts/exercise-sources.json').read_text())
for source in manifest['files']:
    path = Path(source['path'])
    if not path.exists():
        if source['key'].startswith('aaos-'):
            request = urllib.request.Request(source['url'], headers={'User-Agent': 'HumanAtlas/1.0 source snapshot'})
            with urllib.request.urlopen(request, timeout=60) as response:
                content = response.read()
        else:
            url = source['url']
            result = None
            while url:
                request = urllib.request.Request(url, headers={'User-Agent': 'HumanAtlas/1.0 source snapshot'})
                for attempt in range(4):
                    try:
                        with urllib.request.urlopen(request, timeout=60) as response:
                            page = json.load(response)
                        break
                    except Exception:
                        if attempt == 3:
                            raise
                        time.sleep(2 ** attempt)
                if result is None:
                    result = {**page, 'results': list(page['results'])}
                else:
                    result['results'].extend(page['results'])
                url = page['next']
            assert len(result['results']) == result['count'], 'Incomplete source pagination'
            result.update(next=None, previous=None)
            content = json.dumps(result).encode()
        assert hashlib.sha256(content).hexdigest() == source['sha256'], f"Live source differs from the pinned snapshot: {source['key']}"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source['sha256'], f"Snapshot changed: {path}"
    print(source['key'], 'verified')
