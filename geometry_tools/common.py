"""Stable artifact identities, atomic output and local configuration."""
from __future__ import annotations
from pathlib import Path
import hashlib, json, os, tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get('ATLAS_GEOMETRY_CACHE', ROOT / '.local-models/geometry-workbench')).resolve()
VERSION = 'geometry-workbench/0.1.0'


def digest_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def json_bytes(value) -> bytes:
    def default(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,Path):return str(v)
        if hasattr(v,'isoformat'):return v.isoformat()
        raise TypeError(type(v).__name__)
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False,default=default) + '\n').encode()


def identity(*values) -> str:
    h = hashlib.sha256()
    for value in values:
        if isinstance(value, np.ndarray):
            a = np.ascontiguousarray(value)
            h.update(json_bytes({'shape': a.shape, 'dtype': a.dtype.str})); h.update(a.tobytes())
        else:
            h.update(json_bytes(value))
        h.update(b'\0')
    return h.hexdigest()


def atomic_json(path, value, *, overwrite=False):
    return atomic_bytes(path, json_bytes(value), overwrite=overwrite)


def atomic_bytes(path, data, *, overwrite=False):
    path = Path(path).resolve(); path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        if path.read_bytes() == data:
            return path
        raise FileExistsError(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.'+path.name, delete=False) as stream:
        tmp = Path(stream.name)
        try:
            stream.write(data); stream.flush(); os.fsync(stream.fileno())
        except BaseException:
            tmp.unlink(missing_ok=True); raise
    try:
        if overwrite:
            os.replace(tmp, path)
        else:
            # Atomic no-clobber publication, including concurrent writers.
            try:
                os.link(tmp, path)
            except FileExistsError:
                if path.read_bytes() != data:
                    raise
    finally:
        tmp.unlink(missing_ok=True)
    return path


def atomic_npz(path, **arrays):
    path = Path(path).resolve(); path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix='.npz', delete=False) as stream:
        tmp = Path(stream.name)
        try:
            np.savez_compressed(stream, **arrays)
            stream.flush(); os.fsync(stream.fileno())
        except BaseException:
            tmp.unlink(missing_ok=True); raise
    try:
        os.link(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    return path


def owned_path(path, root=CACHE):
    path, root = Path(path).resolve(), Path(root).resolve()
    if not path.is_relative_to(root):
        raise ValueError('Artifact lies outside the declared workbench directory')
    return path
