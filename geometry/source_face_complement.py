"""Compatibility import; maintained implementation is in geometry_tools.algorithms."""
from pathlib import Path as _Path
import sys as _sys
_root = str(_Path(__file__).resolve().parents[1])
if _root not in _sys.path: _sys.path.insert(0, _root)
from geometry_tools._legacy import forward as _forward
_forward(globals(), 'source_face_complement')
