"""Small import bridge for frozen scripts using the former geometry directory."""
from importlib import import_module
import sys


def forward(namespace, name):
    """Expose the same canonical objects, including private helpers and CLIs."""
    legacy_name = namespace['__name__']
    implementation = import_module('geometry_tools.algorithms.' + name)
    if legacy_name == '__main__' and callable(getattr(implementation, 'main', None)):
        raise SystemExit(implementation.main())
    metadata = {'__name__', '__package__', '__loader__', '__spec__', '__file__', '__cached__', '__path__'}
    # Explicit file loaders keep the originally created module object, so expose
    # attributes there as well as aliasing ordinary imports in sys.modules.
    namespace.update({k: v for k, v in vars(implementation).items() if k not in metadata})
    sys.modules[legacy_name] = implementation
    return implementation
