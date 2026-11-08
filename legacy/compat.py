"""Small standard-library backports used only by the Python 3.8 edition.

Importing this module does not change the interpreter. The legacy entry point
calls install_compatibility before loading the shared scientific backend.
"""
from __future__ import annotations

from pathlib import PurePath


def _is_relative_to(path, *other):
    """Python 3.9's lexical containment helper, without resolving symlinks."""
    if not other:
        raise TypeError("is_relative_to() requires at least one path")
    try:
        path.relative_to(type(path)(*other))
    except ValueError:
        return False
    return True


def install_compatibility():
    """Add only APIs absent from Python 3.8; leave modern runtimes intact."""
    if not hasattr(PurePath, "is_relative_to"):
        PurePath.is_relative_to = _is_relative_to

