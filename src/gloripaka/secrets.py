"""Secure credential access — the TabPFN API token, never hard-coded.

The token lives OUTSIDE the repo at ``~/.config/priorlabs/tabpfn_token`` (mode 600),
or in the ``TABPFN_TOKEN`` environment variable. Source code and the git tree never
contain it; ``.gitignore`` additionally blocks any ``*_token`` / ``*api_key*`` file
from being committed by accident.
"""

from __future__ import annotations

import os
from pathlib import Path

TOKEN_FILE = Path.home() / ".config" / "priorlabs" / "tabpfn_token"


def _clean(raw: str) -> str:
    """Normalize a pasted token: drop markdown-escape backslashes, quotes, whitespace.

    Prior Labs keys look like ``tabpfn_sk_<...>``; when pasted from a markdown preview
    the underscores often arrive escaped as ``tabpfn\\_sk\\_...``, which silently breaks
    auth with a 401. Stripping backslashes makes ingestion robust to that.
    """
    return raw.replace("\\", "").strip().strip('"').strip("'").strip()


def get_tabpfn_token() -> str | None:
    """Return the TabPFN API token from env or the secure file, else ``None``.

    Precedence: ``TABPFN_TOKEN`` env var first (useful for CI/one-off overrides),
    then the on-disk secure file. The value is normalized (see :func:`_clean`).
    """
    env = os.environ.get("TABPFN_TOKEN")
    if env and env.strip():
        return _clean(env)
    if TOKEN_FILE.exists():
        text = _clean(TOKEN_FILE.read_text())
        return text or None
    return None


def masked_token() -> str:
    """A safe-to-log fingerprint of the token (never the token itself)."""
    tok = get_tabpfn_token()
    if not tok:
        return "<none>"
    return f"{tok[:3]}…{tok[-3:]} (len {len(tok)})"
