"""Wrapper for :class:`pyopenms.Feature` exposing dict-like meta access."""
from __future__ import annotations

from typing import Optional

import pyopenms as oms

from ._meta_mapping import MetaInfoMappingMixin


class Py_Feature(MetaInfoMappingMixin):
    """A light-weight wrapper around :class:`pyopenms.Feature`."""

    def __init__(self, native_feature: Optional[oms.Feature] = None):
        self._feature = native_feature if native_feature is not None else oms.Feature()

    # ----------------------------- Meta access -----------------------------

    def _meta_object(self) -> oms.MetaInfoInterface:
        return self._feature

    # ----------------------------- Delegation -----------------------------

    @property
    def native(self) -> oms.Feature:
        """Return the underlying :class:`pyopenms.Feature`."""

        return self._feature

    def __getattr__(self, name: str):  # pragma: no cover - simple delegation
        return getattr(self._feature, name)

    def __setattr__(self, name: str, value):  # pragma: no cover - simple delegation
        if name == "_feature":
            object.__setattr__(self, name, value)
            return
        setattr(self._feature, name, value)

    def __repr__(self) -> str:  # pragma: no cover - debug helper
        uid = getattr(self._feature, "getUniqueId", lambda: None)()
        rt = getattr(self._feature, "getRT", lambda: None)()
        mz = getattr(self._feature, "getMZ", lambda: None)()
        return f"Py_Feature(uid={uid}, rt={rt}, mz={mz})"
