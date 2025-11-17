"""Mixin providing dict-like access to pyOpenMS meta values."""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Iterable

import pyopenms as oms


_MISSING = object()


class MetaInfoMappingMixin:
    """Expose :class:`pyopenms.MetaInfoInterface` objects as mappings."""

    _META_KEY_TYPES = (str, bytes)

    def _meta_object(self) -> oms.MetaInfoInterface:
        """Return the native MetaInfoInterface implementation."""
        raise NotImplementedError

    # ----------------------------- Helpers -----------------------------

    @staticmethod
    def _coerce_key(key: Any) -> str:
        if isinstance(key, str):
            return key
        if isinstance(key, bytes):
            return key.decode("utf-8")
        raise TypeError(
            "Meta value keys must be str or bytes, "
            f"got {type(key).__name__}"
        )

    def _meta(self) -> oms.MetaInfoInterface:
        return self._meta_object()

    # ---------------------------- Mapping API ----------------------------

    def __getitem__(self, key: Any) -> Any:
        normalized = self._coerce_key(key)
        if not self._meta().metaValueExists(normalized):
            raise KeyError(normalized)
        return self._meta().getMetaValue(normalized)

    def __setitem__(self, key: Any, value: Any) -> None:
        normalized = self._coerce_key(key)
        self._meta().setMetaValue(normalized, value)

    def __delitem__(self, key: Any) -> None:
        normalized = self._coerce_key(key)
        if not self._meta().metaValueExists(normalized):
            raise KeyError(normalized)
        self._meta().removeMetaValue(normalized)

    def __contains__(self, key: object) -> bool:  # pragma: no cover - trivial
        try:
            normalized = self._coerce_key(key)
        except TypeError:
            return False
        return bool(self._meta().metaValueExists(normalized))

    def get(self, key: Any, default: Any = None) -> Any:
        normalized = self._coerce_key(key)
        if self._meta().metaValueExists(normalized):
            return self._meta().getMetaValue(normalized)
        return default

    def setdefault(self, key: Any, default: Any = None) -> Any:
        normalized = self._coerce_key(key)
        if not self._meta().metaValueExists(normalized):
            self._meta().setMetaValue(normalized, default)
            return default
        return self._meta().getMetaValue(normalized)

    def pop(self, key: Any, default: Any = _MISSING) -> Any:
        normalized = self._coerce_key(key)
        if not self._meta().metaValueExists(normalized):
            if default is _MISSING:
                raise KeyError(normalized)
            return default
        value = self._meta().getMetaValue(normalized)
        self._meta().removeMetaValue(normalized)
        return value

    def update(self, other: Any = None, **kwargs: Any) -> None:
        if other is not None:
            items: Iterable[tuple[Any, Any]]
            if isinstance(other, Mapping):
                items = other.items()
            else:
                items = other
            for key, value in items:
                self[key] = value
        for key, value in kwargs.items():
            self[key] = value

    def keys(self):  # pragma: no cover - simple proxy
        return list(self._meta().getMetaValues().keys())

    def items(self):  # pragma: no cover - simple proxy
        return list(self._meta().getMetaValues().items())

    def values(self):  # pragma: no cover - simple proxy
        return list(self._meta().getMetaValues().values())

    def clear(self) -> None:  # pragma: no cover - simple proxy
        self._meta().clearMetaInfo()
