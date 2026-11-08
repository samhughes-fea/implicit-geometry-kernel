"""ctypes binding to ``libgeos_c``. WKB in, WKB out.

Locate the library with ``GEOS_LIBRARY_PATH``, then conda prefixes, then the
system loader (``find_library("geos_c")``). A later C++ kernel links the same C API.
"""

from __future__ import annotations

import atexit
import ctypes
import ctypes.util
import os
import sys
from pathlib import Path

_EMPTY_POLYGON_WKB = bytes([1, 3, 0, 0, 0, 0, 0, 0, 0])

GEOS_GEOMETRYCOLLECTION = 7
GEOSBUF_CAP_ROUND = 1
GEOSBUF_JOIN_ROUND = 1
GEOSBUF_JOIN_MITRE = 2

_GEOS_C = None
_HANDLE = None
_ERR: list[str] = []

_Notice = ctypes.CFUNCTYPE(None, ctypes.c_char_p)


@_Notice
def _on_notice(_msg: bytes | None) -> None:
    return


@_Notice
def _on_error(msg: bytes | None) -> None:
    _ERR.append(msg.decode("utf-8", "replace") if msg else "GEOS error")


def _candidate_paths() -> list[Path]:
    out: list[Path] = []
    env = os.environ.get("GEOS_LIBRARY_PATH")
    if env:
        out.append(Path(env))
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        root = Path(prefix)
        out.extend(
            [
                root / "Library" / "bin" / "geos_c.dll",
                root / "bin" / "geos_c.dll",
                root / "lib" / "libgeos_c.so",
                root / "lib" / "libgeos_c.so.1",
                root / "lib" / "libgeos_c.dylib",
            ]
        )
    found = ctypes.util.find_library("geos_c")
    if found:
        out.append(Path(found))
    return out


def _load_lib():
    last: Exception | None = None
    for path in _candidate_paths():
        if not path.exists() and path.suffix != "":
            continue
        try:
            if sys.platform == "win32" and path.parent.is_dir():
                os.add_dll_directory(str(path.parent))
            return ctypes.CDLL(str(path))
        except OSError as exc:
            last = exc
            continue
    found = ctypes.util.find_library("geos_c")
    if found:
        try:
            return ctypes.CDLL(found)
        except OSError as exc:
            last = exc
    raise OSError(
        "Could not load libgeos_c. Install GEOS (conda-forge geos, or libgeos-dev) "
        f"or set GEOS_LIBRARY_PATH to geos_c.dll / libgeos_c.so. Last error: {last}"
    )


def _bind(lib) -> None:
    p = ctypes.c_void_p
    lib.GEOS_init_r.restype = p
    lib.GEOS_finish_r.argtypes = [p]
    lib.GEOSContext_setNoticeHandler_r.argtypes = [p, _Notice]
    lib.GEOSContext_setErrorHandler_r.argtypes = [p, _Notice]
    lib.GEOSGeomFromWKB_buf_r.argtypes = [p, ctypes.c_char_p, ctypes.c_size_t]
    lib.GEOSGeomFromWKB_buf_r.restype = p
    lib.GEOSGeomToWKB_buf_r.argtypes = [p, p, ctypes.POINTER(ctypes.c_size_t)]
    lib.GEOSGeomToWKB_buf_r.restype = ctypes.POINTER(ctypes.c_ubyte)
    lib.GEOSFree_r.argtypes = [p, p]
    lib.GEOSGeom_destroy_r.argtypes = [p, p]
    lib.GEOSBufferParams_create_r.argtypes = [p]
    lib.GEOSBufferParams_create_r.restype = p
    lib.GEOSBufferParams_destroy_r.argtypes = [p, p]
    lib.GEOSBufferParams_setEndCapStyle_r.argtypes = [p, p, ctypes.c_int]
    lib.GEOSBufferParams_setJoinStyle_r.argtypes = [p, p, ctypes.c_int]
    lib.GEOSBufferParams_setMitreLimit_r.argtypes = [p, p, ctypes.c_double]
    lib.GEOSBufferParams_setQuadrantSegments_r.argtypes = [p, p, ctypes.c_int]
    lib.GEOSBufferWithParams_r.argtypes = [p, p, p, ctypes.c_double]
    lib.GEOSBufferWithParams_r.restype = p
    lib.GEOSIntersection_r.argtypes = [p, p, p]
    lib.GEOSIntersection_r.restype = p
    lib.GEOSDifference_r.argtypes = [p, p, p]
    lib.GEOSDifference_r.restype = p
    lib.GEOSUnaryUnion_r.argtypes = [p, p]
    lib.GEOSUnaryUnion_r.restype = p
    lib.GEOSMakeValid_r.argtypes = [p, p]
    lib.GEOSMakeValid_r.restype = p
    lib.GEOSGeom_createCollection_r.argtypes = [p, ctypes.c_int, ctypes.POINTER(p), ctypes.c_uint]
    lib.GEOSGeom_createCollection_r.restype = p
    lib.GEOSisValidReason_r.argtypes = [p, p]
    lib.GEOSisValidReason_r.restype = ctypes.c_void_p
    lib.GEOSisValid_r.argtypes = [p, p]
    lib.GEOSisValid_r.restype = ctypes.c_char


def _handle():
    global _GEOS_C, _HANDLE
    if _HANDLE is not None:
        return _HANDLE
    lib = _load_lib()
    _bind(lib)
    _GEOS_C = lib
    h = lib.GEOS_init_r()
    if not h:
        raise OSError("GEOS_init_r failed")
    lib.GEOSContext_setNoticeHandler_r(h, _on_notice)
    lib.GEOSContext_setErrorHandler_r(h, _on_error)
    _HANDLE = h
    atexit.register(_finish)
    return h


def _finish() -> None:
    global _HANDLE
    if _HANDLE is None or _GEOS_C is None:
        return
    _GEOS_C.GEOS_finish_r(_HANDLE)
    _HANDLE = None


def _from_wkb(h, wkb: bytes):
    if not wkb:
        wkb = _EMPTY_POLYGON_WKB
    g = _GEOS_C.GEOSGeomFromWKB_buf_r(h, wkb, len(wkb))
    if not g:
        raise RuntimeError(_ERR[-1] if _ERR else "GEOSGeomFromWKB_buf_r failed")
    return g


def _to_wkb(h, geom) -> bytes:
    if not geom:
        return _EMPTY_POLYGON_WKB
    size = ctypes.c_size_t()
    buf = _GEOS_C.GEOSGeomToWKB_buf_r(h, geom, ctypes.byref(size))
    if not buf:
        return _EMPTY_POLYGON_WKB
    try:
        return ctypes.string_at(buf, size.value)
    finally:
        _GEOS_C.GEOSFree_r(h, buf)


def _destroy(h, geom) -> None:
    if geom:
        _GEOS_C.GEOSGeom_destroy_r(h, geom)


def buffer(
    wkb: bytes,
    distance: float,
    *,
    quad_segs: int = 8,
    cap_style: int = GEOSBUF_CAP_ROUND,
    join_style: int = GEOSBUF_JOIN_ROUND,
    mitre_limit: float = 5.0,
) -> bytes:
    h = _handle()
    geom = _from_wkb(h, wkb)
    params = _GEOS_C.GEOSBufferParams_create_r(h)
    out = None
    try:
        _GEOS_C.GEOSBufferParams_setEndCapStyle_r(h, params, int(cap_style))
        _GEOS_C.GEOSBufferParams_setJoinStyle_r(h, params, int(join_style))
        _GEOS_C.GEOSBufferParams_setMitreLimit_r(h, params, float(mitre_limit))
        _GEOS_C.GEOSBufferParams_setQuadrantSegments_r(h, params, int(quad_segs))
        out = _GEOS_C.GEOSBufferWithParams_r(h, geom, params, float(distance))
        return _to_wkb(h, out)
    finally:
        if params:
            _GEOS_C.GEOSBufferParams_destroy_r(h, params)
        _destroy(h, geom)
        _destroy(h, out)


def _binary(op_name: str, a: bytes, b: bytes) -> bytes:
    h = _handle()
    ga = _from_wkb(h, a)
    gb = _from_wkb(h, b)
    out = None
    try:
        out = getattr(_GEOS_C, op_name)(h, ga, gb)
        return _to_wkb(h, out)
    finally:
        _destroy(h, ga)
        _destroy(h, gb)
        _destroy(h, out)


def intersection(a: bytes, b: bytes) -> bytes:
    return _binary("GEOSIntersection_r", a, b)


def difference(a: bytes, b: bytes) -> bytes:
    return _binary("GEOSDifference_r", a, b)


def unary_union(wkbs: bytes | list[bytes]) -> bytes:
    h = _handle()
    if isinstance(wkbs, (bytes, bytearray)):
        geom = _from_wkb(h, bytes(wkbs))
        out = None
        try:
            out = _GEOS_C.GEOSUnaryUnion_r(h, geom)
            return _to_wkb(h, out)
        finally:
            _destroy(h, geom)
            _destroy(h, out)
    items = [bytes(w) for w in wkbs]
    if not items:
        return _EMPTY_POLYGON_WKB
    geoms = [_from_wkb(h, w) for w in items]
    arr = (ctypes.c_void_p * len(geoms))(*geoms)
    coll = _GEOS_C.GEOSGeom_createCollection_r(h, GEOS_GEOMETRYCOLLECTION, arr, len(geoms))
    out = None
    try:
        if not coll:
            raise RuntimeError(_ERR[-1] if _ERR else "GEOSGeom_createCollection_r failed")
        out = _GEOS_C.GEOSUnaryUnion_r(h, coll)
        return _to_wkb(h, out)
    finally:
        _destroy(h, coll)
        _destroy(h, out)


def make_valid(wkb: bytes) -> bytes:
    h = _handle()
    geom = _from_wkb(h, wkb)
    out = None
    try:
        out = _GEOS_C.GEOSMakeValid_r(h, geom)
        return _to_wkb(h, out)
    finally:
        _destroy(h, geom)
        _destroy(h, out)


def is_valid(wkb: bytes) -> bool:
    h = _handle()
    geom = _from_wkb(h, wkb)
    try:
        return bool(_GEOS_C.GEOSisValid_r(h, geom))
    finally:
        _destroy(h, geom)


def is_valid_reason(wkb: bytes) -> str:
    h = _handle()
    geom = _from_wkb(h, wkb)
    try:
        ptr = _GEOS_C.GEOSisValidReason_r(h, geom)
        if not ptr:
            return "Unknown"
        try:
            return ctypes.cast(ptr, ctypes.c_char_p).value.decode("utf-8", "replace")
        finally:
            _GEOS_C.GEOSFree_r(h, ptr)
    finally:
        _destroy(h, geom)


__all__ = [
    "GEOSBUF_CAP_ROUND",
    "GEOSBUF_JOIN_MITRE",
    "GEOSBUF_JOIN_ROUND",
    "buffer",
    "difference",
    "intersection",
    "is_valid",
    "is_valid_reason",
    "make_valid",
    "unary_union",
]
