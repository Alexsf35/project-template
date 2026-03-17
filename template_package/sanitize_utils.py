# sanitize_utils.py
from __future__ import annotations

from typing import Dict, Any, Iterable, Optional
import math


# ----------------------------
# Missing / NA detection
# ----------------------------

def _is_na(v: Any) -> bool:
    """
    Treat None, empty strings, 'nan'/'na'/'null'/'none' (any case),
    and float('nan') as missing.
    """
    if v is None:
        return True
    # float NaN
    if isinstance(v, float):
        try:
            if math.isnan(v):
                return True
        except Exception:
            pass
    # string-ish empties
    if isinstance(v, str) and v.strip().lower() in {"", "nan", "na", "null", "none"}:
        return True
    return False


# ----------------------------
# Primitive coercers
# ----------------------------

def to_int_or_none(v: Any) -> Optional[int]:
    if _is_na(v):
        return None
    try:
        if isinstance(v, float):
            return int(v) if v.is_integer() else int(round(v))
        if isinstance(v, str):
            s = v.strip()
            if "." in s:
                f = float(s)
                return int(f) if f.is_integer() else int(round(f))
            return int(s)
        return int(v)
    except Exception:
        return None


def to_float_or_none(v: Any) -> Optional[float]:
    if _is_na(v):
        return None
    try:
        f = float(v)
        return None if (isinstance(f, float) and math.isnan(f)) else f
    except Exception:
        return None


def to_bool_or_none(v: Any) -> Optional[bool]:
    if _is_na(v):
        return None
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"true", "t", "yes", "y", "1"}:
        return True
    if s in {"false", "f", "no", "n", "0"}:
        return False
    return None


# ----------------------------
# CSV-safe string
# ----------------------------

def safe_str(s: Any) -> str:
    """
    Make strings safe for Neo4j CSV import:
      - remove CR/LF/TAB
      - normalize curly quotes/primes
      - remove all double quotes (")
      - replace ASCII apostrophe (') with prime (′)
      - trim surrounding spaces
    """
    if s is None:
        return ""
    s = str(s)
    s = s.replace("\r", " ").replace("\n", " ").replace("\t", " ").strip()
    s = (
        s.replace("’", "'")
         .replace("‘", "'")
         .replace("“", '"')
         .replace("”", '"')
         .replace("′", "'")
         .replace("″", '"')
    )
    s = s.replace('"', "")
    s = s.replace("'", "′")  # U+2032
    return s


# ----------------------------
# Dict-level utilities
# ----------------------------

def normalize_types(
    props: Dict[str, Any],
    int_fields: Iterable[str] = (),
    float_fields: Iterable[str] = (),
    bool_fields: Iterable[str] = (),
    *,
    in_place: bool = True,
) -> Dict[str, Any]:
    """
    Coerce selected keys to int/float/bool (for Neo4j :long/:double/:boolean).
    """
    target = props if in_place else dict(props)
    for k in int_fields:
        if k in target:
            target[k] = to_int_or_none(target[k])
    for k in float_fields:
        if k in target:
            target[k] = to_float_or_none(target[k])
    for k in bool_fields:
        if k in target:
            target[k] = to_bool_or_none(target[k])
    return target


def sanitize_strings(
    props: Dict[str, Any],
    *,
    exclude: Iterable[str] = (),
    in_place: bool = True,
) -> Dict[str, Any]:
    """
    Apply safe_str() to all values except keys in `exclude`
    (use exclude to skip already-typed fields).
    """
    target = props if in_place else dict(props)
    exclude_set = set(exclude)
    for k, v in list(target.items()):
        if k in exclude_set:
            continue
        target[k] = safe_str(v)
    return target


def lowercase_bool_values(
    props: Dict[str, Any],
    *,
    in_place: bool = True,
) -> Dict[str, Any]:
    """
    Convert Python bools to Neo4j-friendly lowercase literals.
    Avoids importing "True"/"False" strings which Neo4j misreads.
    """
    target = props if in_place else dict(props)
    for k, v in list(target.items()):
        if isinstance(v, bool):
            target[k] = "true" if v else "false"
    return target


def coerce_by_type_map(
    props: Dict[str, Any],
    type_map: Dict[str, str],
    *,
    in_place: bool = True,
) -> Dict[str, Any]:
    """
    One-shot coercion if you prefer a single mapping:
      type_map = {"left_end_position":"int","centisome_position":"float","interrupted":"bool"}
    Unknown types fall back to string sanitization.
    """
    target = props if in_place else dict(props)
    for k, t in type_map.items():
        if k not in target:
            continue
        v = target[k]
        if t == "int":
            target[k] = to_int_or_none(v)
        elif t == "float":
            target[k] = to_float_or_none(v)
        elif t == "bool":
            target[k] = to_bool_or_none(v)
        else:
            target[k] = safe_str(v)
    return target


def sanitize_row(
    props: Dict[str, Any],
    int_fields: Iterable[str] = (),
    float_fields: Iterable[str] = (),
    bool_fields: Iterable[str] = (),
    *,
    in_place: bool = True,
) -> Dict[str, Any]:
    """
    Convenience: normalize selected types, then safe-string everything else.
    """
    target = props if in_place else dict(props)
    normalize_types(target, int_fields, float_fields, bool_fields, in_place=True)
    sanitize_strings(target, exclude=set(int_fields) | set(float_fields) | set(bool_fields), in_place=True)
    lowercase_bool_values(target, in_place=True)
    return target
