"""Response sanitization for external API data.

Prevents toxic data flow attacks by cleaning API responses
before returning them to LLM context.
"""

from __future__ import annotations

import re

# Max field length to prevent context overflow
MAX_FIELD_LENGTH = 10_000
MAX_LIST_ITEMS = 500

# Patterns that could be prompt injection attempts
_INJECTION_PATTERNS = [
    re.compile(r"<\s*(?:system|user|assistant)\s*>", re.IGNORECASE),
    re.compile(r"\b(?:ignore|forget|disregard)\s+(?:all\s+)?(?:previous|above|prior)\b", re.IGNORECASE),
    re.compile(r"\byou\s+are\s+now\b", re.IGNORECASE),
    re.compile(r"\b(?:act|pretend|behave)\s+as\b", re.IGNORECASE),
]


def sanitize_string(value: str, max_length: int = MAX_FIELD_LENGTH) -> str:
    """Truncate and flag suspicious strings from external APIs."""
    if not isinstance(value, str):
        return str(value)[:max_length]
    truncated = value[:max_length]
    for pattern in _INJECTION_PATTERNS:
        if pattern.search(truncated):
            truncated = pattern.sub("[FILTERED]", truncated)
    return truncated


def sanitize_dict(data: dict, max_length: int = MAX_FIELD_LENGTH) -> dict:
    """Recursively sanitize a dict from external API response."""
    result = {}
    for key, value in data.items():
        skey = sanitize_string(str(key), 256)
        if isinstance(value, str):
            result[skey] = sanitize_string(value, max_length)
        elif isinstance(value, dict):
            result[skey] = sanitize_dict(value, max_length)
        elif isinstance(value, list):
            result[skey] = sanitize_list(value, max_length)
        else:
            result[skey] = value
    return result


def sanitize_list(data: list, max_length: int = MAX_FIELD_LENGTH) -> list:
    """Sanitize a list from external API response, truncating if too long."""
    truncated = data[:MAX_LIST_ITEMS]
    result = []
    for item in truncated:
        if isinstance(item, str):
            result.append(sanitize_string(item, max_length))
        elif isinstance(item, dict):
            result.append(sanitize_dict(item, max_length))
        elif isinstance(item, list):
            result.append(sanitize_list(item, max_length))
        else:
            result.append(item)
    return result
