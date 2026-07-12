import hashlib
import json
from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import Any

CACHE_DIR = Path(".harmonizome_cache")
CACHE_DIR.mkdir(exist_ok=True)


def cache_to_file(func: Callable[..., Any]) -> Callable[..., Any]:
    """Persist pure JSON-serializable function results between calls."""

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        key = f"{func.__name__}_{args}_{kwargs}"
        filename = CACHE_DIR / f"{hashlib.md5(key.encode()).hexdigest()}.json"
        if filename.exists():
            with filename.open("r", encoding="utf-8") as file_handle:
                return json.load(file_handle)
        result = func(*args, **kwargs)
        with filename.open("w", encoding="utf-8") as file_handle:
            json.dump(result, file_handle)
        return result

    return wrapper
