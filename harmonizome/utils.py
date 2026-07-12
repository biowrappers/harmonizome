import hashlib
import json
from functools import wraps
from pathlib import Path

CACHE_DIR = Path(".harmonizome_cache")
CACHE_DIR.mkdir(exist_ok=True)


def cache_to_file(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        key = f"{func.__name__}_{args}_{kwargs}"
        filename = CACHE_DIR / f"{hashlib.md5(key.encode()).hexdigest()}.json"
        if filename.exists():
            with open(filename, "r") as f:
                return json.load(f)
        result = func(*args, **kwargs)
        with open(filename, "w") as f:
            json.dump(result, f)
        return result

    return wrapper
