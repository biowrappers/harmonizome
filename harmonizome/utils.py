import hashlib
import json
import os
from functools import wraps

CACHE_DIR = ".harmonizome_cache"
os.makedirs(CACHE_DIR, exist_ok=True)


def cache_to_file(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        key = f"{func.__name__}_{args}_{kwargs}"
        filename = os.path.join(
            CACHE_DIR, hashlib.md5(key.encode()).hexdigest() + ".json"
        )
        if os.path.exists(filename):
            with open(filename, "r") as f:
                return json.load(f)
        result = func(*args, **kwargs)
        with open(filename, "w") as f:
            json.dump(result, f)
        return result

    return wrapper
