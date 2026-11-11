"""Registered input FITS files: locate them on local disk (find_input) or download them explicitly (fetch).

Nothing is downloaded unless you run ``deimos fetch``. The registration CSV may carry a ``sha256``
column; when it does, a file whose digest differs is refused (find_input) or deleted (fetch).
"""
import hashlib, os, pathlib, urllib.error, urllib.parse, urllib.request
from .util import sha

TIMEOUT = 60                              # seconds per network operation (connect / read), not for the whole file
SCHEMES = ('https', 'http', 'file')
FETCH_HINT = 'run `deimos fetch <batch.json>` to download it, or pass --input-dir / DEIMOS_INPUT_DIR'


class InputNotFound(FileNotFoundError):
    pass


class InputMismatch(ValueError):
    pass


def find_input(row, input_dirs):
    """Return (path, sha256) of the registered file, verifying the registration's sha256 when present."""
    name = row['file_name']
    for folder in input_dirs:
        path = pathlib.Path(folder) / name
        if path.is_file():
            digest = sha(path)
            pinned = row.get('sha256')
            if isinstance(pinned, str) and pinned and digest != pinned:
                raise InputMismatch(f'{name}: sha256 {digest} does not match the registration pin {pinned}')
            return path, digest
    raise InputNotFound(f'{name} not found in input_dirs {[str(d) for d in input_dirs]}; {FETCH_HINT}')


def _pin(row):
    pinned = row.get('sha256')
    return pinned if isinstance(pinned, str) and pinned else None


def _download(url, part, timeout):
    """Stream url into part; return the sha256 of what was written."""
    scheme = urllib.parse.urlparse(url).scheme
    if scheme not in SCHEMES:
        raise ValueError(f'unsupported URL scheme {scheme!r} (expected one of {SCHEMES})')
    h = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=timeout) as response, open(part, 'wb') as out:
        for block in iter(lambda: response.read(1 << 20), b''):
            out.write(block); h.update(block)
        out.flush(); os.fsync(out.fileno())
    return h.hexdigest()


MIRROR = 'https://archive.stsci.edu/hlsps/'


def candidate_urls(url):
    """The registration URL, then (for MAST HLSP download URLs) the archive.stsci.edu/hlsps mirror path; both are sha256-checked."""
    urls = [url]
    if 'uri=mast:HLSP/' in url:
        urls.append(MIRROR + url.split('uri=mast:HLSP/', 1)[1])
    return urls


def fetch_one(row, dest, timeout=TIMEOUT):
    """Download one registered file into dest. Returns 'fetched' or 'skipped'; raises on any failure.

    The bytes go to <name>.part first and are renamed into place only after the sha256 pin (when the registration
    has one) and the optional `bytes` column match; on any failure the .part file is removed, so nothing unverified
    ever sits under its final name. If the first URL fails, the mirror path is tried.
    """
    name = row['file_name']
    if name != pathlib.Path(name).name or not name.endswith('.fits'):
        raise ValueError(f'unsafe file_name {name!r} in the registration (must be a bare *.fits name)')
    pinned = _pin(row); final = dest / name; part = dest / (name + '.part')
    if final.is_file() and (pinned is None or sha(final) == pinned):
        return 'skipped'
    url = row.get('url')
    if not isinstance(url, str) or not url:
        raise ValueError(f'{name}: the registration has no url for this file')
    size = row.get('bytes'); errors = []
    for candidate in candidate_urls(url):
        try:
            digest = _download(candidate, part, timeout)
            if pinned and digest != pinned:
                raise InputMismatch(f'{name}: downloaded sha256 {digest} does not match the registration pin {pinned}; '
                                    'the file was discarded (corrupt transfer, or the archive file changed)')
            if size == size and size not in (None, '') and part.stat().st_size != int(size):
                raise InputMismatch(f'{name}: downloaded size {part.stat().st_size} != registered bytes {size}; discarded')
            os.replace(part, final)
            return 'fetched'
        except urllib.error.URLError as e:
            errors.append(f'{candidate}: {e.reason}')
        except (OSError, ValueError) as e:
            errors.append(f'{candidate}: {e}')
        finally:
            part.unlink(missing_ok=True)
    raise OSError(f'{name}: could not fetch ({"; ".join(errors)}). Check your network connection and retry; finished files are kept and skipped next time')


def fetch(rows, dest, timeout=TIMEOUT, log=None):
    """Fetch every row into dest (created if needed). Returns dict(fetched, skipped, failed) with file names; failed holds
    dict(file, error). A failure on one file does not stop the others."""
    dest = pathlib.Path(dest); dest.mkdir(parents=True, exist_ok=True)
    result = dict(fetched=[], skipped=[], failed=[])
    for i, row in enumerate(rows, 1):
        name = row.get('file_name')
        if log:
            log(f'[{i}/{len(rows)}] {name}')
        try:
            result[fetch_one(row, dest, timeout)].append(name)
        except (OSError, ValueError) as e:
            result['failed'].append(dict(file=name, error=str(e)))
    return result
