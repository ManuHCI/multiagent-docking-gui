"""Download a small, pinned official archive with retries and verify its hash."""
import argparse
import hashlib
import io
from pathlib import Path
import tarfile
import time
import urllib.request

HASHES = {
    "linux-64": "e7274528ceb9c20d048a428d6c22d7e02e268f8ffb762c4c365422347c8b8ba2",
    "osx-64": "49c2375bc7553287404f3e348494f3dc7f7987bcd99fe996d2099a654b384b95",
    "osx-arm64": "bd5b8d3f151c3ad2d92f0fb918806543d713c7f67895c020db4041097a7f004d",
    "win-64": "df531329cb3dd59d93f569bb597c715cc8ac239144718b5c530ca325d5d2e42d",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=HASHES, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    url = f"https://github.com/mamba-org/micromamba-releases/releases/download/2.3.3-0/micromamba-{args.platform}.tar.bz2"
    for attempt in range(1, 11):
        try:
            request_url = url + f"?download=1&studio_attempt={attempt}"
            with urllib.request.urlopen(urllib.request.Request(request_url, headers={"User-Agent": "DockingStudio-build"}), timeout=90) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != HASHES[args.platform]:
                raise ValueError("The build tool archive does not match its published SHA256.")
            break
        except (OSError, TimeoutError) as error:
            if attempt == 10:
                raise
            print(f"Download attempt {attempt} failed: {error}. Retrying.", flush=True)
            time.sleep(5)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:bz2") as archive:
        members = [m for m in archive.getmembers() if m.isfile() and Path(m.name).name in {"micromamba", "micromamba.exe"}]
        if len(members) != 1:
            raise ValueError("The verified archive has an unexpected binary layout.")
        args.destination.parent.mkdir(parents=True, exist_ok=True)
        args.destination.write_bytes(archive.extractfile(members[0]).read())
        args.destination.chmod(0o755)
    print("Verified build tool ready:", args.destination, flush=True)


if __name__ == "__main__":
    main()
