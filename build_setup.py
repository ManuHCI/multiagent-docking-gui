"""Wrap the verified standalone bundle in one per-user Windows setup EXE."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
NSIS_HASH = "c7d27f780ddb6cffb4730138cd1591e841f4b7edb155856901cdf5f214394fa1"


def compiler() -> Path:
    local = ROOT / ".build-tools" / "nsis-3.11" / "makensis.exe"
    if local.is_file():
        return local
    archive = ROOT / ".build-tools" / "nsis-3.11.zip"
    if archive.is_file():
        if hashlib.sha256(archive.read_bytes()).hexdigest() != NSIS_HASH:
            raise RuntimeError("NSIS archive does not match the official 3.11 SHA256.")
        with zipfile.ZipFile(archive) as bundle:
            base = (ROOT / ".build-tools").resolve()
            for member in bundle.infolist():
                target = (base / member.filename).resolve()
                if not target.is_relative_to(base):
                    raise RuntimeError("Unsafe path in NSIS archive.")
            bundle.extractall(base)
        if local.is_file():
            return local
    found = shutil.which("makensis")
    if found:
        return Path(found)
    raise RuntimeError("Install NSIS 3.11+ or put the official nsis-3.11.zip in .build-tools. See DESKTOP.md.")


def quote_nsis(value: str) -> str:
    return value.replace("$", "$$").replace('"', '$\\"')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-unverified", action="store_true",
                        help="Allow an intermediate build lacking a successful diagnostic manifest.")
    args = parser.parse_args()
    try:
        import json
        payload = ROOT / "dist" / "InSilicoTeachingSuite"
        manifest_path = payload / "build-manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not args.allow_unverified and not (manifest.get("verification") or {}).get("can_start"):
            raise RuntimeError("Run python build_desktop.py and pass packaged diagnostics first.")
        for relative, record in manifest["files"].items():
            path = (payload / relative).resolve()
            if not path.is_relative_to(payload.resolve()) or not path.is_file():
                raise RuntimeError(f"Missing or invalid payload file: {relative}")
            if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
                raise RuntimeError(f"Payload changed since verification: {relative}")
        includes = ROOT / "build" / "desktop" / "uninstall-files.nsh"
        includes.parent.mkdir(parents=True, exist_ok=True)
        files = sorted(p for p in payload.rglob("*") if p.is_file())
        lines = ['Delete "$INSTDIR\\' + quote_nsis(str(p.relative_to(payload))) + '"' for p in files]
        directories = sorted((p for p in payload.rglob("*") if p.is_dir()),
                             key=lambda p: len(p.parts), reverse=True)
        lines.extend('RMDir "$INSTDIR\\' + quote_nsis(str(p.relative_to(payload))) + '"' for p in directories)
        includes.write_text("\n".join(lines) + "\n", encoding="utf-8-sig")
        output = ROOT / "dist" / "InSilicoTeachingSuite-1.3.0-Setup.exe"
        executable_hash = hashlib.sha256((payload / "InSilicoTeachingSuite.exe").read_bytes()).hexdigest()
        subprocess.run([str(compiler()), "/V3", "/INPUTCHARSET", "UTF8", f"/DPAYLOAD={payload}",
                        f"/DEXE_SHA256={executable_hash}",
                        f"/DOUTPUT={output}", f"/DUNINSTALL_FILES={includes}",
                        str(ROOT / "installer" / "desktop.nsi")],
                       cwd=ROOT / "installer", check=True)
        digest = hashlib.sha256(output.read_bytes()).hexdigest()
        output.with_suffix(".exe.sha256").write_text(f"{digest}  {output.name}\n", encoding="ascii")
        print(f"Setup executable: {output}\nSize: {output.stat().st_size / 1e6:.1f} MB\nSHA256: {digest}")
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"Setup build failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
