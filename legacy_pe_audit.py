"""Static Win7 API/architecture screening; not a substitute for an OS trial."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pefile

# These documented APIs first appeared after Windows 7. Dynamic GetProcAddress
# probes can provide a fallback; an eager import cannot and is rejected here.
NEWER_REQUIRED_APIS = {
    "GetSystemTimePreciseAsFileTime": "Windows 8",
    "GetCurrentThreadStackLimits": "Windows 8",
    "GetOverlappedResultEx": "Windows 8",
    "GetProcessInformation": "Windows 8",
    "SetProcessInformation": "Windows 8",
    "CreateFile2": "Windows 8",
    "WaitOnAddress": "Windows 8",
    "WakeByAddressSingle": "Windows 8",
    "WakeByAddressAll": "Windows 8",
    "DiscardVirtualMemory": "Windows 8.1",
    "OfferVirtualMemory": "Windows 8.1",
    "ReclaimVirtualMemory": "Windows 8.1",
    "SetThreadDescription": "Windows 10",
    "GetThreadDescription": "Windows 10",
    "GetSystemCpuSetInformation": "Windows 10",
    "GetTempPath2W": "Windows 10",
    "GetTempPath2A": "Windows 10",
    "VirtualAlloc2": "Windows 10",
    "VirtualAlloc2FromApp": "Windows 10",
    "MapViewOfFile3": "Windows 10",
}


def audit_payload(folder):
    folder = Path(folder)
    records, issues, providers = [], [], {}
    for path in sorted(folder.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".exe", ".dll", ".pyd"}:
            continue
        relative = str(path.relative_to(folder))
        pe = pefile.PE(str(path), fast_load=True)
        pe.parse_data_directories(directories=[pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"],
                                               pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_EXPORT"],
                                               pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_DELAY_IMPORT"]])
        subsystem = (pe.OPTIONAL_HEADER.MajorSubsystemVersion, pe.OPTIONAL_HEADER.MinorSubsystemVersion)
        if pe.FILE_HEADER.Machine != 0x8664:
            issues.append({"file": relative, "problem": "Native file is not Windows x64."})
        # The subsystem version constrains process executables, not the OS
        # requirements of DLL imports. Microsoft's officially redistributable
        # Win7-compatible UCRT API-set forwarding DLLs report subsystem 10.0.
        # DLL compatibility is screened through their eager imported APIs.
        if path.suffix.lower() == ".exe" and subsystem > (6, 1):
            issues.append({"file": relative, "problem": "PE subsystem minimum is newer than Windows 7.",
                           "minimum": list(subsystem)})
        imports = {}
        for descriptor in getattr(pe, "DIRECTORY_ENTRY_IMPORT", []):
            dll = descriptor.dll.decode("ascii", "replace")
            names = [item.name.decode("ascii", "replace") if item.name else "ordinal:" + str(item.ordinal)
                     for item in descriptor.imports]
            imports[dll] = names
            if dll.lower() in {"kernel32.dll", "kernelbase.dll", "api-ms-win-core-synch-l1-2-0.dll"}:
                for name in names:
                    if name in NEWER_REQUIRED_APIS:
                        issues.append({"file": relative, "dll": dll, "api": name,
                                       "requires": NEWER_REQUIRED_APIS[name]})
        exports = {item.name.decode("ascii", "replace") for item in
                   getattr(getattr(pe, "DIRECTORY_ENTRY_EXPORT", None), "symbols", []) if item.name}
        if path.name.lower() not in providers or path.parent == folder:
            providers[path.name.lower()] = (relative, exports)
        records.append({"file": relative, "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "machine": "x64" if pe.FILE_HEADER.Machine == 0x8664 else hex(pe.FILE_HEADER.Machine),
                        "subsystem_version": list(subsystem), "imports": imports})
        pe.close()
    # Check eager imports against bundled VC providers as well as API version
    # names. A C++ DLL built with a newer toolset might need newer exports even
    # when its filename still says MSVCP140.dll.
    for record in records:
        for dll, names in record["imports"].items():
            if dll.lower().startswith(("msvcp140", "vcruntime140", "concrt140")):
                if dll.lower() not in providers:
                    issues.append({"file": record["file"], "problem": "Required VC runtime DLL missing.", "dll": dll})
                    continue
                provider, exports = providers[dll.lower()]
                for name in names:
                    if not name.startswith("ordinal:") and name not in exports:
                        issues.append({"file": record["file"], "dll": dll, "api": name,
                                       "problem": "Required symbol absent from the bundled legacy VC runtime.",
                                       "provider": provider})
    return {"passed": not issues, "method": "Static PE x64/subsystem/newer-API/VC-export screening",
            "limitation": "Does not prove compatibility. Windows 7/8/8.1 execution trials are still required.",
            "old_windows_execution_tested": False, "files_checked": len(records), "issues": issues, "files": records}
