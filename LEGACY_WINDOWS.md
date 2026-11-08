# Docking Studio — Experimental Legacy Beta

This is a separate native desktop edition intended for **Windows 7 SP1,
Windows 8 and Windows 8.1, 64-bit**. It keeps the Home, Projects, Reference
and About sections and the five-step docking workflow. Python and the
scientific libraries are included. Docking does not require WebView2,
Anaconda or a separately installed Python.

**Compatibility limit:** this build has not been executed on Windows 7,
Windows 8 or Windows 8.1. Its runtime is selected for older Windows and is
screened for known incompatible required APIs, but that is not an operating
system execution test. It is an experimental beta; test it on the exact
teaching computer before relying on it. The modern Windows 10/11 edition
remains the preferred download for those systems.

## Install

1. Use the Legacy Windows x64 setup file from the GitHub release.
2. Windows 7 needs **Service Pack 1** and the DLL-directory API update
   **KB2533623**, or a later Microsoft update containing it. Setup checks for
   the required API. It does not install or change Windows updates.
3. Run setup and open **In-Silico Docking Studio Legacy Beta**.

The Microsoft VC++ 2019 runtime and original Windows SDK app-local Universal
CRT files are included beside the executable. No system VC runtime is
replaced. A modern 3D-capable browser is needed only for optional HTML pose
viewing; the native interface, preparation, docking and CSV/ZIP exports do
not depend on a browser.

New legacy studies use compact internal folder and file names to stay within
classic Windows path limits. Full project, target and drug names remain in
the interface and exported records. Keep project storage paths short.

Windows XP, Vista and 32-bit Windows are not included. The application and
installer are not publisher-code-signed. SHA256 checksums accompany the
downloads. Local files and the offline example work without internet;
new protein and drug retrieval require internet and working TLS certificates.

## Separate saved projects

The legacy edition must use its own data folder:

    %LOCALAPPDATA%\InSilicoDockingStudioLegacy\Projects

It installs separately from the modern edition and does not move, convert or
delete modern saved work. Export a complete study ZIP before deleting or
transferring a project. Uninstall retains projects outside the app folder.

## Build and verification

Build with official **CPython 3.8.10 x64**, Tcl/Tk 8.6 and the pinned versions
in `requirements-legacy.txt`, in a separate build prefix. The builder expects
the signed VC++ 2019 14.29 runtime files and original Windows SDK 10.0.10240
app-local UCRT files in `.build-tools/legacy-vc-runtime/`.

```text
python build_legacy.py --setup --zip
```

The builder runs packaged dependency checks, scientific preparation/scoring,
real two-target by two-drug docking, native-window checks, and a static PE
architecture/subsystem/API/VC-export screen. It records native-file hashes,
dependency versions and verification limitations in `build-manifest.json`.

Primary compatibility references:

- [Python Windows requirements](https://docs.python.org/3.8/using/windows.html)
- [PyInstaller 5.13.2](https://pypi.org/project/pyinstaller/5.13.2/)
- [Microsoft Universal CRT deployment](https://learn.microsoft.com/en-us/cpp/windows/universal-crt-deployment)
- [Microsoft C++ runtime deployment](https://learn.microsoft.com/en-us/cpp/windows/deployment-in-visual-cpp)

Developed by **Dr Manu Kumar Shetty, Professor, Department of Pharmacology,
Maulana Azad Medical College, New Delhi**, and **Dr Pratibha Nadig, Founder,
Akshaisarathi Healthcare**.

For medical education. Docking scores and poses do not establish potency,
clinical benefit or safety. The older dependency stack no longer receives
the same maintenance as the current edition.
