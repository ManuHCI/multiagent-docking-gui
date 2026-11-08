# Legacy Beta verification — 8 November 2026

Version: 1.3.0-legacy.1. Intended systems: Windows 7 SP1, Windows 8 and
Windows 8.1, 64-bit. **Execution on these original operating systems has
not been verified.** This is an experimental compatibility build.

The build uses official Python 3.8.10, native Tk controls, RDKit 2022.9.5,
NumPy 1.22.4, SciPy 1.7.3, Pandas 1.3.5, OpenMM 8.1.1, PDBFixer 1.9 and
Meeko 0.5.1. It has no NiceGUI, WebView2, .NET or browser dependency for its
interface, preparation or docking. Optional HTML 3D viewing uses a browser.

## Completed checks

- All 61 source regression tests passed on the current Windows host.
- Python 3.8 native UI acceptance passed nine checks: four sections, logo,
  complete credits, references, saved two-by-two inputs, native 2D molecules,
  real preparation and three-drug docking, CSV/ZIP/3D exports, reopening,
  confirmed deletion preserving other work and older-project record viewing.
- The frozen EXE passed dependency checks, all five scientific checks, real
  two-target by two-drug docking and native-window startup checks.
- Static PE screening checked 408 native files with no identified unsupported
  required APIs, architecture errors or missing VC-runtime exports. This is
  not an exhaustive old-OS execution test.
- All 52 app-local VC++ 2019/UCRT DLL providers were copied from authenticated
  Microsoft redistributables and matched their recorded hashes.
- The final installer was installed in a deep path containing spaces. All
  2,508 payload hashes matched after fresh installation and upgrade. The
  in-use update guard returned 10 without changing files; the stopped retry
  succeeded. Installed scientific, four-pair docking and native UI checks passed.
- An initial path-limit failure was reproduced at 261 characters. Compact
  internal names fixed it; the same deep-folder test then passed while full
  project and drug labels remained in records.
- Silent uninstall removed the QA application and preserved the real saved
  project and notes files byte-for-byte. Modern installation and user data
  were untouched.
- The installer checksum parser passed eight compiled cases covering older
  spaced/tabbed/uppercase certutil output and malformed-digest rejection.

Execution testing used **Windows 11, NT build 26200**. No Windows 7, 8 or 8.1
machine or licensed VM was available for a target-OS trial. Test the beta on
the exact teaching computers before relying on it. Windows 7 also requires
SP1 and the DLL-directory update KB2533623 or a later replacement.

The release includes an installer, portable ZIP, source ZIP, SHA256 files
and a compatibility JSON report. The existing modern release tag and
Windows/macOS installers remain unchanged. Legacy source is tagged
`v1.3.0-legacy.1`.
