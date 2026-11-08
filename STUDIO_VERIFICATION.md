# Docking Studio v1.3.0 — 8 November 2026

The current installer is `dist/InSilicoTeachingSuite-1.3.0-Setup.exe`.
Size: 163,963,296 bytes. SHA256:

    593c7f1484499b5126aa1765a054fc9d9c10f566e6831c55c10ac73e8cedac10

The fixed menu is Home, Projects, Reference, About. Home starts a new
drug–target docking project; its five numbered steps are kept within the
project. Projects lists both modern studies and earlier individual projects,
shows disk usage, resumes saved work and deletes only the selected project
after confirmation. Reference supplies the scientific catalogue and original
protein/drug download links. About includes the logo and both complete credits.
The root address and previous workshop address lead to the same Home screen.

## Completed verification

- All 56 source unit tests passed, including eight project-library tests.
  Bounded deletion rejects storage-root/outside paths, links/junctions and
  running work, and preserves other projects. Listing is read-only.
- Source and packaged UI acceptance each passed all 12 checks. The packaged
  test used a 1280 × 860 viewport and completed in 71.8s
  with zero browser errors or server tracebacks.
- Acceptance exercised fixed navigation, complete credits, references,
  older-project reopening/closing/deletion, a saved two-target × two-drug
  matrix, retained current step, real preparation/3D, fresh three-drug Vina
  docking, CSV/ZIP exports, and cancel/confirm deletion preserving another
  calculated project.
- The package passed diagnostics, scientific preparation/scoring and real
  two-target × two-drug docking with Python/Anaconda paths removed.
- The new setup rejected replacement of a running copy without modifying
  payload/version metadata, then installed successfully once it was closed.
  All 3,623 installed payload hashes matched.
  The native window opened, downloads were enabled, and its server stopped
  on close. Isolated uninstall preserved saved user data.
- The actual installed app was updated, reopened and checked against the
  verified payload. All 201 existing saved
  files were preserved. Its live server shows the four sections, logo,
  references and complete developer credits. No existing projects were deleted.

## Storage review

The installed project library contained six saved projects using 22,410,757
bytes (21.4 MB). Reading their summaries and sizes took less than 0.1 seconds.
Stored project files use disk space; molecular data is loaded when opened.
Repeated library scans during input editing were removed. Full size scans run
only in Projects, in a worker thread. No upload staging files were present.

Evidence is under `.verification/studio/`; the distribution's
`build-manifest.json` records source/payload hashes and scientific checks.
Earlier installers and source backups were retained. Tests were performed
on this Windows build computer; a separate participant-PC trial remains
outside this acceptance. Windows 10/11 x64 and Edge WebView2 are required.

## Public distribution

The public Windows setup includes the full AutoDock Vina Apache-2.0 licence
and combined 3Dmol notices. Its final setup was reinstalled in an isolated
folder and passed the in-use upgrade guard, native window/download/shutdown
checks, all 3,623 payload hashes and saved-data preservation. The installed
local copy also received these notices without changing its executable or
projects. The sanitized public source passed all 56 tests on Windows.

Native Mac builds are recorded in the repository's Build macOS applications
workflow, with separate Apple silicon and Intel build manifests and checksums.
The public Mac minimum is macOS 15; binaries are ad-hoc signed and not
Apple-notarized.
