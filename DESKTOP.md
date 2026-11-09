# Windows desktop edition: Docking Studio v1.4.0

Docking Studio opens in its own window at Home. Its menu always reads
Home · Projects · Reference · About. Home starts a new docking project;
Projects resumes saved work, shows disk usage and deletes selected projects.
Deletion removes local project files after confirmation. Export a complete ZIP
first to retain a copy. Saved project folders use disk space; they are not
loaded into working memory until opened.
It includes the scientific Python runtime, AutoDock Vina, local 3D viewer,
logo and prepared example inputs. Students do not need Python or Anaconda.
There are no quiz questions or required prediction/discussion forms.

## Install and open

1. Use Windows 10/11 x64 with Microsoft Edge WebView2 Runtime. WebView2 is a
   stated Windows prerequisite rather than a bundled system runtime.
2. Close running copies before updating, then run
   dist/InSilicoDockingStudio-1.4.0-Windows-Setup.exe. The installer refuses an in-use
   update and verifies the copied executable before updating version metadata.
3. Open **In-Silico Teaching Suite** from its shortcut. It starts Docking
   Studio at /study.
4. Select **Start new project**, then follow **Protein targets → Drug molecules → Prepare and inspect →
   Run docking → Results and export**.

Add multiple structures and drugs to run all target–drug combinations.
Inspect each target's chain selection, preparation and binding box before
docking. Each target retains its own results and checks. The offline EGFR
example supplies reference inputs; its docking runs when selected.

**Reference** opens original protein/drug sources, the scientific
catalogue and a downloadable methods/reference record. **About** displays
the logo and both developers. Previous individual projects open from **Projects** with their original
steps. Closing the desktop window stops its local server.

## Portable folder and saved work

Open InSilicoTeachingSuite.exe inside the complete dist/InSilicoTeachingSuite/
folder. Keep _internal/ beside it; the EXE alone is insufficient. Extract an
optional portable ZIP before use. The earlier v1.1.0 setup is retained.

Desktop studies default to:

    %LOCALAPPDATA%\InSilicoTeachingSuite\Projects

A study contains study.json and separate project folders for its targets.
These retain prepared structures, poses, checks, versions and settings.
The CSV preserves target–drug identity; the complete ZIP adds structures,
methods, checks, source links and references. Copy the whole study to
transfer it. Existing individual projects and old workshop.json notes remain;
the old /workshop route redirects to the practical workspace.

Settings, writable example copies and logs are kept outside the application
folder. Check logs/desktop.log, logs/runtime.log or logs/startup-error.log in
the desktop data folder for diagnostics. Uninstall is designed to retain
separately stored studies. INSILICO_DATA_DIR selects an isolated data folder
without moving existing work.

## Offline use and interpretation

Local PDB/mmCIF files, supplied SMILES/SDF drugs and bundled example inputs
can be prepared and docked offline. New structure retrieval, PubChem lookup
and external services require internet. WebGL supports the 3D viewer.

A completed calculation can have a failed scientific check. Scores do not
establish potency, clinical benefit or safety. Raw scores across different
proteins do not establish selectivity. Prepared structures, simple reports,
SwissADME input and selected-pair docking and MD code exports are available
in the main study workflow. Membrane MD requires a separately prepared lipid system.

## Build and verify

Use a Windows Python environment containing the working scientific packages,
PyInstaller and pywebview. Project-local tools may live in .build-tools/.

    python build_desktop.py
    python build_setup.py

The first command bundles and checks dist/InSilicoTeachingSuite/; the second
creates the v1.4.0 setup EXE. Use python build_desktop.py --zip for a complete
portable ZIP. build-manifest.json records hashes, versions and packaged
checks. NSIS 3.11 or newer is needed; build_setup.py can verify and extract
the official nsis-3.11.zip from .build-tools/.

The desktop build runs packaged diagnostics, scientific calculations and a
real two-target by two-drug docking check. `verify_setup.py` exercises an
isolated setup without desktop shortcut integration; `verify_studio_ui.py`
checks the practical interface, credits, source links and downloads.
Evidence is recorded under `.verification/studio/` and in the distribution's
`build-manifest.json`. A separate participant-PC installation, normal
shortcut check and class rehearsal remain necessary.
