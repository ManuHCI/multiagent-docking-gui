# Installation and startup

## Windows desktop edition

Use Windows 10/11 x64 with Microsoft Edge WebView2 Runtime. Students do not
need Python, Anaconda or individual scientific-package installation.

The current setup is v1.3.0 and installs Docking Studio application v1.3.0:

    dist\InSilicoTeachingSuite-1.3.0-Setup.exe

Close all running copies before updating. Setup refuses an in-use upgrade
and verifies the copied executable before changing version metadata.
Run setup and open **In-Silico Teaching Suite**.
Docking Studio starts at Home with a fixed Home, Projects, Reference, About menu.
A new project follows five practical steps: protein targets,
drug molecules, preparation, docking, and results/export. There are no quizzes
or question forms. **About** contains the logo and both developers;
**Reference** contains citations and original source links.

For portable use, open dist/InSilicoTeachingSuite/InSilicoTeachingSuite.exe
with its complete folder and _internal/ directory intact. Extract a ZIP first.
The earlier v1.1.0 setup is retained separately.

Desktop work defaults to %LOCALAPPDATA%\InSilicoTeachingSuite\Projects.
The adjacent data folder contains settings, example copies and logs.
INSILICO_DATA_DIR selects an isolated data folder. Close the app window to
stop its server. See [DESKTOP.md](DESKTOP.md).

## Source edition

Run run.bat and open <http://127.0.0.1:8111/study>. The Windows launcher
validates Python environments. The previously verified scientific interpreter
is C:\Users\admin\anaconda3\python.exe; the separate default Python 3.14
should not be used. The included Vina executable is 1.2.7.

Use run.bat check for local diagnostics. The older run.bat workshop shortcut
remains compatible because /workshop redirects to /study.

For a custom environment, set its full executable path in Command Prompt:

    set "INSILICO_PYTHON=D:\my environment\python.exe"
    run.bat

This takes precedence over an activated conda environment. An invalid
selection reports a diagnostic. Clear INSILICO_PYTHON to restore discovery.

## Source edition on a new computer

1. Install Miniconda and open Anaconda Prompt in this application folder.
2. Create and activate the environment:

       conda env create -f environment.yml
       conda activate insilico

3. Obtain the official Windows Vina executable from
   <https://github.com/ccsb-scripps/AutoDock-Vina/releases>. Rename it
   vina.exe and place it beside app.py, or add its folder to PATH.
4. Check and start:

       python -m core.selftest --offline
       python app.py

5. Open <http://127.0.0.1:8111/study>.

Linux/macOS source use needs a platform-appropriate Vina executable; the
Windows launcher and Windows desktop EXE are not used there.

## Offline studies and storage

Copy the complete application, including assets/vendor/ and assets/workshop/.
The bundled 3Dmol.js viewer uses local assets and browser WebGL. It needs no
CDN connection or Python py3Dmol package.

Use local PDB/mmCIF files and SMILES/SDF drug inputs to prepare and dock
offline. PubChem lookup and new PDB downloads require internet. The bundled
EGFR example provides a ready input set. build_workshop.py is an authoring
tool, not a class installation step.

In the source edition, settings.json selects the project root. Each batch
study keeps study.json and one project folder per target. Copy the entire
study to transfer it. Existing project folders and workshop notes are retained.

Restart the source server after updating files. INSILICO_PORT selects a
different port; the default interface is local to this computer. Shared
hosting has no login system, so use individual copies for private work.

## Building and acceptance

    python build_desktop.py
    python build_setup.py

Use python build_desktop.py --zip for a complete portable ZIP. The desktop
build includes packaged scientific checks and a file-hash manifest. The setup
output for this revision is the v1.3.0 EXE named above.

Version 1.3.0 passed source, packaged and isolated installer acceptance on
this computer. See STUDIO_VERIFICATION.md and .verification/studio/ for
results. Rehearse installation, normal shortcut creation, a multi-target/
multi-drug run, exports and reopening a saved study on a separate participant
computer before teaching.

build_installer.py remains the older constructor packaging route; use the
desktop commands above for native Docking Studio.
