# In-Silico Docking Studio

Molecular docking software for medical education, with a clear five-step
workflow and local AutoDock Vina calculations. No quiz questions or written
answers are required. An LLM or Ollama installation is not needed.

## Download and install

**[Download version 1.4.0](https://github.com/ManuHCI/multiagent-docking-gui/releases/tag/v1.4.0)**

Download only one installer for your computer. You do not need all files in Assets.

| Computer | Installer |
| --- | --- |
| Windows 10 or 11, 64 bit | [Download Windows setup](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.0/InSilicoDockingStudio-1.4.0-Windows-Setup.exe) |
| Mac with an Apple chip, macOS 15 or newer | [Download Apple Mac installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.0/InSilicoDockingStudio-1.4.0-macOS-arm64.dmg) |
| Intel Mac, macOS 15 or newer | [Download Intel Mac installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.0/InSilicoDockingStudio-1.4.0-macOS-x86_64.dmg) |

The installers include Python and the scientific tools. On Windows, run
setup and open the installed shortcut; [Microsoft Edge WebView2 Runtime](https://developer.microsoft.com/en-us/microsoft-edge/webview2/) is
required. On Mac, open the DMG and drag the app into Applications. The checksum
file is optional and is not needed for installation.

The public binaries are not signed with a publisher certificate or notarized
by Apple. Windows may display a publisher prompt. On Mac, first launch may
require approval in System Settings > Privacy & Security. See
[Mac installation](MACOS.md) and [Windows installation](DESKTOP.md).

![Home screen](docs/home.png)

## Older Windows — experimental legacy beta

A separate native edition is available for **Windows 7 SP1, Windows 8 and
Windows 8.1, 64-bit**:

**[Download the legacy beta installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.3.0/InSilicoDockingStudio-1.3.0-Legacy-Windows-x64-Setup.exe)**

It includes Python 3.8 and older scientific libraries and does not require
WebView2. It installs separately from the modern edition. Windows 7 requires
SP1 and the DLL-directory update KB2533623 or a later update containing it.

**This is an experimental compatibility build, not a verified Windows 7/8/8.1
release.** Its scientific workflow and installer were tested on the current
Windows host, and the native libraries were screened for known incompatible
APIs. Actual execution on the original operating systems remains unverified.
See [legacy installation and limitations](LEGACY_WINDOWS.md). The maintained
Windows 10/11 edition remains the preferred choice on those systems.

## Four clear sections

- **Home:** Start a new drug and target docking project or try the offline example.
- **Projects:** Reopen saved work, see disk sizes, or delete a selected project after confirmation.
- **Reference:** Scientific references, protein/drug download links and methods export.
- **About:** Logo and full developer credits.

Inside a project, follow **Targets → Drugs → Prepare → Dock → Results**.
Add multiple proteins and drugs to calculate every target–drug combination.
Each target keeps its own preparation, binding box, results and scientific
checks. Review poses in 3D and export a CSV or complete project ZIP with
structures, methods, parameters, sources and references.

## Structure views, reports and research code

Prepare and inspect shows the protein before and after preparation, the
prepared drug, added atoms, removed waters, retained components and unresolved
gaps. Use Show atoms to inspect details. The drug coordinates shown here are
the coordinates passed to docking. The original and prepared files can be downloaded.

Results provides a simple text report and a printable HTML report. Pose
recovery and score comparison are explained separately. A failed comparison
does not overturn experimental evidence for an established drug and target.

Select a protein and drug to export its docking inputs, actual search settings,
Python code and Colab notebook. The selected-pair rerun keeps the original
comparison checks in the report; it does not rerun those controls. Differences
between platforms and processor counts may affect scores.

SwissADME is available through a link with prepared SMILES to copy or download.
The app does not automatically obtain SwissADME predictions. Its local chemical
descriptors are not measured pharmacokinetics.

MD export supplies the selected docked drug with chemical bond orders, the
protein, a starting OpenMM protocol and a Colab notebook. No simulation is run
inside the desktop app. Soluble proteins can use the water protocol after
review of chemistry and missing regions. Membrane proteins such as SV2A need
a separately prepared and parameterised lipid system. The membrane runner
stops until that system is provided. Review the README inside each package.

Local PDB/mmCIF files, SMILES/SDF molecules and the bundled EGFR example work
offline. New RCSB/PubChem retrieval requires internet. A completed calculation
can still have a failed scientific check; the outcome stays visible.

## Saved work

Projects are separate from the installed application and are retained on
update or uninstall:

- Windows: `%LOCALAPPDATA%\InSilicoTeachingSuite\Projects`
- Mac: `~/Library/Application Support/InSilicoTeachingSuite/Projects`

Saved projects occupy disk space; their molecular data is loaded when opened.
Export a complete ZIP before deleting a project or moving it between machines.
Earlier individual projects remain available from Projects. Your own projects
and settings are excluded from this public source repository.

## Developers

**Dr Manu Kumar Shetty**, Professor, Department of Pharmacology,
Maulana Azad Medical College, New Delhi.

**Dr Pratibha Nadig**, Founder, Akshaisarathi Healthcare.

## Source and builds

```sh
git clone https://github.com/ManuHCI/multiagent-docking-gui.git
cd multiagent-docking-gui
conda env create -f environment.yml
conda activate insilico
python app.py
```

For Mac source use `environment-macos.yml` instead. See [MACOS.md](MACOS.md).
The Windows source environment also needs the official bundled `vina.exe`.
For a desktop source window, install pywebview and run `python desktop.py`.
Use [DESKTOP.md](DESKTOP.md) for Windows packaging. GitHub's **Build macOS
applications** workflow builds and verifies Apple silicon and Intel apps.

The previous local-LLM application remains in the
[archived branch](https://github.com/ManuHCI/multiagent-docking-gui/tree/archive/pre-docking-studio-20261108)
and its original version tags. This repository now opens Docking Studio.

## Interpretation and licence

For medical education and exploratory in-silico studies. Docking scores and
poses do not establish potency, selectivity, clinical benefit or safety.
Compare scores within each target's setup and review the recorded checks.

Project source is MIT licensed; scientific dependencies retain their own
licences. Full Vina and 3Dmol notices accompany the bundled assets. Use
**Reference** and the exported methods to cite the tools and data actually used.
