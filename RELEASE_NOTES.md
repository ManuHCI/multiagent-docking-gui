# In Silico Docking Studio 1.4.0

This app helps students and teachers prepare proteins and drug molecules,
run molecular docking, inspect structures in 3D and save reports.

Developed by Dr Manu Kumar Shetty, Professor, Department of Pharmacology,
Maulana Azad Medical College, New Delhi, and Dr Pratibha Nadig, Founder,
Akshaisarathi Healthcare.

## Download only one installer

| Your computer | Download |
| --- | --- |
| Windows 10 or 11, 64 bit | [Windows installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.0/InSilicoDockingStudio-1.4.0-Windows-Setup.exe) |
| Mac with an Apple chip, macOS 15 or newer | [Apple Mac installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.0/InSilicoDockingStudio-1.4.0-macOS-arm64.dmg) |
| Mac with an Intel processor, macOS 15 or newer | [Intel Mac installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.0/InSilicoDockingStudio-1.4.0-macOS-x86_64.dmg) |

On Windows, open the setup file and follow the instructions. Microsoft Edge
WebView2 Runtime is required. On Mac, open the DMG and drag the app into
Applications. Your saved projects are kept when you update.

You do not need the checksum file or source code to install the app.

## What is new

1. See the protein before and after preparation, next to each other.
2. View and download the prepared drug structure.
3. Read what preparation changed and what each change means.
4. Save a simple docking report or a printable report.
5. Download docking code, actual inputs and a Google Colab notebook for the selected pair.
6. Open SwissADME and copy or download the prepared drug SMILES for its ADME report.
7. Download MD code, chemical structures and a Colab notebook for the selected pair.

The UI colours are retained. Docking now uses the prepared drug coordinates
shown in the viewer. Pose recovery and score comparison are explained
separately. A failed score comparison does not overturn experimental
evidence for an established drug and target.

## MD and external services

The desktop exports MD code; it does not run a simulation. The soluble
protein protocol requires review of missing regions, protonation and retained
components. Membrane proteins such as SV2A require a separately built lipid
system with appropriate parameters. Membrane mode stops until those inputs
are supplied. A stable simulation alone does not establish binding.

SwissADME runs on its own website. The app does not generate or invent
SwissADME predictions. Local descriptors are labelled as calculated properties.

## Installation notes

These applications do not have a publisher signing certificate or Apple
notarisation. Windows may show an unknown publisher message. On Mac, first
launch may need approval in System Settings > Privacy & Security.

Windows 7, 8 and 8.1 users can find the separate experimental legacy edition
in [version 1.3.0](https://github.com/ManuHCI/multiagent-docking-gui/releases/tag/v1.3.0).
That edition has not been verified on those original operating systems.

For medical education and exploratory studies. Docking and ADME predictions
do not establish potency, clinical benefit, dose or safety.
