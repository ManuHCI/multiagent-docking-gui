# In Silico Docking Studio 1.4.1

This app helps students and teachers prepare proteins and drug molecules,
run molecular docking, inspect structures in 3D and save reports.

Developed by Dr Manu Kumar Shetty, Professor, Department of Pharmacology,
Maulana Azad Medical College, New Delhi, and Dr Pratibha Nadig, Founder,
Akshaisarathi Healthcare.

## Download only one installer

| Your computer | Download |
| --- | --- |
| Windows 10 or 11, 64 bit | [Windows installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.1/InSilicoDockingStudio-1.4.1-Windows-Setup.exe) |
| Mac with an Apple chip, macOS 15 or newer | [Apple Mac installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.1/InSilicoDockingStudio-1.4.1-macOS-arm64.dmg) |
| Mac with an Intel processor, macOS 15 or newer | [Intel Mac installer](https://github.com/ManuHCI/multiagent-docking-gui/releases/download/v1.4.1/InSilicoDockingStudio-1.4.1-macOS-x86_64.dmg) |

On Windows, open the setup file and follow the instructions. Microsoft Edge
WebView2 Runtime is required. On Mac, open the DMG and drag the app into
Applications. Your saved projects are kept when you update.

You do not need the checksum file or source code to install the app.

## A correction, first

Preparing a protein did not remove its water. The step reported the waters
as removed and listed each one under "Components removed", but they stayed
in the prepared structure, were given hydrogen atoms, and appeared in the
3D view and in the downloaded file. On the bundled EGFR structure, 20
waters went in and 60 water atoms came out.

Docking results were not affected. The receptor that AutoDock Vina reads is
built separately and has always had water removed, so scores, poses and
checks from version 1.4.0 stand. What was wrong was the picture, the
downloaded protein file and the report.

The tick box that keeps waters now works as well, so a structural water in
the pocket can be kept on purpose.

## What else is new

**The downloaded code now explains itself.** The Colab notebooks for
docking and for MD have a plain-English explanation before every cell, and
the few lines that actually perform the docking are shown and explained
instead of running out of sight. What a search box is, what exhaustiveness
and the random seed do, how to tell whether the search has settled, and
what the results file contains.

**See the result in 3D.** Each downloaded code package now contains
`complex.pdb`, holding the protein and the docked drug together, with a
ready-made script for ChimeraX. Open `view_in_chimerax.cxc` and the drug
appears in its pocket with the contacting amino acids labelled and the
hydrogen bonds drawn. ChimeraX is free for academic use. A script for the
older UCSF Chimera is included too. Keep the files in one folder.

**Which amino acids the drug is touching.** After docking, a table lists
every residue within 4.5 angstroms of the drug, with the distance and the
kind of contact. The same list is written to `contacts.txt` in the download.
These are distances: a hydrogen bond also needs the right angle, so polar
contacts are called possible, and the ChimeraX script works out the real
ones.

**How to read the score.** A new panel converts the score into the range of
binding strengths it is actually consistent with, and gives the ligand
efficiency for each drug. A ten-fold change in potency is worth only 1.42
kcal/mol, and AutoDock Vina's error is about 2.5, so two scores a few
tenths apart are the same score. The panel says so.

**A glossary** of 24 terms under Reference, each in plain English, each
ending with the part that is most often misread.

## Known limitation in this version

The built-in worked examples were generated before the water correction, so
the offline EGFR example still shows water in its prepared protein. Your
own projects are unaffected. The examples will be rebuilt in the next
release.

## MD and external services

The desktop exports MD code; it does not run a simulation. The soluble
protein protocol requires review of missing regions, protonation and
retained components. Membrane proteins such as SV2A require a separately
built lipid system with appropriate parameters. Membrane mode stops until
those inputs are supplied. A stable simulation alone does not establish
binding.

SwissADME runs on its own website. The app does not generate or invent
SwissADME predictions. Local descriptors are labelled as calculated
properties.

## Installation notes

These applications do not have a publisher signing certificate or Apple
notarisation. Windows may show an unknown publisher message. On Mac, first
launch may need approval in System Settings > Privacy & Security.

Windows 7, 8 and 8.1 users can find the separate experimental legacy edition
in [version 1.3.0](https://github.com/ManuHCI/multiagent-docking-gui/releases/tag/v1.3.0).
That edition has not been verified on those original operating systems.

For medical education and exploratory studies. Docking and ADME predictions
do not establish potency, clinical benefit, dose or safety.
