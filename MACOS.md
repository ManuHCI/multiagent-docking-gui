# Install on a Mac

Download the Mac installer from the [public release page](https://github.com/ManuHCI/multiagent-docking-gui/releases/tag/v1.4.0).

Choose **macOS-arm64.dmg** for an Apple silicon Mac (M1 or newer), or
**macOS-x86_64.dmg** for an Intel Mac. macOS 15 or newer is required.
Open the DMG, drag **In-Silico Docking Studio** into **Applications**, then
open it from Applications. Python, Anaconda and AutoDock Vina are bundled.

The app uses native macOS WebKit. Its menu is **Home · Projects · Reference · About**.
The docking steps, logo, complete developer credits and source links match
the Windows edition. Local files and the example work offline; new protein
or drug downloads require internet.

## First launch

These public builds are ad-hoc signed and are not notarized by Apple.
macOS may block a downloaded app on its first launch. The release includes
SHA256 checksums. After checking the download and deciding to allow it,
use **System Settings → Privacy & Security → Open Anyway** if macOS offers
that option. See [Apple's instructions](https://support.apple.com/en-us/102445).

## Saved projects

Projects are saved separately from the application:

    ~/Library/Application Support/InSilicoTeachingSuite/Projects

Updating or removing the application keeps these projects. The Projects
section shows their sizes and provides a confirmation before deletion.
Export a complete study ZIP to keep a copy or move work to another computer.
Windows project paths are rebased when a complete project folder is opened
on a Mac.

## Build from source

Build each architecture on a matching Mac. From a terminal with conda:

```sh
conda env create -f environment-macos.yml
conda activate insilico-macos
python desktop.py
python build_macos.py --arch arm64
# On an Intel Mac: python build_macos.py --arch x86_64
```

The GitHub **Build macOS applications** workflow uses native Apple silicon
and Intel runners. Each build runs source tests, packaged dependency checks,
scientific calculations, a real two-target by two-drug docking study and a
native window startup/shutdown check before creating DMG and ZIP downloads.
Final file hashes and verification results accompany each build.

An Apple Developer ID and Apple notarization are not included in ad-hoc
builds. Adding those requires the publisher's Apple signing credentials.
