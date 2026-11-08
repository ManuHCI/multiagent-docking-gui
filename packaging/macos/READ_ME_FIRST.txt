In-Silico Docking Studio for macOS

Choose the arm64 download for an Apple Silicon Mac (M1 or newer), or the
x86_64 download for an Intel Mac. macOS 15 or newer is required.

Open the DMG and drag In-Silico Docking Studio into Applications, then open
the app from Applications. No Python, Anaconda, Vina installation or browser
setup is needed. The system's native WebKit displays the application.

The menu is always Home, Projects, Reference, About. Home starts a new
drug-target docking project. Projects resumes work, shows disk sizes and
deletes selected projects after confirmation. Reference lists sources and
protein/drug downloads. About contains the logo and complete developer credits.

Projects are saved separately in:
    ~/Library/Application Support/InSilicoTeachingSuite/Projects

Keep that folder when updating or removing the app. Moving the app to Trash
does not delete your projects. Export a complete study ZIP before sharing or
deleting a project. Preparation and docking with local structures or the
bundled EGFR example work offline; new RCSB/PubChem downloads need internet.

Signing status is recorded in the adjacent release-manifest.json.
Builds without an Apple Developer ID are ad-hoc signed, not notarized by
Apple. macOS may block the first launch of such a downloaded build. Review
the app's source and published checksum before deciding to allow it in
System Settings > Privacy & Security. No global security settings need to
be disabled. A locally tested build is not a guarantee of Apple notarization.

For medical education. Docking scores and poses do not establish clinical
benefit, potency or safety.

Developed by Dr Manu Kumar Shetty, Professor, Department of Pharmacology,
Maulana Azad Medical College, New Delhi, and Dr Pratibha Nadig, Founder,
Akshaisarathi Healthcare.
