# Setup Guide

Full step-by-step setup for the Multi-Agent Molecular Docking GUI, from a fresh Windows machine. For Linux/Mac, the steps are the same but substitute platform-specific download paths.

Tested on Windows 11, Python 3.12 and 3.13.

---

## Step 1: Install Python (if not already installed)

Download Python 3.10 or newer from https://www.python.org/downloads/

**Important:** during install, **tick "Add Python to PATH"** on the first screen. Without this, the `py` launcher won't work.

Verify after install (open a new command prompt):

```bash
py -3 --version
```

Should print something like `Python 3.13.1`.

---

## Step 2: Install Ollama

1. Download from https://ollama.com/download
2. Run the installer. No admin rights needed on Windows.
3. Ollama starts as a background service automatically on port 11434.

Verify in a new command prompt:

```bash
ollama --version
ollama list
```

`ollama list` will be empty initially -- that's expected.

---

## Step 3: Pull a language model

Start with a small fast model:

```bash
ollama pull llama3.2:3b
```

This downloads ~2 GB. One-time only.

Sanity check:

```bash
ollama run llama3.2:3b
```

Type `hello`, press Enter, confirm you get a reply, then `/bye`.

**Alternative models** (pick based on RAM):

| RAM | Model | Pull command |
|---|---|---|
| 4-8 GB | llama3.2:3b (default) | `ollama pull llama3.2:3b` |
| 8-16 GB | llama3.2:8b | `ollama pull llama3.2:8b` |
| 16 GB+ | qwen2.5:14b | `ollama pull qwen2.5:14b` |

If you pull a different model, also edit `DEFAULT_MODEL` at the top of `docking_pipeline.py` to match.

---

## Step 4: Install AutoDock Vina executable

The `vina` Python package is painful to install on Windows (requires Boost compilation). This tool uses the Vina **executable** instead.

1. Go to https://github.com/ccsb-scripps/AutoDock-Vina/releases
2. Download the latest Windows release asset (e.g. `vina_1.2.5_windows_x86_64.exe`)
3. Rename to `vina.exe` (drop the version suffix)
4. Put it in a permanent folder, e.g. `C:\Tools\vina\vina.exe`
5. Add `C:\Tools\vina` to your Windows PATH:
   - Press Win key -> type "env" -> "Edit the system environment variables"
   - Environment Variables -> under "User variables", find `Path` -> Edit
   - New -> paste `C:\Tools\vina`
   - OK on every dialog
6. **Close and reopen your command prompt** (PATH changes don't apply to existing windows)

Verify:

```bash
vina --version
```

Should print `AutoDock Vina v1.2.x`.

---

## Step 5: (Optional) Install Open Babel

Only needed if your workflow converts between chemistry file formats outside what RDKit and Meeko handle.

1. Download from https://openbabel.org/docs/Installation/install.html
2. Run installer (adds itself to PATH).

Verify: `obabel --version`

---

## Step 6: Clone this repository

```bash
git clone https://github.com/ManuHCI/multiagent-docking-gui.git
cd multiagent-docking-gui
```

If you don't have git, you can also download the repository as a ZIP from GitHub and extract.

---

## Step 7: Install Python dependencies

```bash
py -3 -m pip install -r requirements.txt
```

Takes 2-5 minutes. Installs PySide6, rdkit, meeko, scipy, gemmi, aisuite + Ollama SDK, and plotting libraries.

---

## Step 8: Run the GUI

```bash
py -3 docking_gui.py
```

You should see the app window. The log panel will show:

```
[ok] Ollama server reachable at http://localhost:11434
Project root: ...
Ready. Enter ligand + PDB IDs and click Run.
```

Test run: enter `Levetiracetam` and `6CM4`, click Run Docking. First LLM call takes 30-90 seconds (Ollama loads the model into RAM). Subsequent calls are fast. Full pipeline: ~5 minutes.

Output: `MultiAgent_Project/reports/` contains the Markdown report; `MultiAgent_Project/figures/` has plots.

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'X'`

A dependency is missing. Reinstall:

```bash
py -3 -m pip install -r requirements.txt
```

If a specific package keeps failing, install it individually:

```bash
py -3 -m pip install <package-name>
```

### `AttributeError: module 'ast' has no attribute 'NameConstant'`

You're on Python 3.14, which removed a deprecated API that some libraries still use. Options:

- Upgrade the offending library: `py -3 -m pip install --upgrade docstring_parser`
- Or install Python 3.12 alongside and use `py -3.12` instead of `py -3`

### Ollama: `404 Not Found` at `/api/chat`

The model you specified isn't pulled. Check with `ollama list`, then either pull it with `ollama pull <n>` or edit `DEFAULT_MODEL` in `docking_pipeline.py` to match what you have.

### Ollama: `500 Internal Server Error` or `timed out`

Two causes:

**Cause 1: out of memory.** Large models (7B+, especially 14B+) need a lot of free RAM. Switch to a smaller model:

```bash
ollama pull llama3.2:3b
```

Then edit `DEFAULT_MODEL = "ollama:llama3.2:3b"`.

**Cause 2: timeout too short.** First LLM call has to load the model into RAM (cold start can take 30-120 seconds). The pipeline includes a 600-second timeout setting in `Client()`, but if your model is very large, you may need longer. Edit `docking_pipeline.py`:

```python
CLIENT = Client({"ollama": {"timeout": 1200}})   # 20 minutes
```

### Docking step fails: `vina: command not found`

Step 4 didn't complete correctly. Make sure:

1. `vina.exe` is renamed (no version in the filename)
2. Its parent folder is in PATH
3. You opened a **new** command prompt after editing PATH
4. `vina --version` works from that new command prompt

### Pipeline runs but gives weird results

Not a bug -- a characteristic of small local LLMs. Try:

1. A larger model (8B+) for better planning and report quality
2. A model specifically trained for structured output: `ollama pull phi4-mini`
3. Manual verification of results against the raw Vina log files in `MultiAgent_Project/docking/`

---

## Environment notes

### Python version compatibility

- Python 3.10-3.13: **fully supported**
- Python 3.14: **known issues** with some scientific libraries (docstring_parser, possibly others). If you can, use 3.12 or 3.13 for this tool.

### Windows-specific

- Use `py -3` or `py -3.12` launcher, not bare `python`
- PATH changes require a fresh command prompt to take effect
- Windows Defender may briefly quarantine `vina.exe` on first run -- click "Allow"

### Linux / Mac

- Use `python3` in place of `py -3`
- Install Vina via your package manager (`apt install autodock-vina` on Ubuntu)
- PATH is set via `~/.bashrc` or `~/.zshrc`

---

## Offline use

After one-time setup (steps 1-7), the tool can run **fully offline**. The only network requests during a docking run are:

- PubChem lookup for a new ligand (one-shot, a few KB)
- RCSB download for a new receptor PDB (one-shot, <1 MB per structure)

All LLM inference happens locally. No data leaves your machine otherwise. This makes the tool suitable for institutional environments with restricted internet or for air-gapped computers (once you manually transfer the pulled Ollama model).

---

## Updating

```bash
cd multiagent-docking-gui
git pull
py -3 -m pip install -r requirements.txt --upgrade
```
