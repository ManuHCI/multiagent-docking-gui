"""Use the same verified environment installer shipped in research downloads."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "assets/templates/bootstrap_micromamba.py"), run_name="__main__")
