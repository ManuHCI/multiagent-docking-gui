"""Execute exported code against genuine prepared and docked example inputs."""
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import traceback
import zipfile

ROOT = Path(__file__).resolve().parent


def main():
    from core import research, study
    output = ROOT / ".verification" / "research-code"
    output.mkdir(parents=True, exist_ok=True)
    report = {"passed": False, "scope": "Execution checks only. MD limits minimisation to 200 iterations and uses 30 integration steps, not a scientific trajectory.", "checks": []}
    try:
        with tempfile.TemporaryDirectory(prefix="export_code_", dir=output) as temporary:
            case = study.create_example_study(Path(temporary))
            case.data["drugs"] = case.data["drugs"][:1]
            case.save()
            study.run_preparation(case)
            study.run_docking(case, exhaustiveness=1, n_decoys=1)
            target, drug = case.data["targets"][0], case.data["drugs"][0]
            for kind, script, extra in [("docking", "run_docking.py", []), ("md", "run_md.py", ["--smoke"])]:
                package = research.export_pair(case, target["id"], drug["id"], kind)
                folder = Path(temporary) / kind
                folder.mkdir()
                with zipfile.ZipFile(io.BytesIO(package)) as archive:
                    archive.extractall(folder)
                environment = os.environ.copy()
                environment["OPENMM_CPU_THREADS"] = "2"
                environment["PYTHONUNBUFFERED"] = "1"
                print("Executing exported " + kind + " code", flush=True)
                probe = "import faulthandler,runpy,sys; faulthandler.dump_traceback_later(180,repeat=True); sys.argv=sys.argv[1:]; runpy.run_path(sys.argv[0],run_name='__main__')"
                command = [sys.executable, "-u", "-c", probe, str(folder / script), *extra]
                subprocess.run(command, cwd=folder, env=environment, check=True,
                               timeout=int(os.environ.get("INSILICO_EXPORT_CHECK_TIMEOUT", "900")))
                expected = folder / "output" / ("results.json" if kind == "docking" else "final.pdb")
                if not expected.is_file():
                    raise AssertionError("The exported code did not produce " + expected.name)
                if kind == "md":
                    run_record = json.loads((folder / "output/run.json").read_text())
                    assert run_record["smoke_test_only"] and run_record["minimization_max_iterations"] == 200
                    assert (folder / "output/state.chk").is_file()
                    subprocess.run(command, cwd=folder, env=environment, check=True, timeout=120)
                    report["checks"].append("MD saved its actual system, checkpoint and protocol and resumed the completed execution check")
                report["checks"].append(kind + " code executed with exported inputs")
            report["passed"] = True
    except Exception as error:
        report.update(error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        (output / "verification.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("Exported docking and MD code execution checks passed.")


if __name__ == "__main__":
    main()
