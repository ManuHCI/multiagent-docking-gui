"""Offline acceptance of the four-section docking studio and real calculations."""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/".build-tools"))

from desktop import free_port
from playwright.sync_api import sync_playwright


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable",type=Path)
    args = parser.parse_args()
    output = ROOT/".verification"/"studio"
    output.mkdir(parents=True,exist_ok=True)
    temporary = ROOT/".qa-temp"
    temporary.mkdir(exist_ok=True)
    case = Path(tempfile.mkdtemp(prefix="studio-ui-",dir=temporary))
    user = case/"student data"
    from core.project import Project
    legacy_folder = user/"Projects"/"Legacy archive"
    legacy = Project.create(legacy_folder, "Legacy archive")
    legacy.close()
    port = free_port()
    base = f"http://127.0.0.1:{port}"
    env = os.environ.copy()
    env["INSILICO_DATA_DIR"] = str(user)
    env.pop("INSILICO_PROJECTS",None)
    if args.executable:
        command = [str(args.executable.resolve())]
        for name in ("PYTHONPATH","PYTHONHOME","CONDA_PREFIX","CONDA_DEFAULT_ENV","RDBASE","OPENMM_PLUGIN_DIR"):
            env.pop(name,None)
        env["PATH"] = r"C:\Windows\System32;C:\Windows"
    else:
        command = [sys.executable,str(ROOT/"desktop.py")]
    checks = []
    browser_errors = []
    report = {"passed":False,"standalone":bool(args.executable),"checks":checks,
              "browser_errors":browser_errors,"user_data_dir":str(user.resolve())}
    started = time.monotonic()

    def no_questions(page):
        text = page.locator("body").inner_text()
        prohibited = ("Your prediction", "What did you observe", "Explain in your own words", "discussed the question")
        if any(phrase.lower() in text.lower() for phrase in prohibited):
            raise AssertionError("The practical workspace contains a quiz prompt")
        if "?" in text:
            location = text.index("?")
            raise AssertionError("Question text remains in the workflow: " + text[max(0,location-100):location+30])
        if page.locator(".q-notification.bg-negative:visible").count() or "parent element this slot belongs to has been deleted" in text:
            raise AssertionError("An application error notification is visible: " + text[-700:])

    def screenshot(page,name,full_page=False):
        page.evaluate("window.scrollTo(0,0)")
        page.wait_for_timeout(450)
        page.screenshot(path=output/name,full_page=full_page,animations="disabled")

    def navigation(page, section):
        header = page.locator("header")
        header.get_by_text("In-Silico Docking Studio",exact=True).wait_for()
        expected = ["Home", "Projects", "Reference", "About"]
        positions = []
        for label in expected:
            button = header.get_by_role("button",name=label,exact=True)
            button.wait_for()
            positions.append(button.bounding_box()["x"])
        if positions != sorted(positions):
            raise AssertionError("The fixed navigation is not in Home, Projects, Reference, About order")
        text = header.inner_text().lower()
        if any(label in text for label in ("guided docking", "advanced workflow", "regular workflow")):
            raise AssertionError("Competing workflow navigation remains in the header")
        header.get_by_role("button",name=section,exact=True).click()
        header.get_by_text("In-Silico Docking Studio",exact=True).wait_for()

    with (output/"server.log").open("w",encoding="utf-8") as log:
        server = subprocess.Popen(command+["--serve","--port",str(port)],cwd=case,
                                  env=env,stdout=log,stderr=log,
                                  creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0))
        browser = None
        try:
            deadline = time.monotonic()+90
            while True:
                if server.poll() is not None:
                    raise RuntimeError(f"The study server exited with {server.returncode}")
                try:
                    with urllib.request.urlopen(base+"/assets/logo.png",timeout=1):
                        break
                except OSError:
                    if time.monotonic()>deadline:
                        raise TimeoutError("The docking studio did not start")
                    time.sleep(.2)
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch(channel="msedge",headless=True)
                page = browser.new_page(viewport={"width":1280,"height":860},device_scale_factor=1)
                page.set_default_timeout(30000)
                page.on("pageerror",lambda error:browser_errors.append(str(error)))
                page.route("**/*",lambda route:route.continue_() if
                           route.request.url.startswith((base,"data:","blob:")) else route.abort())
                page.goto(base+"/")
                page.wait_for_url(base+"/study")
                page.get_by_role("button",name="Use the offline EGFR example",exact=True).wait_for()
                navigation(page,"Home")
                logo = page.locator("header img[alt='Akshaisarathi Healthcare logo']").first
                logo.wait_for()
                logo.evaluate("image => image.decode()")
                if (logo.bounding_box() or {}).get("width",0)<100:
                    raise AssertionError("The product logo is not visibly sized")
                no_questions(page)
                screenshot(page,"home.png",True)
                checks.append("canonical Home with fixed four-section navigation and visible supplied logo")

                navigation(page,"About")
                page.get_by_text("Dr Manu Kumar Shetty, Professor, Department of Pharmacology, Maulana Azad Medical College, New Delhi",exact=True).wait_for()
                page.get_by_text("Dr Pratibha Nadig, Founder, Akshaisarathi Healthcare",exact=True).wait_for()
                screenshot(page,"about.png")
                checks.append("both requested developer credits")

                navigation(page,"Reference")
                page.get_by_role("link",name="RCSB Protein Data Bank",exact=True).wait_for()
                page.get_by_role("link",name="PubChem drug and compound structures",exact=True).wait_for()
                if not page.locator("a[href='https://doi.org/10.1021/acs.jcim.1c00203']").count():
                    raise AssertionError("The Vina reference DOI is missing")
                with page.expect_download() as download_info:
                    page.get_by_role("button",name="Download references and methods",exact=True).click()
                download_info.value.save_as(output/"references-catalogue.txt")
                text = (output/"references-catalogue.txt").read_text(encoding="utf-8")
                if "AutoDock Vina" not in text or "https://www.rcsb.org/" not in text:
                    raise AssertionError("The downloaded references are incomplete")
                checks.append("protein/drug download links, Vina DOI and downloadable references")

                navigation(page,"Projects")
                legacy_card = page.locator('[data-testid="project-card"]').filter(
                    has=page.get_by_text("Legacy archive",exact=True))
                legacy_card.wait_for()
                legacy_card.get_by_role("button",name="Open project",exact=True).click()
                page.wait_for_url(re.compile(re.escape(base)+r"/legacy\?project=Legacy(?:%20|\+)archive"))
                page.get_by_text("Legacy archive",exact=True).first.wait_for()
                navigation(page,"Projects")
                page.wait_for_url(re.compile(re.escape(base)+r"/study\?section=projects"))
                legacy_card.wait_for()
                legacy_card.get_by_role("button",name="Delete project",exact=True).click()
                page.get_by_role("dialog").get_by_role("button",name="Delete project",exact=True).click()
                deadline = time.monotonic()+10
                while legacy_folder.exists() and time.monotonic()<deadline:
                    page.wait_for_timeout(100)
                if legacy_folder.exists():
                    raise AssertionError("The closed legacy project could not be deleted")
                legacy_card.wait_for(state="detached")
                checks.append("older individual projects reopen with the same navigation and close before deletion")
                navigation(page,"Home")
                page.get_by_role("button",name="Start new project",exact=True).click()
                page.get_by_label("Study name",exact=True).fill("Disposable matrix")

                # Exercise a genuine two-target/two-drug input matrix without
                # contacting online protein services during an offline test.
                for pdb_id,label in [("1M17","EGFR input"),("3ERT","ESR1 input")]:
                    page.get_by_label("PDB entry ID",exact=True).fill(pdb_id)
                    page.get_by_label("Target label",exact=True).fill(label)
                    page.get_by_role("button",name="Add protein target",exact=True).click()
                    page.get_by_text(label,exact=True).wait_for()
                page.get_by_role("button",name="Continue to drug molecules",exact=True).click()
                page.get_by_label("SMILES and drug name",exact=True).fill("CC(=O)Oc1ccccc1C(=O)O Aspirin\nCCO Ethanol")
                page.get_by_role("button",name="Add listed molecules",exact=True).click()
                page.get_by_text("2 drug molecules",exact=True).wait_for()
                page.get_by_role("button",name="Save study and prepare structures",exact=True).click()
                page.get_by_role("button",name="Prepare proteins and drugs",exact=True).wait_for()
                manifests = list((user/"Projects").glob("*/study.json"))
                matrix = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
                if not any(study["name"] == "Disposable matrix" and len(study["targets"])==2 and len(study["drugs"])==2 for study in matrix):
                    raise AssertionError("The two-by-two study input matrix was not persisted")
                disposable = next(path.parent for path in manifests if json.loads(path.read_text(encoding="utf-8"))["name"]=="Disposable matrix")
                no_questions(page)
                checks.append("two targets and two drugs saved as a four-pair study")
                navigation(page,"Reference")
                navigation(page,"About")
                navigation(page,"Home")
                page.get_by_role("button",name="Continue project",exact=True).click()
                page.get_by_test_id("project-step").filter(has_text="Step 3 of 5").wait_for()
                page.get_by_role("button",name="Prepare proteins and drugs",exact=True).wait_for()
                checks.append("active project and current step survive global navigation")
                navigation(page,"Home")

                page.get_by_role("button",name="Use the offline EGFR example",exact=True).click()
                page.get_by_role("button",name="Prepare proteins and drugs",exact=True).wait_for()
                page.get_by_role("button",name="Prepare proteins and drugs",exact=True).click()
                page.get_by_role("button",name="Continue to docking",exact=True).wait_for(state="visible",timeout=90000)
                deadline = time.monotonic()+90
                while not page.get_by_role("button",name="Continue to docking",exact=True).is_enabled():
                    if time.monotonic()>deadline:
                        raise TimeoutError("The prepared example did not become ready for docking")
                    page.wait_for_timeout(250)
                page.frame_locator("iframe").locator("#status").filter(has_text="3D ready").wait_for(timeout=30000)
                no_questions(page)
                page.evaluate("window.scrollTo(0,0)")
                screenshot(page,"preparation.png",True)
                checks.append("real offline molecule preparation, per-target binding region and local 3D")

                page.get_by_role("button",name="Continue to docking",exact=True).click()
                page.get_by_label("Search effort",exact=True).fill("1")
                page.get_by_label("Comparison molecules per target",exact=True).fill("1")
                page.get_by_role("button",name="Run docking study",exact=True).click()
                page.get_by_role("button",name="Download results CSV",exact=True).wait_for(timeout=180000)
                no_questions(page)
                with page.expect_download() as download_info:
                    page.get_by_role("button",name="Download results CSV",exact=True).click()
                download_info.value.save_as(output/"results.csv")
                csv_text = (output/"results.csv").read_text(encoding="utf-8-sig")
                rows = list(csv.DictReader(io.StringIO(csv_text)))
                if len(rows)!=3 or any(row["status"]!="completed" or not math.isfinite(float(row["score_kcal_mol"])) for row in rows):
                    raise AssertionError("The real example did not produce three complete finite docking scores")
                page.frame_locator("iframe").locator("#status").filter(has_text="3D ready").wait_for(timeout=30000)
                page.evaluate("window.scrollTo(0,0)")
                screenshot(page,"results.png",True)
                checks.append("fresh three-drug docking with finite Vina scores and current checks")

                with page.expect_download(timeout=60000) as download_info:
                    page.get_by_role("button",name="Download complete study",exact=True).click()
                download_info.value.save_as(output/"complete-study.zip")
                with zipfile.ZipFile(output/"complete-study.zip") as package:
                    if not {"results.csv","methods_and_checks.txt","references.txt","study.json"}.issubset(package.namelist()):
                        raise AssertionError("The complete study export is missing its methods or references")
                    methods = package.read("methods_and_checks.txt").decode("utf-8")
                    if "AutoDock Vina" not in methods or "Search box:" not in methods:
                        raise AssertionError("The exported methods omit actual docking provenance")
                checks.append("CSV, full project package, actual methods and references exported")
                for index in range(1,6):
                    page.get_by_role("button",name=re.compile(rf"^{index}\. ")).click()
                    page.get_by_test_id("project-step").filter(has_text=f"Step {index} of 5").wait_for()
                    no_questions(page)
                checks.append("all five workflow steps remain navigable")

                # Saved projects are listed only when that section is opened.
                # Delete a disposable input matrix, while retaining the real
                # calculated example and its complete result files.
                navigation(page,"Projects")
                cards = page.locator('[data-testid="project-card"]')
                disposable_card = cards.filter(has=page.get_by_text("Disposable matrix",exact=True))
                example_card = cards.filter(has=page.get_by_text("EGFR example",exact=True))
                disposable_card.wait_for()
                example_card.wait_for()
                if not re.search(r"\d+(?:\.\d+)?\s*(?:B|KB|MB|GB|KiB|MiB|GiB)\b",disposable_card.inner_text()):
                    raise AssertionError("Saved projects do not display their disk usage")
                screenshot(page,"projects.png",True)
                disposable_card.get_by_role("button",name="Open project",exact=True).click()
                page.get_by_test_id("project-step").filter(has_text="Step 3 of 5").wait_for()
                navigation(page,"Projects")
                example_card.get_by_role("button",name="Open project",exact=True).click()
                page.get_by_test_id("project-step").filter(has_text="Step 5 of 5").wait_for()
                navigation(page,"Projects")
                disposable_card.get_by_role("button",name="Delete project",exact=True).click()
                dialog = page.get_by_role("dialog")
                dialog.get_by_text("Disposable matrix",exact=False).wait_for()
                dialog.get_by_role("button",name="Cancel",exact=True).click()
                if not (disposable/"study.json").is_file():
                    raise AssertionError("Cancelling deletion removed a project")
                disposable_card.get_by_role("button",name="Delete project",exact=True).click()
                dialog.get_by_role("button",name="Delete project",exact=True).click()
                deadline = time.monotonic()+10
                while disposable.exists() and time.monotonic()<deadline:
                    page.wait_for_timeout(100)
                disposable_card.wait_for(state="detached")
                if disposable.exists():
                    raise AssertionError("Confirmed deletion did not remove the selected project's files")
                remaining = list((user/"Projects").glob("*/study.json"))
                if len(remaining)!=1 or json.loads(remaining[0].read_text(encoding="utf-8"))["name"]!="EGFR example":
                    raise AssertionError("Deleting one project affected the retained example")
                navigation(page,"Home")
                page.get_by_role("button",name="Continue project",exact=True).click()
                page.get_by_test_id("project-step").filter(has_text="Step 5 of 5").wait_for()
                page.get_by_role("button",name="Download results CSV",exact=True).wait_for()
                no_questions(page)
                checks.append("Projects shows disk usage, reopens saved steps and deletes only the confirmed project")

                page.goto(base+"/workshop")
                page.wait_for_url(base+"/study")
                page.get_by_role("button",name="Use the offline EGFR example",exact=True).wait_for()
                checks.append("legacy workshop redirects to the practical workspace")
                if browser_errors:
                    raise AssertionError("Browser errors: "+"; ".join(browser_errors))
                log.flush()
                if "Traceback (most recent call last)" in (output/"server.log").read_text(encoding="utf-8",errors="replace"):
                    raise AssertionError("An unhandled application exception was recorded in the server log")
                browser.close()
                browser = None
            report["passed"] = True
            print(f"Docking Studio UI acceptance passed: {output}",flush=True)
            return 0
        except Exception as error:
            report["error"] = f"{type(error).__name__}: {error}"
            try:
                page.screenshot(path=output/"failure.png",full_page=True,animations="disabled")
                (output/"failure-page.txt").write_text(page.locator("body").inner_text(),encoding="utf-8")
            except Exception:
                pass
            raise
        finally:
            report["seconds"] = round(time.monotonic()-started,3)
            (output/"ui-acceptance.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
            if server.poll() is None:
                if os.name == "nt":
                    subprocess.run([str(Path(os.environ.get("SystemRoot",r"C:\Windows"))/"System32"/"taskkill.exe"),
                                    "/PID",str(server.pid),"/T","/F"],
                                   stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT,timeout=20,
                                   creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0))
                else:
                    server.terminate()
                try:
                    server.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    server.kill()


if __name__=="__main__":
    raise SystemExit(main())
