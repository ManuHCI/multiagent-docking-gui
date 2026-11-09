"""Native desktop entry point. The bundled app never needs a system Python."""
from __future__ import annotations

import argparse
import base64
import html
import json
import logging
import multiprocessing
import os
from pathlib import Path
import runpy
import signal
import socket
import subprocess
import sys
import threading
import time
import traceback
import urllib.request

from core.runtime import configure, data_root, resource_root

TITLE = "In-Silico Docking Studio"
VERSION = "1.4.0"


def webview_backend() -> tuple[str | None, str]:
    if sys.platform == "darwin":
        return "cocoa", "Cocoa WebKit"
    if sys.platform == "win32":
        return "edgechromium", "Edge WebView2"
    return None, "System WebView"


def desktop_origin(port: int) -> str:
    """Use a local hostname for WebKit while keeping the server loopback-only."""
    host = "localhost" if sys.platform == "darwin" else "127.0.0.1"
    return f"http://{host}:{port}"


def native_page_state(window) -> dict:
    """Read primitive JSON through Cocoa rather than its native object bridge."""
    script = r"""JSON.stringify({
        url: location.href,
        document_title: document.title,
        ready_state: document.readyState,
        product_title: !!document.body && document.body.innerText.includes('In-Silico Docking Studio'),
        navigation: Array.from(document.querySelectorAll('[data-testid="primary-navigation"] button'))
            .map(button => {
                const label = button.querySelector('.q-btn__content .block');
                return (label ? label.innerText : button.innerText).trim();
            }),
        start_project: Array.from(document.querySelectorAll('button'))
            .some(button => button.innerText.includes('Start new project')),
        body_excerpt: document.body ? document.body.innerText.slice(0, 700) : '',
        html_excerpt: document.documentElement.outerHTML.slice(0, 700)
    })"""
    # Cocoa run_js returns the WebKit primitive string directly. It avoids
    # depending on pywebview's injected object-conversion bridge for this
    # acceptance probe; the application itself uses NiceGUI's own connection.
    raw = window.run_js(script) if sys.platform == "darwin" else window.evaluate_js(script)
    state = json.loads(raw) if isinstance(raw, str) else raw
    if not isinstance(state, dict):
        raise ValueError(f"The native page probe returned {type(raw).__name__}: {str(raw)[:160]}")
    return state


def native_page_ready(state: dict) -> bool:
    return (state.get("url", "").endswith("/study")
            and state.get("navigation") == ["Home", "Projects", "Reference", "About"]
            and bool(state.get("start_project"))
            and state.get("ready_state") in {"interactive", "complete"})


def server_process_options() -> dict:
    """Give the local server its own process group on macOS/POSIX."""
    if sys.platform == "win32":
        return {"creationflags": getattr(subprocess, "CREATE_NO_WINDOW", 0)}
    return {"start_new_session": True}


def stop_server(process: subprocess.Popen) -> None:
    """Stop this launcher's server and its docking children only."""
    if sys.platform == "win32":
        if process.poll() is None:
            taskkill = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "taskkill.exe"
            try:
                subprocess.run([str(taskkill), "/PID", str(process.pid), "/T", "/F"],
                    capture_output=True, timeout=10,
                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            except (OSError, subprocess.SubprocessError):
                pass
            if process.poll() is None:
                process.terminate()
    else:
        # start_new_session makes the child's PID its group ID. Using that
        # recorded ID also reaches Vina if the server has already exited.
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=8)
    except subprocess.TimeoutExpired:
        if sys.platform == "win32":
            process.kill()
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.wait(timeout=5)
    if sys.platform != "win32":
        # A child could ignore TERM after its parent exits; finish the same
        # owned group rather than leaving a docking process behind.
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def serve(port: int) -> None:
    os.environ["INSILICO_HOST"] = "127.0.0.1"
    os.environ["INSILICO_PORT"] = str(port)
    os.environ.setdefault("NICEGUI_STORAGE_PATH", str(data_root() / "Session"))
    runpy.run_module("app", run_name="__main__")


def diagnostic(output: str = "") -> int:
    from dataclasses import asdict
    from core.selftest import print_report, run
    report = run(check_network=False)
    if sys.stdout is not None:
        print_report(report)
    result = {"version": VERSION, "can_start": report.can_start,
              "machine": report.machine,
              "findings": [asdict(item) for item in report.findings]}
    if output:
        Path(output).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if report.can_start else 1


def splash(message: str = "Opening your docking workspace…", error: bool = False) -> str:
    try:
        encoded_logo = base64.b64encode((resource_root() / "assets" / "logo.png").read_bytes()).decode("ascii")
        logo = '<img alt="Akshaisarathi Healthcare logo" style="width:360px;max-width:100%;margin-bottom:25px" src="data:image/png;base64,' + encoded_logo + '">'
    except OSError:
        logo = '<div class="mark">✚</div>'
    return '''<!doctype html><html><head><meta charset="utf-8"><style>
    body{margin:0;background:#f1f6f8;color:#173448;font:17px 'Segoe UI',sans-serif;
    display:grid;place-items:center;min-height:100vh}main{max-width:610px;padding:48px}
    .mark{font-size:42px;color:#087f83}small{color:#087f83;letter-spacing:2px;font-weight:700}
    h1{font-size:36px;line-height:1.2}p{line-height:1.7;color:#506474}
    .bar{height:5px;background:#d5e6e9;border-radius:9px;overflow:hidden;margin-top:32px}
    .bar:after{content:'';display:block;width:40%;height:100%;background:#087f83;
    animation:move 1.4s infinite alternate}@keyframes move{to{transform:translateX(150%)}}
    </style></head><body><main>''' + logo + '''<br><small>MEDICAL EDUCATION</small>
    <h1>In-Silico<br>Docking Studio</h1><p>Add targets and drug molecules.<br>
    Prepare, dock and review results with simple steps.</p><p>''' + html.escape(message) + '''</p>''' + (
        "" if error else '<div class="bar"></div>') + '</main></body></html>'


def open_desktop(port: int, test_output: str = "") -> int:
    import webview
    # Worksheets and exported MD packages use the native Save dialog.
    webview.settings["ALLOW_DOWNLOADS"] = True

    root = data_root()
    logs = root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log_path = logs / "desktop.log"
    command = ([sys.executable] if getattr(sys, "frozen", False)
               else [sys.executable, str(Path(__file__).resolve())])
    env = os.environ.copy()
    env["INSILICO_DATA_DIR"] = str(root)
    # Keep a copied developer environment from changing the private runtime.
    if getattr(sys, "frozen", False):
        for name in ("PYTHONHOME", "PYTHONPATH"):
            env.pop(name, None)
    output = log_path.open("a", encoding="utf-8")
    process = subprocess.Popen(
        command + ["--serve", "--port", str(port)], cwd=resource_root(), env=env,
        stdout=output, stderr=subprocess.STDOUT,
        **server_process_options())
    stopped = False
    test_passed = False
    stop_lock = threading.Lock()

    def stop(*_args):
        nonlocal stopped
        # WebView2 dispatches the closed callback on another thread. The
        # main-thread finally block must wait for that cleanup to finish.
        with stop_lock:
            if stopped:
                return
            stopped = True
            stop_server(process)
            output.close()

    window = webview.create_window(
        TITLE, html=splash(), width=1280, height=860, min_size=(960, 640),
        background_color="#f1f6f8", text_select=True)
    window.events.closed += stop
    url = desktop_origin(port)

    def ready():
        nonlocal test_passed
        deadline = time.monotonic() + 120
        while not stopped and time.monotonic() < deadline:
            if process.poll() is not None:
                break
            try:
                with urllib.request.urlopen(url + "/assets/logo.png", timeout=1) as response:
                    if response.status == 200:
                        window.load_url(url + "/study")
                        if test_output:
                            last_state = {}
                            last_error = ""
                            test_deadline = time.monotonic() + 45
                            while not stopped and time.monotonic() < test_deadline:
                                try:
                                    last_state = native_page_state(window)
                                    if native_page_ready(last_state):
                                        test_passed = True
                                        break
                                except Exception as error:
                                    last_error = f"{type(error).__name__}: {error}"
                                time.sleep(0.25)
                            page_events = {name: getattr(window.events, name).is_set()
                                           for name in ("shown", "loaded", "before_load")}
                            current_url = ""
                            if page_events["loaded"]:
                                try:
                                    current_url = window.get_current_url() or ""
                                except Exception as error:
                                    last_error = f"{type(error).__name__}: {error}"
                            screenshot = ""
                            screenshot_error = ""
                            if sys.platform == "darwin":
                                capture = Path(test_output).with_suffix(".png")
                                try:
                                    result = subprocess.run(["/usr/sbin/screencapture", "-x", str(capture)],
                                        capture_output=True, text=True, timeout=10)
                                    if result.returncode == 0 and capture.is_file():
                                        screenshot = str(capture)
                                    else:
                                        screenshot_error = (result.stderr or f"exit {result.returncode}")[:500]
                                except (OSError, subprocess.SubprocessError) as error:
                                    screenshot_error = str(error)
                            Path(test_output).write_text(json.dumps(
                                {"passed": test_passed, "native_window": True,
                                 "engine": webview_backend()[1], "route": "/study",
                                 "downloads_enabled": bool(webview.settings["ALLOW_DOWNLOADS"]),
                                 "page_state": last_state, "page_events": page_events,
                                 "current_url": current_url, "last_probe_error": last_error,
                                 "screenshot": screenshot, "screenshot_error": screenshot_error}, indent=2), encoding="utf-8")
                            window.destroy()
                        return
            except (OSError, ValueError):
                time.sleep(0.2)
        if not stopped:
            window.load_html(splash(
                "The learning workspace could not open. Close this window and try again. "
                f"The diagnostic log is saved at {log_path}", error=True))
            if test_output:
                Path(test_output).write_text(json.dumps({"passed": False, "log": str(log_path)}), encoding="utf-8")
                window.destroy()

    try:
        webview.start(ready, gui=webview_backend()[0], private_mode=False,
                      storage_path=str(root / "window"), debug=False)
    finally:
        stop()
        if test_output and Path(test_output).is_file():
            result = json.loads(Path(test_output).read_text(encoding="utf-8"))
            result["server_stopped"] = process.poll() is not None
            Path(test_output).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if not test_output or (test_passed and process.poll() is not None) else 1


def main() -> int:
    if getattr(sys, "frozen", False):
        logs = data_root() / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        stream = (logs / "runtime.log").open("a", encoding="utf-8", buffering=1)
        sys.stdout = sys.stderr = stream
    configure()
    parser = argparse.ArgumentParser(description=TITLE)
    parser.add_argument("--serve", "--headless", dest="serve", action="store_true")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--diagnose", action="store_true")
    parser.add_argument("--verify-output", default="", metavar="JSON_FILE")
    parser.add_argument("--verify-science", default="", metavar="JSON_FILE")
    parser.add_argument("--verify-study", default="", metavar="JSON_FILE")
    parser.add_argument("--window-test-output", default="", metavar="JSON_FILE", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.verify_study:
        from core.study_checks import run_study_checks
        return run_study_checks(args.verify_study)
    if args.verify_science:
        from core.desktop_checks import run_scientific_checks
        return run_scientific_checks(args.verify_science)
    if args.diagnose or args.verify_output:
        return diagnostic(args.verify_output)
    if args.serve:
        serve(args.port or free_port())
        return 0
    return open_desktop(args.port or free_port(), args.window_test_output)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    try:
        raise SystemExit(main())
    except Exception:
        error_root = data_root() / "logs"
        error_root.mkdir(parents=True, exist_ok=True)
        (error_root / "startup-error.log").write_text(traceback.format_exc(), encoding="utf-8")
        command_mode = any(arg in sys.argv for arg in
                           ("--serve", "--headless", "--diagnose", "--verify-output", "--verify-science", "--verify-study", "--window-test-output"))
        if os.name == "nt" and not command_mode:
            import ctypes
            ctypes.windll.user32.MessageBoxW(
                0, "The teaching app could not open. Please reinstall the app if files "
                "were moved, and check that Microsoft Edge WebView2 Runtime is installed.\n\n"
                f"Details: {error_root / 'startup-error.log'}", TITLE, 0x10)
        if sys.stderr is not None:
            traceback.print_exc()
        raise SystemExit(1)
