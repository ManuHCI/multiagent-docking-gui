"""
The one place this software talks to the internet.

Everything goes through here so that four things are true everywhere:

  * Responses are cached to disk. The specification requires the
    software to open and teach with no connection, so anything fetched
    once stays available. A workshop room whose wifi fails mid-session
    keeps working.
  * Failures are sentences, not exceptions. A resident should read
    "PubChem did not answer - you can carry on offline", never a
    traceback.
  * Timeouts are short. A hung request on a slow college connection is
    worse than a failed one, because the user cannot tell the
    difference between slow and broken.
  * Every fetch can report what it fetched and when, for the
    provenance record.

No third-party HTTP library: urllib is in the standard library, which
is one less thing the installer has to carry and one less thing to
break.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

USER_AGENT = "InSilicoTeachingSuite/0.1 (academic teaching tool)"
DEFAULT_TIMEOUT = 15.0


@dataclass
class Fetched:
    """A response, or a plain-English explanation of why there isn't one."""

    ok: bool
    url: str
    status: int = 0
    text: str = ""
    from_cache: bool = False
    retrieved: str = ""
    problem: str = ""          # what the user reads
    suggestion: str = ""       # what they can do

    def json(self) -> Any:
        if not self.ok:
            return None
        try:
            return json.loads(self.text)
        except json.JSONDecodeError:
            return None


class Web:
    """Fetches, caches, and never raises."""

    def __init__(self, cache_dir: Path, offline: bool = False) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline

    # -- cache ------------------------------------------------------------

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha256(url.encode()).hexdigest()[:20]
        return self.cache_dir / f"{digest}.json"

    def _read_cache(self, url: str) -> Fetched | None:
        path = self._cache_path(url)
        if not path.exists():
            return None
        try:
            blob = json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            return None
        return Fetched(
            ok=True,
            url=url,
            status=blob.get("status", 200),
            text=blob.get("text", ""),
            from_cache=True,
            retrieved=blob.get("retrieved", ""),
        )

    def _write_cache(self, url: str, status: int, text: str, when: str) -> None:
        try:
            self._cache_path(url).write_text(
                json.dumps({"url": url, "status": status,
                            "text": text, "retrieved": when}),
                encoding="utf-8",
            )
        except Exception:  # noqa: BLE001
            pass   # a cache that cannot be written is not a reason to fail

    # -- fetching ---------------------------------------------------------

    def get(
        self,
        url: str,
        service: str = "that service",
        timeout: float = DEFAULT_TIMEOUT,
        use_cache: bool = True,
        retries: int = 1,
    ) -> Fetched:
        if use_cache:
            cached = self._read_cache(url)
            if cached and self.offline:
                return cached

        if self.offline:
            return Fetched(
                ok=False, url=url,
                problem=f"This is the first time you have asked for this from "
                        f"{service}, and the software is running offline.",
                suggestion="Connect to the internet and try again.",
            )

        last: Exception | None = None
        for attempt in range(retries + 1):
            try:
                req = urllib.request.Request(
                    url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"}
                )
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                    when = datetime.now(timezone.utc).isoformat(timespec="seconds")
                    if use_cache:
                        self._write_cache(url, resp.status, text, when)
                    return Fetched(ok=True, url=url, status=resp.status,
                                   text=text, retrieved=when)

            except urllib.error.HTTPError as exc:
                # The server answered, it just said no. Retrying will not help.
                if exc.code == 404:
                    return Fetched(
                        ok=False, url=url, status=404,
                        problem=f"{service} has no record of that.",
                        suggestion="Check the identifier you entered.",
                    )
                if exc.code == 429:
                    return Fetched(
                        ok=False, url=url, status=429,
                        problem=f"{service} is asking us to slow down.",
                        suggestion="Wait a minute and try again.",
                    )
                return Fetched(
                    ok=False, url=url, status=exc.code,
                    problem=f"{service} returned an error ({exc.code}).",
                    suggestion="Try again shortly, or carry on without it.",
                )

            except Exception as exc:  # noqa: BLE001
                last = exc
                if attempt < retries:
                    time.sleep(1.5)

        # Could not reach it at all. Fall back to anything cached.
        cached = self._read_cache(url) if use_cache else None
        if cached:
            cached.problem = (
                f"Could not reach {service}, so an earlier copy is being "
                f"used (fetched {cached.retrieved})."
            )
            return cached

        return Fetched(
            ok=False, url=url,
            problem=f"Could not reach {service}.",
            suggestion=(
                "Check your internet connection. You can carry on with the "
                "worked examples offline."
            ),
            text=type(last).__name__ if last else "",
        )

    def get_json(self, url: str, service: str = "that service",
                 **kwargs: Any) -> tuple[Any, Fetched]:
        result = self.get(url, service=service, **kwargs)
        return result.json(), result

    def post(
        self,
        url: str,
        data: bytes,
        service: str = "that service",
        content_type: str = "application/json",
        timeout: float = DEFAULT_TIMEOUT,
    ) -> Fetched:
        """
        POST, for GraphQL. Not cached -- the body varies.

        One thing this gets right that the first version of the
        project's own diagnostic got wrong: a service that answers with
        an error is reported as having ANSWERED. "Could not reach it"
        and "it refused" are different problems with different fixes,
        and conflating them sends you looking in the wrong place. The
        body of an error response is kept, because for GraphQL that is
        where the server tells you exactly what was wrong.
        """
        if self.offline:
            return Fetched(
                ok=False, url=url,
                problem=f"{service} needs a connection and the software is "
                        "running offline.",
                suggestion="Connect to the internet and try again.",
            )

        req = urllib.request.Request(
            url, data=data, method="POST",
            headers={"User-Agent": USER_AGENT,
                     "Content-Type": content_type,
                     "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return Fetched(
                    ok=True, url=url, status=resp.status,
                    text=resp.read().decode("utf-8", errors="replace"),
                    retrieved=datetime.now(timezone.utc).isoformat(
                        timespec="seconds"),
                )
        except urllib.error.HTTPError as exc:
            body = ""
            try:
                body = exc.read().decode("utf-8", errors="replace")
            except Exception:  # noqa: BLE001
                pass
            # Answered, but refused. Keep the body: it says why.
            return Fetched(
                ok=False, url=url, status=exc.code, text=body,
                problem=f"{service} refused the request (HTTP {exc.code}).",
                suggestion="This is usually a fault in the request, not the "
                           "connection.",
            )
        except Exception as exc:  # noqa: BLE001
            return Fetched(
                ok=False, url=url,
                problem=f"Could not reach {service}.",
                suggestion="Check your internet connection.",
                text=type(exc).__name__,
            )
