#!/usr/bin/env python3
"""Download only the verified BGC↔MS2 spectral links from PoDP JSON metadata.

This intentionally skips genome FASTA/GBK downloads and only fetches the MS2 files
from MassIVE when the BGC-MS2 entry has a real MS2 pointer and non-empty verification.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.podp_online_file_locator import collect_verified_bgc_ms2_manifest, is_verified_bgc_ms2_link

DEFAULT_MANIFEST = ROOT / "data" / "PoDP" / "verified_bgc_ms2_manifest.json"
OUT_ROOT = ROOT / "data" / "PoDP" / "verified_bgc_ms2_downloads"


def safe_filename(value: str) -> str:
    cleaned = value.strip().replace("/", "_").replace(" ", "_")
    return "".join(ch for ch in cleaned if ch.isalnum() or ch in "._-") or "unknown"


def resolve_downloadable_url(url: str) -> str:
    """Rewrite a raw ftp://massive.ucsd.edu MS2 URI to MassIVE's HTTPS gateway.

    Plain FTP is blocked in this environment (and increasingly elsewhere), but
    MassIVE mirrors every FTP path through ProteoSAFe's DownloadResultFile
    endpoint over HTTPS, so we translate rather than hit the FTP port.
    """
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme == "ftp" and parsed.hostname == "massive.ucsd.edu":
        # unquote first: urlsplit leaves percent-escapes (e.g. "%20") intact in
        # .path, so encoding it again with urlencode below would double-escape them.
        file_param = "f." + urllib.parse.unquote(parsed.path.lstrip("/"))
        query = urllib.parse.urlencode({"forceDownload": "true", "file": file_param})
        return f"https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?{query}"
    return url


def download_url(url: str, destination: Path, retries: int = 6, backoff: float = 5.0) -> tuple[bool, str | None]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return True, None
    resolved = resolve_downloadable_url(url)
    last_error: str | None = None
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(resolved, destination)
            return True, None
        except urllib.error.HTTPError as exc:
            last_error = str(exc)
            if exc.code == 429 and attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
                continue
            return False, last_error
        except Exception as exc:  # pragma: no cover - network/runtime issue
            return False, str(exc)
    return False, last_error


def iter_manifest(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        return data.get("projects", data.get("verified_bgc_ms2_links", []))
    return data


def main() -> int:
    manifest_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_MANIFEST
    out_root = Path(sys.argv[2]) if len(sys.argv) > 2 else OUT_ROOT

    manifest = iter_manifest(manifest_path)
    rows: list[dict[str, str]] = []

    for project in manifest:
        project_id = str(project.get("project_id") or project.get("project_file", "unknown")).strip()
        project_dir = out_root / safe_filename(project_id)
        project_dir.mkdir(parents=True, exist_ok=True)

        project_links = project.get("verified_bgc_ms2_links") or []
        if not project_links:
            continue

        for link in project_links:
            if not is_verified_bgc_ms2_link(link):
                continue
            ms2_url = link.get("MS2_URL")
            if not ms2_url:
                continue
            dest_name = safe_filename(Path(ms2_url.split("?")[0]).name)
            status, error = download_url(str(ms2_url), project_dir / dest_name)
            time.sleep(1.5)
            rows.append({
                "project_id": project_id,
                "kind": "ms2",
                "ms2_url": str(ms2_url),
                "status": "ok" if status else "failed",
                "error": error or "",
            })

    report_path = out_root / "download_report.csv"
    with report_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["project_id", "kind", "ms2_url", "status", "error"])
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps({
        "manifest": str(manifest_path),
        "output_root": str(out_root),
        "projects_processed": len(manifest),
        "download_rows": len(rows),
        "report": str(report_path),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
