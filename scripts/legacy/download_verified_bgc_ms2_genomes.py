#!/usr/bin/env python3
"""Download the genome GBK files that pair with the verified BGC-MS2 MS2 downloads.

Each verified BGC-MS2 link is matched to exactly one genome via
`genome_metabolome_links` (see podp_online_file_locator.collect_verified_bgc_ms2_genome_links),
so this script:

1. Downloads one GBK per unique (project_id, genome_accession) into
   data/PoDP/verified_bgc_ms2_downloads/<project_id>/genes/<accession>.gbk
   -- next to that project's MS2 files, downloaded by download_verified_bgc_ms2.py.
2. Writes data/PoDP/verified_bgc_ms2_downloads/genome_ms2_link_manifest.csv, a
   flat table with one row per verified BGC-MS2 link giving the exact MS2 file
   and genome file that belong together.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.podp_online_file_locator import collect_verified_bgc_ms2_genome_manifest

OUT_ROOT = ROOT / "data" / "PoDP" / "verified_bgc_ms2_downloads"
PODP_DIR = ROOT / "data" / "PoDP" / "json_descriptors"
LINK_MANIFEST_PATH = OUT_ROOT / "genome_ms2_link_manifest.csv"

NCBI_DELAY_SECONDS = 0.4  # stay under NCBI's unauthenticated 3-requests/sec limit


def safe_filename(value: str) -> str:
    cleaned = value.strip().replace("/", "_").replace(" ", "_")
    return "".join(ch for ch in cleaned if ch.isalnum() or ch in "._-") or "unknown"


MIN_VALID_BYTES = 10  # eutils returns a near-empty 200 OK for records with no sequence data


def download_url(url: str, destination: Path, retries: int = 5, backoff: float = 3.0) -> tuple[bool, str | None]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size >= MIN_VALID_BYTES:
        return True, None
    last_error: str | None = None
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(url, destination)
            if destination.stat().st_size < MIN_VALID_BYTES:
                destination.unlink()
                return False, "empty response (no sequence data at this accession)"
            return True, None
        except urllib.error.HTTPError as exc:
            last_error = str(exc)
            if exc.code in (429, 500, 502, 503) and attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
                continue
            return False, last_error
        except Exception as exc:  # pragma: no cover - network/runtime issue
            return False, str(exc)
    return False, last_error


def resolve_ncbi_uid(accession: str) -> str | None:
    """Resolve a WGS master accession (no direct sequence) to its Entrez UID via esearch."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db=nuccore&term={accession}&retmode=json"
    )
    try:
        with urllib.request.urlopen(url) as resp:
            data = json.loads(resp.read())
        ids = data.get("esearchresult", {}).get("idlist", [])
        return ids[0] if ids else None
    except Exception:
        return None


def download_ncbi_sequence(accession: str, out_dir: Path) -> dict[str, str]:
    accession = accession.strip()

    def efetch_url(id_value: str, rettype: str) -> str:
        return (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=nuccore&id={id_value}&rettype={rettype}&retmode=text"
        )

    gbk_file = out_dir / f"{safe_filename(accession)}.gbk"

    gbk_ok, gbk_err = download_url(efetch_url(accession, "gb"), gbk_file)
    time.sleep(NCBI_DELAY_SECONDS)

    # Some accessions (e.g. WGS master records for multi-contig assemblies) aren't
    # resolvable by accession directly but resolve fine once translated to a UID.
    if not gbk_ok:
        uid = resolve_ncbi_uid(accession)
        time.sleep(NCBI_DELAY_SECONDS)
        if uid:
            gbk_ok, gbk_err = download_url(efetch_url(uid, "gb"), gbk_file)
            time.sleep(NCBI_DELAY_SECONDS)

    return {
        "gbk_path": str(gbk_file) if gbk_ok else "",
        "gbk_error": gbk_err or "",
    }


def main() -> int:
    manifest = collect_verified_bgc_ms2_genome_manifest(PODP_DIR)

    genome_cache: dict[tuple[str, str], dict[str, str]] = {}
    link_rows: list[dict[str, str]] = []

    for project in manifest:
        project_id = project["project_id"]
        project_dir = OUT_ROOT / safe_filename(project_id)
        genes_dir = project_dir / "genes"

        for link in project["verified_bgc_ms2_genome_links"]:
            ms2_url = link.get("MS2_URL") or ""
            ms2_filename = safe_filename(Path(ms2_url.split("?")[0]).name) if ms2_url else ""
            ms2_local_path = str(project_dir / ms2_filename) if ms2_filename else ""

            accession = link.get("genome_accession")
            genome_label = link.get("genome_label") or ""
            mibig_number = (link.get("BGC_ID") or {}).get("MIBiG_number") or ""

            gbk_path = gbk_error = ""
            if accession:
                cache_key = (project_id, accession)
                if cache_key not in genome_cache:
                    print(f"Fetching genome {accession} for project {project_id} ...")
                    genome_cache[cache_key] = download_ncbi_sequence(accession, genes_dir)
                result = genome_cache[cache_key]
                gbk_path = result["gbk_path"]
                gbk_error = result["gbk_error"]

            link_rows.append({
                "project_id": project_id,
                "mibig_number": mibig_number,
                "link_type": link.get("link") or "",
                "ms2_filename": ms2_filename,
                "ms2_local_path": ms2_local_path,
                "genome_label": genome_label,
                "genome_accession": accession or "",
                "genome_gbk_path": gbk_path,
                "gbk_error": gbk_error,
            })

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    with LINK_MANIFEST_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "project_id", "mibig_number", "link_type", "ms2_filename", "ms2_local_path",
            "genome_label", "genome_accession", "genome_gbk_path", "gbk_error",
        ])
        writer.writeheader()
        writer.writerows(link_rows)

    n_ok = sum(1 for r in link_rows if r["genome_gbk_path"])
    print(json.dumps({
        "unique_genomes_fetched": len(genome_cache),
        "verified_bgc_ms2_links": len(link_rows),
        "links_with_genome_file": n_ok,
        "link_manifest": str(LINK_MANIFEST_PATH),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
