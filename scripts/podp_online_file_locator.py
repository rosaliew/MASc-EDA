#!/usr/bin/env python3
"""Map Paired omics Data Platform JSON metadata to online resource URLs.

The PoDP project JSON files do not contain the raw sequence and MS2 files directly,
but they do contain the identifiers needed to locate them online:

- metabolomics.project.GNPSMassIVE_ID / MaSSIVE_URL
- genomes[*].genome_ID.*_accession or JGI_Genome_ID

This module normalizes those identifiers into reusable URLs and exposes a small
CLI that scans a PoDP directory.
"""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from pathlib import Path
from typing import Any


def normalize_identifier(value: Any) -> list[str]:
    """Split a raw metadata value into a list of clean identifiers."""
    if value is None:
        return []

    if isinstance(value, (int, float)):
        value = str(value)

    if not isinstance(value, str):
        return []

    raw = value.strip()
    if not raw:
        return []

    parts = re.split(r"[,;\s]+", raw)
    cleaned: list[str] = []
    for part in parts:
        part = part.strip().strip("()[]{}'")
        if part and not part.startswith("http"):
            cleaned.append(part.rstrip("."))
    return cleaned


def _coerce_accession(value: Any) -> list[str]:
    return normalize_identifier(value)


def make_ncbi_urls(accession: str) -> dict[str, str]:
    """Return a few useful NCBI query and download endpoints for an accession."""
    accession = accession.strip()
    polished = accession.rstrip(".")
    urls = {
        "ncbi_assembly_url": f"https://www.ncbi.nlm.nih.gov/assembly/?term={polished}",
        "ncbi_nuccore_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{polished}",
        "ncbi_fasta_url": (
            "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?"
            f"db=nuccore&id={polished}&rettype=fasta&retmode=text"
        ),
    }
    return urls


def make_jgi_urls(jgi_id: str) -> dict[str, str]:
    """Return basic JGI landing-page URLs for a JGI genome identifier."""
    jgi_id = jgi_id.strip()
    return {
        "jgi_genome_url": f"https://genome.jgi.doe.gov/portal/{jgi_id}.html",
        "jgi_download_url": f"https://genome.jgi.doe.gov/portal/{jgi_id}/download/",
    }


def is_ncbi_accession(value: str) -> bool:
    """Return True for the NCBI accession patterns used by PoDP genome records."""
    accession = value.strip()
    if not accession:
        return False
    if accession.startswith(("GCA_", "GCF_", "GCF", "GCA")):
        return True
    if re.fullmatch(r"[A-Z]{1,2}[A-Z0-9]*\d+[A-Z0-9]*\.?\d*", accession, flags=re.IGNORECASE):
        return True
    if re.fullmatch(r"NZ_[A-Z0-9]+\.\d+", accession, flags=re.IGNORECASE):
        return True
    return False


def find_metabolomics_files(payload: dict[str, Any]) -> list[str]:
    """Recursively in the project JSON find all known metabolomics file URIs."""
    files: list[str] = []
    seen: set[str] = set()

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "metabolomics_file" and isinstance(value, str):
                    url = value.strip()
                    if url and url not in seen:
                        files.append(url)
                        seen.add(url)
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(payload)
    return files


def extract_project_links(project_json: dict[str, Any]) -> dict[str, Any]:
    """Extract MassIVE and genome identifiers and translate them into online links."""
    metabolomics = project_json.get("metabolomics", {}) or {}
    project = metabolomics.get("project", {}) or {}

    massive_id = project.get("GNPSMassIVE_ID") or project.get("GNPSMassIVE_ID ")
    massive_url = project.get("MaSSIVE_URL") or project.get("MassIVE_URL")

    if massive_id:
        massive_id = str(massive_id).strip()
    if not massive_url and massive_id:
        massive_url = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task={massive_id}"

    genome_entries = project_json.get("genomes", []) or []
    genome_links: list[dict[str, Any]] = []

    for genome in genome_entries:
        if not isinstance(genome, dict):
            continue

        genome_id = genome.get("genome_ID", {}) or {}
        label = genome.get("genome_label") or genome.get("genome_name")
        genome_type = genome_id.get("genome_type") or genome.get("genome_type")

        identifiers: list[str] = []
        for key in (
            "GenBank_accession",
            "RefSeq_accession",
            "ENA_NCBI_accession",
            "NCBI_accession",
            "JGI_Genome_ID",
        ):
            identifiers.extend(_coerce_accession(genome_id.get(key)))

        # Some PoDP entries keep the accession as a top-level field, not nested
        for key in ("GenBank_accession", "RefSeq_accession", "ENA_NCBI_accession", "JGI_Genome_ID"):
            identifiers.extend(_coerce_accession(genome.get(key)))

        unique_identifiers = []
        seen_ids: set[str] = set()
        for item in identifiers:
            if item and item not in seen_ids:
                unique_identifiers.append(item)
                seen_ids.add(item)

        genome_entry: dict[str, Any] = {
            "genome_label": label,
            "genome_type": genome_type,
            "accessions": unique_identifiers,
            "genome_accession": unique_identifiers[0] if unique_identifiers else None,
        }

        if unique_identifiers:
            first = unique_identifiers[0]
            if is_ncbi_accession(first):
                genome_entry.update(make_ncbi_urls(first))
            elif re.fullmatch(r"[A-Za-z0-9]+", first) and len(first) >= 5:
                genome_entry["ncbi_assembly_url"] = f"https://www.ncbi.nlm.nih.gov/assembly/?term={first}"
                genome_entry["ncbi_nuccore_url"] = f"https://www.ncbi.nlm.nih.gov/nuccore/{first}"

        jgi_ids = [
            item for item in unique_identifiers if item.isdigit() and len(item) >= 6
        ]
        if jgi_ids:
            genome_entry.update(make_jgi_urls(jgi_ids[0]))

        genome_links.append(genome_entry)

    result = {
        "massive_id": massive_id,
        "massive_url": massive_url,
        "massive_files": find_metabolomics_files(project_json),
        "genomes": genome_links,
    }
    return result


def build_download_plan(project_json: dict[str, Any], project_id: str | None = None, output_dir: str | Path | None = None) -> dict[str, Any]:
    """Build a list of FASTA/GBK and MS2 download actions for a PoDP project."""
    links = extract_project_links(project_json)
    if project_id is None:
        project_id = project_json.get("project_id") or "PoDP_project"

    base_dir = Path(output_dir) if output_dir is not None else Path("data/downloads") / project_id
    genome_requests: list[dict[str, str]] = []
    for genome in links.get("genomes", []):
        accessions = genome.get("accessions") or [genome.get("genome_accession")]
        for accession in accessions:
            if not accession:
                continue
            safe_accession = str(accession).strip().replace("/", "_").replace(" ", "_")
            genome_dir = base_dir / "genomes" / safe_accession
            genome_requests.append({
                "accession": str(accession).strip(),
                "kind": "fasta",
                "url": f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text",
                "output_path": str(genome_dir / f"{safe_accession}.fasta"),
            })
            genome_requests.append({
                "accession": str(accession).strip(),
                "kind": "gbk",
                "url": f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gb&retmode=text",
                "output_path": str(genome_dir / f"{safe_accession}.gbk"),
            })

    ms2_requests: list[dict[str, str]] = []
    for file_url in links.get("massive_files", []):
        url = str(file_url).strip()
        file_name = Path(url.split("?")[0]).name
        dest_path = base_dir / "ms2" / file_name
        ms2_requests.append({
            "kind": "ms2",
            "url": url,
            "output_path": str(dest_path),
        })

    return {
        "project_id": project_id,
        "output_dir": str(base_dir),
        "massive_id": links.get("massive_id"),
        "massive_url": links.get("massive_url"),
        "genome_requests": genome_requests,
        "ms2_requests": ms2_requests,
    }


def download_url(url: str, output_path: str | Path) -> str:
    """Download a URL to a filesystem path, creating parent dirs as needed."""
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, destination)
    return str(destination)


def download_project_data(project_json: dict[str, Any], project_id: str | None = None, output_dir: str | Path | None = None) -> dict[str, Any]:
    """Download the FASTA/GBK and MS2 files for one PoDP project."""
    plan = build_download_plan(project_json, project_id=project_id, output_dir=output_dir)
    downloaded: list[dict[str, str]] = []
    for req in plan["genome_requests"]:
        try:
            file_path = download_url(req["url"], req["output_path"])
            downloaded.append({"kind": req["kind"], "url": req["url"], "output_path": file_path})
        except Exception as exc:  # pragma: no cover - network/download failures are runtime-only
            downloaded.append({"kind": req["kind"], "url": req["url"], "output_path": "FAILED", "error": str(exc)})
    for req in plan["ms2_requests"]:
        try:
            file_path = download_url(req["url"], req["output_path"])
            downloaded.append({"kind": req["kind"], "url": req["url"], "output_path": file_path})
        except Exception as exc:  # pragma: no cover - network/download failures are runtime-only
            downloaded.append({"kind": req["kind"], "url": req["url"], "output_path": "FAILED", "error": str(exc)})
    plan["downloaded_files"] = downloaded
    return plan


def scan_podp_directory(directory: str | Path) -> list[dict[str, Any]]:
    root = Path(directory)
    results: list[dict[str, Any]] = []
    for json_path in sorted(root.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text())
        except json.JSONDecodeError:
            continue
        record = extract_project_links(payload)
        record["project_file"] = json_path.name
        results.append(record)
    return results


def is_confirmed_paired_project(project_json: dict[str, Any]) -> bool:
    """Return True only when a project has both genome IDs and linked MS2 file entries."""
    metabolomics = project_json.get("metabolomics", {}) or {}
    project = metabolomics.get("project", {}) or {}
    massive_id = project.get("GNPSMassIVE_ID") or project.get("GNPSMassIVE_ID ")
    massive_url = project.get("MaSSIVE_URL") or project.get("MassIVE_URL")
    mass_files = find_metabolomics_files(project_json)
    genomes = project_json.get("genomes", []) or []

    has_genome_accession = False
    for genome in genomes:
        if not isinstance(genome, dict):
            continue
        genome_id = genome.get("genome_ID", {}) or {}
        accession_keys = (
            "GenBank_accession",
            "RefSeq_accession",
            "ENA_NCBI_accession",
            "NCBI_accession",
            "JGI_Genome_ID",
        )
        for key in accession_keys:
            val = genome_id.get(key) or genome.get(key)
            if val and normalize_identifier(val):
                has_genome_accession = True
                break
        if has_genome_accession:
            break

    has_massive_link = bool(massive_id or massive_url or mass_files)
    return bool(has_genome_accession and has_massive_link and mass_files)


def collect_confirmed_paired_projects(directory: str | Path) -> list[dict[str, Any]]:
    """Collect JSON projects where genome and MS2 metadata are explicitly linked."""
    root = Path(directory)
    paired: list[dict[str, Any]] = []
    for json_path in sorted(root.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text())
        except json.JSONDecodeError:
            continue
        if not is_confirmed_paired_project(payload):
            continue
        links = extract_project_links(payload)
        paired.append({
            "project_file": json_path.name,
            "project_id": json_path.stem,
            "massive_id": links.get("massive_id"),
            "massive_url": links.get("massive_url"),
            "genome_accessions": [
                g.get("genome_accession") for g in links.get("genomes", []) if g.get("genome_accession")
            ],
            "ms2_file_count": len(links.get("massive_files", [])),
            "ms2_files": links.get("massive_files", []),
        })
    return paired


def write_confirmed_paired_manifest(directory: str | Path, output_path: str | Path) -> list[dict[str, Any]]:
    """Write a JSON manifest for all confirmed paired PoDP projects and return it."""
    paired = collect_confirmed_paired_projects(directory)
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(paired, indent=2, sort_keys=True))
    return paired


def _as_verification_list(value: Any) -> list[str]:
    """Normalize a verification field to a list of strings."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value.strip() else []
    if isinstance(value, list):
        items: list[str] = []
        for item in value:
            if item is None:
                continue
            if isinstance(item, str) and item.strip():
                items.append(item.strip())
        return items
    return []


def is_verified_bgc_ms2_link(link: dict[str, Any]) -> bool:
    """Return True only when a BGC-MS2 link points to a spectrum and has non-empty verification."""
    if not isinstance(link, dict):
        return False
    has_ms2_pointer = bool(link.get("MS2_URL") or link.get("MS2_scan"))
    verification = _as_verification_list(link.get("verification"))
    return has_ms2_pointer and bool(verification)


def collect_verified_bgc_ms2_links(project_json: dict[str, Any]) -> list[dict[str, Any]]:
    """Collect only the BGC-MS2 entries with a real MS2 pointer and a validated evidence note."""
    rows: list[dict[str, Any]] = []
    project_id = project_json.get("_project_id") or project_json.get("project_id") or "unknown"
    for link in project_json.get("BGC_MS2_links", []) or []:
        if not isinstance(link, dict):
            continue
        if not is_verified_bgc_ms2_link(link):
            continue
        row = dict(link)
        row["project_id"] = project_id
        row["verification"] = _as_verification_list(link.get("verification"))
        rows.append(row)
    return rows


def collect_verified_bgc_ms2_manifest(directory: str | Path) -> list[dict[str, Any]]:
    """Collect project-level manifests for verified BGC spectrum links only."""
    root = Path(directory)
    results: list[dict[str, Any]] = []
    for json_path in sorted(root.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text())
        except json.JSONDecodeError:
            continue
        verified_links = collect_verified_bgc_ms2_links(payload)
        if not verified_links:
            continue
        project_id = payload.get("_project_id") or json_path.stem.rsplit(".", 1)[0]
        project_payload = extract_project_links(payload)
        results.append({
            "project_file": json_path.name,
            "project_id": project_id,
            "massive_id": project_payload.get("massive_id"),
            "massive_url": project_payload.get("massive_url"),
            "verified_bgc_ms2_link_count": len(verified_links),
            "verified_bgc_ms2_links": verified_links,
            "ms2_files": [
                link.get("MS2_URL") for link in verified_links if link.get("MS2_URL")
            ],
        })
    return results


def write_verified_bgc_ms2_manifest(directory: str | Path, output_path: str | Path) -> list[dict[str, Any]]:
    """Write a manifest of verified BGC-MS2 links and return it."""
    manifest = collect_verified_bgc_ms2_manifest(directory)
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


_GENOME_ACCESSION_KEYS = (
    "GenBank_accession",
    "RefSeq_accession",
    "ENA_NCBI_accession",
    "NCBI_accession",
    "JGI_Genome_ID",
)


def _genome_accession_for_label(genomes: list[Any], genome_label: str) -> tuple[str, str] | None:
    """Return the first (accession_key, accession) found for a genome_label, or None."""
    for genome in genomes:
        if not isinstance(genome, dict) or genome.get("genome_label") != genome_label:
            continue
        genome_id = genome.get("genome_ID", {}) or {}
        for key in _GENOME_ACCESSION_KEYS:
            ids = normalize_identifier(genome_id.get(key) or genome.get(key))
            if ids:
                return key, ids[0]
    return None


def collect_verified_bgc_ms2_genome_links(project_json: dict[str, Any]) -> list[dict[str, Any]]:
    """Pair each verified BGC-MS2 link with the genome that produced its MS2 spectrum.

    The bridge is `genome_metabolome_links`: each entry there pins one
    `genome_label` to one `metabolomics_file` URL. A verified BGC-MS2 link's
    `MS2_URL` is that same metabolomics_file, so matching on that URL gives an
    unambiguous genome_label, which is then resolved to an NCBI/JGI accession
    via the project's `genomes` block.
    """
    genomes = project_json.get("genomes", []) or []
    gm_links = project_json.get("genome_metabolome_links", []) or []

    file_to_labels: dict[str, set[str]] = {}
    for gm in gm_links:
        if not isinstance(gm, dict):
            continue
        metabolomics_file = gm.get("metabolomics_file")
        genome_label = gm.get("genome_label")
        if metabolomics_file and genome_label:
            file_to_labels.setdefault(metabolomics_file, set()).add(genome_label)

    rows: list[dict[str, Any]] = []
    for link in collect_verified_bgc_ms2_links(project_json):
        ms2_url = link.get("MS2_URL")
        genome_labels = sorted(file_to_labels.get(ms2_url, set())) if ms2_url else []
        for genome_label in genome_labels or [None]:
            accession_info = _genome_accession_for_label(genomes, genome_label) if genome_label else None
            row = dict(link)
            row["genome_label"] = genome_label
            row["genome_accession"] = accession_info[1] if accession_info else None
            row["genome_accession_type"] = accession_info[0] if accession_info else None
            rows.append(row)
    return rows


def collect_verified_bgc_ms2_genome_manifest(directory: str | Path) -> list[dict[str, Any]]:
    """Collect project-level manifests pairing verified BGC-MS2 links to genome accessions."""
    root = Path(directory)
    results: list[dict[str, Any]] = []
    for json_path in sorted(root.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text())
        except json.JSONDecodeError:
            continue
        linked = collect_verified_bgc_ms2_genome_links(payload)
        if not linked:
            continue
        project_id = payload.get("_project_id") or json_path.stem.rsplit(".", 1)[0]
        results.append({
            "project_file": json_path.name,
            "project_id": project_id,
            "verified_bgc_ms2_genome_links": linked,
        })
    return results


def write_verified_bgc_ms2_genome_manifest(directory: str | Path, output_path: str | Path) -> list[dict[str, Any]]:
    """Write a manifest pairing verified BGC-MS2 links to genome accessions and return it."""
    manifest = collect_verified_bgc_ms2_genome_manifest(directory)
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract MassIVE and genome identifiers from PoDP project JSON files."
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default="data/PoDP/json_descriptors",
        help="Path to the directory containing PoDP JSON files (default: data/PoDP/json_descriptors)",
    )
    parser.add_argument(
        "--project",
        help="Optional project JSON file to inspect; if provided, scan only this file.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/downloads",
        help="Base folder for downloaded FASTA/GBK and MS2 outputs (default: data/downloads)",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download the genome FASTA/GBK files and MS2 files into the output folder.",
    )
    parser.add_argument(
        "--paired-only",
        action="store_true",
        help="Only keep projects where a genome accession and MassIVE/MS2 file entries are both confirmed in the JSON.",
    )
    parser.add_argument(
        "--manifest",
        help="Write a JSON manifest for the paired project set to this path.",
    )
    parser.add_argument(
        "--verified-bgc-ms2",
        action="store_true",
        help="Keep only BGC-MS2 links that have an actual MS2 pointer and non-empty verification evidence.",
    )
    parser.add_argument(
        "--verified-bgc-ms2-manifest",
        help="Write a JSON manifest containing only verified BGC-MS2 links to this path.",
    )
    parser.add_argument(
        "--verified-bgc-ms2-genome-manifest",
        help="Write a JSON manifest pairing verified BGC-MS2 links to their genome accessions to this path.",
    )
    args = parser.parse_args()

    if args.project:
        path = Path(args.project)
        payload = json.loads(path.read_text())
        project_id = path.stem
        if args.paired_only and not is_confirmed_paired_project(payload):
            print(json.dumps({"project_id": project_id, "paired": False}, indent=2, sort_keys=True))
            return
        result = extract_project_links(payload)
        if args.download:
            plan = download_project_data(payload, project_id=project_id, output_dir=args.output_dir)
            print(json.dumps({
                "project_id": plan["project_id"],
                "output_dir": plan["output_dir"],
                "genome_requests": len(plan["genome_requests"]),
                "ms2_requests": len(plan["ms2_requests"]),
                "downloaded_files": plan["downloaded_files"][:5],
            }, indent=2, sort_keys=True))
        else:
            print(json.dumps(result, indent=2, sort_keys=True))
        return

    if args.verified_bgc_ms2 or args.verified_bgc_ms2_manifest:
        verified = collect_verified_bgc_ms2_manifest(args.directory)
        if args.verified_bgc_ms2_manifest:
            write_verified_bgc_ms2_manifest(args.directory, args.verified_bgc_ms2_manifest)
        print(json.dumps({"verified_bgc_ms2_project_count": len(verified), "projects": verified}, indent=2, sort_keys=True))
        return

    if args.verified_bgc_ms2_genome_manifest:
        linked = write_verified_bgc_ms2_genome_manifest(args.directory, args.verified_bgc_ms2_genome_manifest)
        print(json.dumps({"verified_bgc_ms2_genome_project_count": len(linked), "projects": linked}, indent=2, sort_keys=True))
        return

    if args.paired_only:
        paired = collect_confirmed_paired_projects(args.directory)
        if args.manifest:
            write_confirmed_paired_manifest(args.directory, args.manifest)
        print(json.dumps({"paired_project_count": len(paired), "projects": paired}, indent=2, sort_keys=True))
        return

    for record in scan_podp_directory(args.directory):
        print(f"Project: {record['project_file']}")
        print(json.dumps({
            "massive_id": record.get("massive_id"),
            "massive_url": record.get("massive_url"),
            "massive_files_count": len(record.get("massive_files", [])),
            "genome_count": len(record.get("genomes", [])),
        }, indent=2, sort_keys=True))
        if args.download:
            project_path = Path(args.directory) / record["project_file"]
            payload = json.loads(project_path.read_text())
            plan = download_project_data(payload, project_id=project_path.stem, output_dir=args.output_dir)
            print(json.dumps({
                "output_dir": plan["output_dir"],
                "genome_requests": len(plan["genome_requests"]),
                "ms2_requests": len(plan["ms2_requests"]),
            }, indent=2, sort_keys=True))
        else:
            for genome in record.get("genomes", [])[:3]:
                print(json.dumps(genome, indent=2, sort_keys=True))
        print()


if __name__ == "__main__":
    main()
