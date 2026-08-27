#!/usr/bin/env python3
"""Re-fetch the genome records that downloaded as stubs or error pages.

Most genome files in the paired-data tree hold no sequence. Two causes:

- an assembly accession (``GCA_``/``GCF_``) was handed to efetch, which does not
  accept them, and the error page was saved with a .gbk suffix;
- a WGS or CON *master* record was fetched, which is a header with no contigs.

Both are recoverable, because the stub or the ground-truth row still carries an
identifier that resolves to a real assembly. Routes tried, in order:

1. the accession is itself an assembly            -> NCBI assembly FTP
2. the stub names one in ``DBLINK  Assembly:``    -> NCBI assembly FTP
3. a BioSample accession is known                 -> esearch -> assembly -> FTP
4. a WGS prefix is known                          -> esearch -> assembly -> FTP
5. the record is a CON contig-join                -> efetch rettype=gbwithparts
6. a bare JGI IMG genome id, or nothing else      -> search by strain name

Route 6 matches on the organism/strain label rather than an accession, so it is
weaker evidence than the rest and is reported separately for review: a label
without a strain qualifier can match a sibling strain.

Downloads land in a staging folder and are only promoted into the paired-data
tree once they parse as GenBank *and* contain sequence, so a bad fetch can never
overwrite a good file.

    python scripts/refetch_genomes.py --plan
    python scripts/refetch_genomes.py --download
    python scripts/refetch_genomes.py --promote
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GROUND_TRUTH = ROOT / "data" / "PoDP" / "ground_truth.csv"
PAIRED_DATA = ROOT / "data" / "PoDP" / "ground_truth_paired_data"
STAGING = ROOT / "data" / "PoDP" / "genome_refetch"
PLAN_CSV = STAGING / "refetch_plan.csv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
#: NCBI allows 3 requests/second without an API key; stay comfortably under it.
THROTTLE = 0.4
_last_call = 0.0


def _fetch(url: str, tries: int = 4) -> bytes:
    """GET with throttling and backoff on NCBI's 429."""
    global _last_call
    for attempt in range(tries):
        wait = THROTTLE - (time.monotonic() - _last_call)
        if wait > 0:
            time.sleep(wait)
        _last_call = time.monotonic()
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503) and attempt < tries - 1:
                time.sleep(2 ** attempt)
                continue
            raise
    raise RuntimeError("unreachable")


def esearch_assembly(term: str) -> tuple[str, str, str] | None:
    """Return (assembly_accession, ftp_path, organism) for a search term, or None.

    Prefers RefSeq (GCF_) over GenBank and the newest version, because a RefSeq
    assembly carries the annotation the .gbff needs to be useful.
    """
    payload = _fetch(f"{EUTILS}/esearch.fcgi?db=assembly"
                     f"&term={urllib.request.quote(term)}&retmode=json")
    uids = json.loads(payload).get("esearchresult", {}).get("idlist", [])
    if not uids:
        return None
    payload = _fetch(f"{EUTILS}/esummary.fcgi?db=assembly&id={','.join(uids)}&retmode=json")
    result = json.loads(payload).get("result", {})

    best = None
    for uid in result.get("uids", []):
        record = result[uid]
        accession = record.get("assemblyaccession", "")
        ftp = record.get("ftppath_refseq") or record.get("ftppath_genbank")
        if not accession or not ftp:
            continue
        rank = (0 if accession.startswith("GCF_") else 1, -_version_of(accession))
        if best is None or rank < best[0]:
            best = (rank, accession, ftp, record.get("organism", ""))
    return (best[1], best[2], best[3]) if best else None


def _version_of(accession: str) -> int:
    match = re.search(r"\.(\d+)$", accession)
    return int(match.group(1)) if match else 0


def assembly_ftp_gbff(ftp_path: str) -> str:
    """Build the _genomic.gbff.gz URL for an assembly FTP directory."""
    https = ftp_path.replace("ftp://", "https://").rstrip("/")
    return f"{https}/{https.rsplit('/', 1)[-1]}_genomic.gbff.gz"


def stub_assembly(path: Path) -> str | None:
    """Assembly accession named in a stub's DBLINK line, if any."""
    if not path.exists():
        return None
    match = re.search(r"Assembly:\s*(GC[AF]_\d+\.\d+)", path.read_text(errors="replace"))
    return match.group(1) if match else None


def stub_wgs_prefix(path: Path) -> str | None:
    if not path.exists():
        return None
    match = re.search(r"^WGS\s+([A-Z]+)\d", path.read_text(errors="replace"), re.M)
    return match.group(1) if match else None


def is_con_record(path: Path) -> bool:
    return path.exists() and bool(
        re.search(r"^CONTIG\s+join", path.read_text(errors="replace"), re.M))


def classify(text: str) -> str:
    """Same vocabulary build_ground_truth.py uses."""
    if not text.startswith("LOCUS"):
        return "not_genbank_error"
    if re.search(r"^ORIGIN", text, flags=re.MULTILINE):
        return "sequence_ok"
    return "metadata_only_no_contigs"


#: An isolate bacterial genome is a few Mb over tens-to-hundreds of contigs.
#: Well past that and the record is a metagenome or a co-assembly, which will
#: pass a "has ORIGIN" check while being useless as one strain's genome.
MAX_ISOLATE_CONTIGS = 5_000
MAX_ISOLATE_BASES = 50_000_000


def classify_file(path: Path) -> str:
    """Classify a record by streaming it.

    Must not be handed a truncated head: ORIGIN sits after the whole FEATURES
    block, which in a real assembly is megabytes in. Judging from the first few
    kilobytes marks every genuine genome as a header-only stub.
    """
    first = True
    with path.open(errors="replace") as handle:
        for line in handle:
            if first:
                if not line.startswith("LOCUS"):
                    return "not_genbank_error"
                first = False
            if line.startswith("ORIGIN"):
                return "sequence_ok"
    return "metadata_only_no_contigs"


def scale_check(path: Path) -> tuple[str, int, int]:
    """Return (flag, contig_count, total_bases) for a downloaded record."""
    contigs = 0
    bases = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("LOCUS"):
                contigs += 1
                match = re.search(r"^LOCUS\s+\S+\s+(\d+) bp", line)
                if match:
                    bases += int(match.group(1))
    if contigs > MAX_ISOLATE_CONTIGS or bases > MAX_ISOLATE_BASES:
        return "metagenome_scale", contigs, bases
    return "isolate_scale", contigs, bases


# --------------------------------------------------------------------------
# planning
# --------------------------------------------------------------------------

def broken_records(rows: list[dict]) -> dict[tuple[str, str], dict]:
    """{(study_id, accession): row} for every genome that is not usable."""
    out: dict[tuple[str, str], dict] = {}
    for row in rows:
        if row["genome_accession"] and row["status_genome"] != "sequence_ok":
            out.setdefault((row["study_id"], row["genome_accession"]), row)
    return out


def plan_one(study: str, accession: str, row: dict) -> dict:
    """Decide how (or whether) this accession can be re-fetched."""
    stub = PAIRED_DATA / study / "genomes" / f"{accession}.gbk"
    entry = {
        "study_id": study, "genome_accession": accession,
        "current_status": row["status_genome"], "genome_label": row["genome_label"],
        "biosample": row["biosample_accession"], "route": "", "resolved_assembly": "",
        "source_url": "", "resolved_organism": "", "strain_verified": "", "note": "",
    }

    if re.match(r"GC[AF]_\d", accession):
        entry["route"], term = "assembly_accession", accession
    elif (linked := stub_assembly(stub)):
        entry["route"], term = "stub_dblink_assembly", linked
    elif row["biosample_accession"]:
        entry["route"], term = "biosample", f"{row['biosample_accession']}[BioSample]"
    elif (prefix := stub_wgs_prefix(stub)):
        entry["route"], term = "wgs_prefix", prefix
    elif is_con_record(stub):
        entry["route"] = "con_gbwithparts"
        entry["source_url"] = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={accession}"
                               "&rettype=gbwithparts&retmode=text")
        return entry
    elif row["genome_label"]:
        # last resort: the strain label. Weaker than any accession-based route --
        # flagged so it can be reviewed rather than trusted silently.
        entry["route"] = "strain_name"
        term = row["genome_label"].split("|")[0].strip()
        entry["note"] = f"matched on strain label {term!r}, not an accession -- review"
    else:
        entry["route"] = "unresolved"
        entry["note"] = "no assembly, BioSample, WGS prefix, CONTIG line or label"
        return entry

    try:
        found = esearch_assembly(term)
    except Exception as exc:                       # network/parse, keep planning
        entry["note"] = f"lookup failed: {type(exc).__name__}"
        return entry
    if not found and entry["route"] != "strain_name" and row["genome_label"]:
        # the accession-based lookup came back empty; fall through to the label
        label = row["genome_label"].split("|")[0].strip()
        try:
            found = esearch_assembly(label)
        except Exception:
            found = None
        if found:
            entry["route"] = "strain_name"
            entry["note"] = (f"{term} found nothing; matched on strain label "
                             f"{label!r} instead -- review")
    if not found:
        entry["note"] = f"no assembly found for {term}"
        entry["route"] = "unresolved"
        return entry
    entry["resolved_assembly"], ftp, entry["resolved_organism"] = found
    entry["source_url"] = assembly_ftp_gbff(ftp)
    return entry


# --------------------------------------------------------------------------
# download / promote
# --------------------------------------------------------------------------

def download_one(entry: dict, staging: Path) -> dict:
    """Fetch to staging and record what actually arrived."""
    target = staging / f"{entry['genome_accession']}.gbk"
    entry = dict(entry, downloaded="", bytes="", new_status="", promoted="")
    if not entry["source_url"]:
        entry["new_status"] = "no_route"
        return entry
    if target.exists() and target.stat().st_size > 0:
        entry.update(downloaded="cached", bytes=target.stat().st_size,
                     new_status=classify_file(target))
        flag, contigs, bases = scale_check(target)
        entry.update(scale=flag, contigs=contigs, bases=bases)
        return entry
    try:
        payload = _fetch(entry["source_url"])
    except Exception as exc:
        entry["new_status"] = f"download_failed: {type(exc).__name__}"
        return entry
    if payload[:2] == b"\x1f\x8b":
        payload = gzip.decompress(payload)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)
    entry.update(downloaded="yes", bytes=len(payload),
                 new_status=classify_file(target))
    flag, contigs, bases = scale_check(target)
    entry.update(scale=flag, contigs=contigs, bases=bases)
    return entry


def promote_one(entry: dict, staging: Path, paired: Path) -> str:
    """Move a validated download into the paired-data tree, keeping a backup."""
    if entry.get("new_status") != "sequence_ok":
        return "skipped_not_sequence"
    if entry.get("strain_verified") == "no":
        return "skipped_wrong_strain"
    if entry.get("scale") == "metagenome_scale":
        # passes "has ORIGIN" but is not one strain's genome; needs a decision
        return "skipped_metagenome_scale"
    source = staging / f"{entry['genome_accession']}.gbk"
    if not source.exists():
        return "skipped_missing"
    destination = paired / entry["study_id"] / "genomes" / f"{entry['genome_accession']}.gbk"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        # keep the stub -- it is the evidence of what the original fetch returned
        shutil.move(str(destination), str(destination.with_suffix(".gbk.stub")))
    shutil.copy2(source, destination)
    return "promoted"


def _strain(label: str) -> str:
    """Normalised strain label, for deciding whether two rows mean one organism."""
    return re.sub(r"[^a-z0-9]+", "", (label or "").split("|")[0].lower())


def verify_strain(plan: list[dict]) -> int:
    """Check that a label-matched assembly really is the strain that was asked for.

    NCBI's assembly search treats a strain token as a loose term: asking for
    "Salinispora arenicola CNH643" happily returns CNS820. The species matches,
    the strain does not, and the result looks entirely plausible. So for the
    label route, require the strain token to appear in the assembly's own
    organism string; anything else is rejected rather than promoted.
    """
    rejected = 0
    for entry in plan:
        if entry.get("route") != "strain_name":
            entry.setdefault("strain_verified", "n/a")
            continue
        label = (entry.get("genome_label") or "").split("|")[0].strip()
        token = label.split()[-1] if label else ""
        organism = entry.get("resolved_organism") or ""
        if token and token.lower() in organism.lower():
            entry["strain_verified"] = "yes"
        else:
            entry["strain_verified"] = "no"
            entry["note"] = ((entry.get("note") or "") +
                             f" || WRONG STRAIN: asked for {token!r}, assembly is "
                             f"{organism.split('(')[0].strip()!r} -- not promoted").strip()
            rejected += 1
    return rejected


def flag_collisions(plan: list[dict]) -> int:
    """Flag distinct source records that resolved to the *same* assembly.

    Two records sharing an assembly is only a problem when they name different
    strains -- one assembly cannot be two organisms, so at most one match is
    right. The same strain cited two ways (a JGI id in one project, the assembly
    accession in another) legitimately resolves to one assembly and is left
    alone. Conflicts are flagged, never dropped: which is correct needs a human.
    """
    seen: dict[str, list[dict]] = {}
    for entry in plan:
        assembly = entry.get("resolved_assembly")
        if assembly:
            seen.setdefault(assembly, []).append(entry)

    flagged = 0
    for entry in plan:
        peers = seen.get(entry.get("resolved_assembly") or "", [])
        conflicting = [e["genome_accession"] for e in peers
                       if e["genome_accession"] != entry["genome_accession"]
                       and _strain(e.get("genome_label", ""))
                       != _strain(entry.get("genome_label", ""))]
        if conflicting:
            entry["collision"] = "|".join(conflicting)
            entry["note"] = ((entry.get("note") or "") +
                             f" || CONFLICT: assembly also matched by "
                             f"{', '.join(conflicting)}, which names a different "
                             f"strain -- at most one is correct").strip()
            flagged += 1
        else:
            entry["collision"] = ""
    return flagged


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarise(rows: list[dict], key: str, title: str) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.get(key) or "(blank)"] = counts.get(row.get(key) or "(blank)", 0) + 1
    print(f"\n{title}")
    for value, count in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"   {value:<44} {count:>4}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan", action="store_true", help="resolve routes only")
    parser.add_argument("--download", action="store_true", help="fetch into staging")
    parser.add_argument("--promote", action="store_true",
                        help="copy validated downloads into the paired-data tree")
    parser.add_argument("--limit", type=int, help="stop after N records (for testing)")
    parser.add_argument("--ground-truth", type=Path, default=GROUND_TRUTH)
    parser.add_argument("--staging", type=Path, default=STAGING)
    args = parser.parse_args()

    if not any((args.plan, args.download, args.promote)):
        parser.error("choose --plan, --download or --promote")

    with args.ground_truth.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    broken = broken_records(rows)
    print(f"genome records needing attention: {len(broken)}")

    if args.plan or not PLAN_CSV.exists():
        plan: list[dict] = []
        for index, ((study, accession), row) in enumerate(sorted(broken.items()), 1):
            if args.limit and index > args.limit:
                break
            entry = plan_one(study, accession, row)
            plan.append(entry)
            print(f"  [{index:>2}/{len(broken)}] {accession:<22} {entry['route']:<22}"
                  f"{entry['resolved_assembly'] or entry['note'][:40]}")
        rejected = verify_strain(plan)
        if rejected:
            print(f"\n!! {rejected} label-matched record(s) are the WRONG STRAIN "
                  f"and will not be promoted")
        flagged = flag_collisions(plan)
        write_csv(PLAN_CSV, plan)
        print(f"\nplan written: {PLAN_CSV}")
        summarise(plan, "route", "routes resolved")
        if flagged:
            print(f"\n!! {flagged} record(s) share a resolved assembly with another "
                  f"-- at most one of each pair is correct:")
            for entry in plan:
                if entry.get("collision"):
                    print(f"   {entry['genome_accession']:<20} {entry['resolved_assembly']:<18}"
                          f"also claimed by {entry['collision']}")
    else:
        with PLAN_CSV.open(newline="", encoding="utf-8") as handle:
            plan = list(csv.DictReader(handle))

    if args.download:
        results = []
        for index, entry in enumerate(plan, 1):
            result = download_one(entry, args.staging)
            results.append(result)
            print(f"  [{index:>2}/{len(plan)}] {result['genome_accession']:<22}"
                  f"{result['new_status']:<28}{result['bytes']}")
        write_csv(PLAN_CSV, results)
        summarise(results, "new_status", "download outcome")
        summarise(results, "scale", "record scale")
        plan = results

    if args.promote:
        for entry in plan:
            entry["promoted"] = promote_one(entry, args.staging, PAIRED_DATA)
        write_csv(PLAN_CSV, plan)
        summarise(plan, "promoted", "promotion outcome")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
