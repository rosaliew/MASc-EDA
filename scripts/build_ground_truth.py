#!/usr/bin/env python3
"""Build tiered ground-truth tables from the PoDP JSON descriptors.

Emits three CSVs into data/PoDP/:

- ``ground_truth.csv``    one row per BGC-MS2 link, BGC <-> spectrum. Tiered.
- ``strain_pairings.csv`` one row per genome-metabolome link, strain <-> file.
                          NOT ground truth and deliberately untiered -- these are
                          declared co-origin pairings, not evidence of any
                          BGC <-> spectrum link. Kept as a separate sheet so they
                          are never mistaken for scoreable ground truth.
The funnel is printed to the console as a build log and rendered live in
``notebooks/podp-2026-08-ground-truth-report.ipynb``; it is not written to disk,
so there is no stale copy to drift from the CSVs.

Tiers describe *biological evidence only*. Data completeness, file status and
Metcalf rankability are separate columns, so a row can be tier 1 and still be
unusable, and both facts stay visible.

    tier 1  verification contains "Experimentally validated ..."
    tier 2  canonical InChIKey match to the cited MIBiG compound -- either the
            full key, or the connectivity layer alone (stereo differs)
    tier 3  neither -- no experimental validation and no structural match
            (kept, with tier_reason, never dropped)

Tier 2 was originally split into whole-SMILES and partial-SMILES matches, but
every whole-SMILES match landed on a link that was already experimentally
validated, so the whole-SMILES tier was always empty and the two were merged.

MIBiG retirement does NOT affect tier. Tier grades the biological claim -- that
this BGC produces this compound -- while retirement grades MIBiG's own record
("missing gene annotations", "spread over multiple contigs"). A knockout
experiment is not undone by a database annotation gap, so retirement is carried
in mibig_status / mibig_quality / mibig_retirement_reason instead, the same way
rankability and file status are kept off the tier axis.

Tier 2 is decided by rdkit's canonical InChIKey, so rdkit is required. A SMILES
string screen is still computed as a diagnostic and reported against the InChIKey
verdict, but it never assigns a tier: its strong calls (identical / permutation)
are exact, while the atom-formula screen is wrong about a third of the time and a
0.60 similarity threshold was refuted on every row it fired on.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.podp_online_file_locator import (  # noqa: E402
    _as_verification_list,
    _genome_accession_for_label,
    normalize_identifier,
)

ROOT = Path(__file__).resolve().parent.parent
PODP_DIR = ROOT / "data" / "PoDP" / "json_descriptors"
DOWNLOAD_DIR = ROOT / "data" / "PoDP" / "ground_truth_paired_data"
MIBIG_VERSION = "4.0"
#: Gitignored and re-downloaded on demand -- a released reference set is a derived
#: artifact reproducible from one URL, not project data, so deleting it must never
#: break the build and it does not belong in version control.
MIBIG_CACHE = ROOT / "data"
MIBIG_DIR = MIBIG_CACHE / f"MiBIG_{MIBIG_VERSION}"
#: MIBiG 3.1, as shipped by NPLinker's own downloader. Kept only so
#: --migration-report can diff the two releases.
MIBIG_31_DIR = (
    Path("/home/rosaliewang/projects/NPLinker-Application")
    / "nplinker-dev" / "podp_exercise" / "mibig"
)
OUT_DIR = ROOT / "data" / "PoDP"

#: SequenceMatcher ratio above which the diagnostic screen reports "partial".
#: Diagnostic only -- it no longer assigns a tier. Every row it fired on was
#: refuted by rdkit, which is why InChIKey is now authoritative.
PARTIAL_THRESHOLD = 0.60

#: Declared-genome count below which Metcalf cannot rank. This is about the
#: strain axis of the co-occurrence statistic, NOT about antiSMASH -- antiSMASH
#: runs fine on one genome. With a single strain the occurrence matrix has one
#: column, every pair scores identically and nothing is rankable.
METCALF_MIN_GENOMES = 5

FILE_KINDS = ("bgc", "genome", "ms2_run", "ms2_spectrum")

#: Nothing extracts single MS2 spectra to their own file yet -- the spectrum is
#: addressed as a scan number inside the run file. So ms2_spectrum is missing on
#: every row, and a plain complete/incomplete split would collapse to one value.
#: Rows missing *only* that get their own status, to keep the column useful.
PENDING_KINDS = ("ms2_spectrum",)


def summarise_status(missing: list[str]) -> str:
    """complete | spectrum_pending | incomplete."""
    if not missing:
        return "complete"
    if set(missing) <= set(PENDING_KINDS):
        return "spectrum_pending"
    return "incomplete"


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

def load_descriptors(podp_dir: Path) -> tuple[dict[str, dict], list[str]]:
    """Return {filename: payload} for every parseable descriptor, plus skips."""
    payloads: dict[str, dict] = {}
    skipped: list[str] = []
    for path in sorted(podp_dir.glob("*")):
        if not path.is_file():
            continue
        if path.suffix != ".json":
            skipped.append(f"{path.name} (not .json)")
            continue
        try:
            payloads[path.name] = json.loads(path.read_text())
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            skipped.append(f"{path.name} ({type(exc).__name__})")
    return payloads, skipped


def _index_mibig_31(cluster: dict) -> dict:
    """Normalise a MIBiG <= 3.1 record (everything nested under 'cluster')."""
    loci = cluster.get("loci") or {}
    return {
        "schema": "3.1",
        "compounds": [
            {
                "name": c.get("compound"),
                "smiles": c.get("chem_struct"),
                "formula": c.get("molecular_formula"),
            }
            for c in (cluster.get("compounds") or [])
        ],
        "locus_accession": loci.get("accession"),
        "locus_start": loci.get("start_coord"),
        "locus_end": loci.get("end_coord"),
        "organism": cluster.get("organism_name"),
        "status": cluster.get("status"),
        "quality": "",
        "completeness": loci.get("completeness"),
        "retirement_reasons": [],
    }


def _index_mibig_40(record: dict) -> dict:
    """Normalise a MIBiG 4.0 record (flat, and loci is a list of locations)."""
    loci = record.get("loci") or []
    first = loci[0] if loci else {}
    location = first.get("location") or {}
    return {
        "schema": "4.0",
        "compounds": [
            {
                "name": c.get("name"),
                "smiles": c.get("structure"),
                "formula": c.get("formula"),
            }
            for c in (record.get("compounds") or [])
        ],
        "locus_accession": first.get("accession"),
        "locus_start": location.get("from"),
        "locus_end": location.get("to"),
        "organism": (record.get("taxonomy") or {}).get("name"),
        "status": record.get("status"),
        "quality": record.get("quality"),
        "completeness": record.get("completeness"),
        # 4.0 ships retired entries rather than deleting them, and says why
        "retirement_reasons": record.get("retirement_reasons") or [],
    }


def ensure_mibig(version: str, cache: Path) -> Path:
    """Return the MIBiG JSON directory, downloading the release if it is absent.

    Same bundle URL NPLinker itself uses. Kept out of git: a 15 MB reference set
    is reproducible from one URL, so it is cached rather than vendored.
    """
    target = cache / f"MiBIG_{version}"
    if target.is_dir() and any(target.glob("BGC*.json")):
        return target

    import io
    import tarfile
    import urllib.request

    url = f"https://dl.secondarymetabolites.org/mibig/mibig_json_{version}.tar.gz"
    print(f"MIBiG {version} not cached -- downloading {url}")
    cache.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=300) as response:
        payload = response.read()
    with tarfile.open(fileobj=io.BytesIO(payload)) as archive:
        archive.extractall(cache)

    # the bundle extracts as mibig_json_<version>/ (or a bare pile of files);
    # normalise either onto the folder name this project uses
    extracted = cache / f"mibig_json_{version}"
    if extracted.is_dir() and not target.is_dir():
        extracted.rename(target)
    if not target.is_dir():
        target.mkdir(parents=True, exist_ok=True)
        for path in cache.glob("BGC*.json"):
            path.rename(target / path.name)
    print(f"MIBiG {version} cached at {target}")
    return target


def build_mibig_index(mibig_dir: Path) -> dict[str, dict]:
    """Return {BGC accession: normalised record}, reading 3.1 or 4.0 layout.

    The two releases are structurally unrelated: 3.1 nests everything under a
    'cluster' key and names compound fields compound/chem_struct/molecular_formula,
    while 4.0 is flat and uses name/structure/formula, turns loci into a list, and
    adds status/quality/retirement_reasons. Detect per file rather than per
    directory so a mixed directory cannot silently half-parse.
    """
    index: dict[str, dict] = {}
    if not mibig_dir.is_dir():
        return index
    for path in sorted(mibig_dir.glob("BGC*.json")):
        try:
            payload = json.loads(path.read_text())
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        if "cluster" in payload:
            entry = _index_mibig_31(payload["cluster"])
        else:
            entry = _index_mibig_40(payload)
        entry["path"] = path
        index[path.stem] = entry
    return index


_MASSIVE_RE = re.compile(r"(MSV\d+)")


def massive_accessions(payload: dict) -> frozenset[str]:
    """Every MassIVE accession a project refers to, declared or embedded in a URL."""
    found: set[str] = set()
    metabolomics = payload.get("metabolomics") or {}
    project = metabolomics.get("project") or {}
    for key in ("GNPSMassIVE_ID",):
        if project.get(key):
            found.add(str(project[key]).strip())
    if metabolomics.get("related_GNPSMassIVE_ID"):
        found.add(str(metabolomics["related_GNPSMassIVE_ID"]).strip())
    for link in (payload.get("genome_metabolome_links") or []):
        match = _MASSIVE_RE.search(str(link.get("metabolomics_file") or ""))
        if match:
            found.add(match.group(1))
    for link in (payload.get("BGC_MS2_links") or []):
        match = _MASSIVE_RE.search(str(link.get("MS2_URL") or ""))
        if match:
            found.add(match.group(1))
    return frozenset(found)


def resolve_studies(payloads: dict[str, dict]) -> dict[str, dict]:
    """Group projects that reference exactly the same MassIVE datasets.

    PoDP sometimes carries one study as two project records -- same PI, same
    depositions, different accession. When the MassIVE sets match *entirely*
    those records describe one body of data, so they share one folder on disk
    named for the main study. A partial overlap is left alone: those projects
    genuinely differ.

    Main study = most ground-truth links, then most strain pairings, then most
    declared genomes, then lexicographically smallest podp_id so the choice is
    deterministic across runs.

    Returns {podp_id: {study_id, study_members, study_role}}.
    """
    by_massive: dict[frozenset[str], list[str]] = {}
    for filename, payload in payloads.items():
        accessions = massive_accessions(payload)
        if not accessions:
            continue
        podp_id = filename[:-5] if filename.endswith(".json") else filename
        by_massive.setdefault(accessions, []).append(podp_id)

    def weight(podp_id: str) -> tuple:
        payload = payloads.get(f"{podp_id}.json", {})
        return (
            -len(payload.get("BGC_MS2_links") or []),
            -len(payload.get("genome_metabolome_links") or []),
            -len(payload.get("genomes") or []),
            podp_id,
        )

    studies: dict[str, dict] = {}
    for members in by_massive.values():
        main = sorted(members, key=weight)[0]
        for podp_id in members:
            studies[podp_id] = {
                "study_id": main,
                "study_members": "|".join(sorted(members)),
                "study_role": "main" if podp_id == main else "merged",
            }
    return studies


_DUPLICATE_RE = re.compile(r"[Dd]uplicate of (BGC\d+)")


def build_redirect_map(index: dict[str, dict]) -> dict[str, str]:
    """{retired accession: active accession} for MIBiG's duplicate retirements.

    MIBiG 4.0 retires one member of a duplicate pair and names the survivor in
    retirement_reasons ("Duplicate of BGC0001230"). Following that redirect
    recovers an active, usually richer record. Only followed when the target
    exists and is itself active -- a chain of retirements is left alone.
    """
    redirects: dict[str, str] = {}
    for accession, entry in index.items():
        for reason in entry.get("retirement_reasons") or []:
            match = _DUPLICATE_RE.search(str(reason))
            if not match:
                continue
            target = match.group(1)
            if index.get(target, {}).get("status") == "active":
                redirects[accession] = target
    return redirects


# --------------------------------------------------------------------------
# SMILES comparison (string-based; rdkit only verifies, later)
# --------------------------------------------------------------------------

_ATOM_RE = re.compile(r"Cl|Br|[BCNOFPSIbcnops]")
_AROMATIC = {"b": "B", "c": "C", "n": "N", "o": "O", "p": "P", "s": "S"}


def smiles_atom_formula(smiles: str | None) -> tuple | None:
    """Permutation-invariant atom multiset parsed straight from a SMILES string.

    Deliberately crude -- it counts element tokens and folds aromatic lowercase
    onto uppercase. It cannot see implicit hydrogens, so it is a *screen*, not a
    formula. Good enough to say "same heavy atoms", which is what tier 2 needs.
    """
    if not smiles:
        return None
    tokens = _ATOM_RE.findall(smiles)
    if not tokens:
        return None
    counts = Counter(_AROMATIC.get(t, t) for t in tokens)
    return tuple(sorted(counts.items()))


def smiles_char_multiset(smiles: str | None) -> tuple | None:
    """Character multiset -- catches the same string written in another order."""
    if not smiles:
        return None
    return tuple(sorted(Counter(smiles).items()))


_INCHIKEY_CACHE: dict[str, str | None] = {}


def inchikey(smiles: str | None) -> str | None:
    """Canonical InChIKey for a SMILES, or None if rdkit cannot parse it."""
    if not smiles:
        return None
    if smiles in _INCHIKEY_CACHE:
        return _INCHIKEY_CACHE[smiles]
    try:
        from rdkit import Chem
    except ImportError:  # pragma: no cover - guarded at startup
        return None
    mol = Chem.MolFromSmiles(smiles)
    key = Chem.MolToInchiKey(mol) if mol is not None else None
    _INCHIKEY_CACHE[smiles] = key
    return key


def _inchikey_verdict(a: str | None, b: str | None) -> str:
    """full (incl. stereo) | skeleton (connectivity only) | none | uncomputable."""
    if not a or not b:
        return "uncomputable"
    if a == b:
        return "full"
    if a.split("-")[0] == b.split("-")[0]:
        return "skeleton"
    return "none"


def compare_structures(podp: str | None, candidates: list[dict]) -> dict:
    """Compare a PoDP SMILES against every compound of the cited MIBiG entry.

    The InChIKey verdict is authoritative and drives the tier. The string screen
    (identical / permutation / formula / partial) is still computed, but only as a
    diagnostic: it is exact on its strong calls and unreliable on its weak ones --
    the atom-formula screen is wrong about a third of the time and the similarity
    threshold was refuted on every row it fired on. It is kept so screen-vs-rdkit
    disagreement stays visible, never to assign a tier.
    """
    usable = [c for c in candidates if c.get("smiles")]
    blank = {
        "inchikey_match": "uncomputable", "inchikey_podp": "", "inchikey_mibig": "",
        "match": "uncomputable", "similarity": None, "compound": None, "smiles": None,
    }
    if not podp or not usable:
        return blank

    podp_key = inchikey(podp)
    podp_chars = smiles_char_multiset(podp)
    podp_formula = smiles_atom_formula(podp)

    rank = ["uncomputable", "none", "skeleton", "full"]
    screen_rank = ["none", "partial", "formula", "permutation", "identical"]
    best = dict(blank, inchikey_match="none", match="none", similarity=0.0,
                inchikey_podp=podp_key or "")

    for cand in usable:
        other = cand["smiles"]
        verdict = _inchikey_verdict(podp_key, inchikey(other))

        ratio = SequenceMatcher(None, podp, other).ratio()
        if podp == other:
            screen = "identical"
        elif podp_chars == smiles_char_multiset(other):
            screen = "permutation"
        elif podp_formula is not None and podp_formula == smiles_atom_formula(other):
            screen = "formula"
        elif ratio >= PARTIAL_THRESHOLD:
            screen = "partial"
        else:
            screen = "none"

        better = rank.index(verdict) > rank.index(best["inchikey_match"]) or (
            verdict == best["inchikey_match"]
            and screen_rank.index(screen) > screen_rank.index(
                best["match"] if best["match"] in screen_rank else "none")
        )
        if better:
            best = {
                "inchikey_match": verdict,
                "inchikey_podp": podp_key or "",
                "inchikey_mibig": inchikey(other) or "",
                "match": screen,
                "similarity": round(ratio, 4),
                "compound": cand.get("name"),
                "smiles": other,
            }
    return best


# --------------------------------------------------------------------------
# tiering
# --------------------------------------------------------------------------

def is_experimentally_validated(link: dict) -> bool:
    return any(
        str(v).startswith("Experimentally validated")
        for v in _as_verification_list(link.get("verification"))
    )


def assign_tier(link: dict, result: dict, mibig_entry: dict) -> tuple[int, str]:
    """Return (tier, reason). Highest-qualifying tier wins; tiers are exclusive.

    Tier 2 is decided by rdkit's canonical InChIKey, not by string similarity.
    A `skeleton` match (connectivity layer, first block) counts: it means the
    molecular graph agrees and only stereo annotation differs, which PoDP and
    MIBiG curate inconsistently. Which of the two fired is recorded in
    inchikey_match, so the line can be moved without re-tiering.
    """
    if is_experimentally_validated(link):
        tags = [v for v in _as_verification_list(link.get("verification"))
                if str(v).startswith("Experimentally validated")]
        return 1, f"experimentally validated: {'; '.join(tags)}"

    verdict = result["inchikey_match"]
    if verdict == "full":
        return 2, "InChIKey match to cited MIBiG compound (full key, stereo included)"
    if verdict == "skeleton":
        return 2, ("InChIKey match to cited MIBiG compound "
                   "(connectivity layer; stereo differs or is unannotated)")

    if verdict == "uncomputable":
        if not link.get("SMILES"):
            return 3, "not experimentally validated; no SMILES on the PoDP link"
        if not mibig_entry:
            return 3, ("not experimentally validated; cited MIBiG entry absent "
                       "from this MIBiG release")
        if not any(c.get("smiles") for c in mibig_entry.get("compounds", [])):
            return 3, ("not experimentally validated; MIBiG entry carries no "
                       "structure to compare against")
        return 3, "not experimentally validated; rdkit could not parse a structure"
    return 3, "not experimentally validated; InChIKeys differ from cited MIBiG entry"


# --------------------------------------------------------------------------
# on-disk file status
# --------------------------------------------------------------------------

def classify_genome_file(path: Path | None) -> str:
    """Classify a downloaded genome record by what it actually contains.

    A record with no contigs is not a genome. WGS/CON master records carry only
    a header and must be refetched with -style withparts or from the assembly
    FTP; efetch error text saved with a .gbk suffix is not GenBank at all.
    """
    if path is None or not path.exists():
        return "absent"
    try:
        head = path.read_text(errors="replace")
    except OSError:
        return "absent"
    if not head.startswith("LOCUS"):
        return "not_genbank_error"
    if re.search(r"^ORIGIN", head, flags=re.MULTILINE):
        return "sequence_ok"
    return "metadata_only_no_contigs"


def index_downloads(download_dir: Path) -> dict[str, dict[str, Path]]:
    """Index what is on disk, keyed by project directory name.

    Registered under both the directory name and its bare UUID, because the
    tree is named by podp_id (``<uuid>.<version>``) after rename_downloads.py
    but was named by bare UUID before it. Recursive, so it does not care whether
    files sit at the project root or under genomes/ and ms2/.
    """
    index: dict[str, dict[str, Path]] = {}
    if not download_dir.is_dir():
        return index
    for project_dir in sorted(download_dir.iterdir()):
        if not project_dir.is_dir():
            continue
        entry = {path.name: path for path in project_dir.rglob("*") if path.is_file()}
        index[project_dir.name] = entry
        index.setdefault(project_dir.name.split(".")[0], entry)
    return index


def _normalise_filename(name: str) -> str:
    """Fold every run of non-alphanumeric characters to a single underscore.

    The download step renamed files to <MassIVE>_<strain>_<original> and
    sanitised punctuation on the way, so `DetoxinP1-MS2.mzXML` landed as
    `..._DetoxinP1_MS2.mzXML`. Comparing raw stems misses those.
    """
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def _squash(name: str) -> str:
    """Drop every non-alphanumeric character. Used only as an exact-match key."""
    return re.sub(r"[^A-Za-z0-9]+", "", name).lower()


def load_rename_manifest(download_dir: Path) -> dict[str, Path]:
    """{squashed original filename: renamed path} from ms2_naming_metadata.csv.

    The download step renamed every MS2 file and recorded the mapping. That
    manifest is authoritative -- reconstructing the rename by pattern-matching
    guesses at transformations it actually made (`-` to `_`, parentheses
    dropped), so read it instead wherever it exists.
    """
    manifest = download_dir / "ms2_naming_metadata.csv"
    if not manifest.exists():
        return {}

    # rename_downloads.py may have moved everything since; follow that hop so a
    # stale path here degrades to "not found" rather than a wrong answer
    moves: dict[Path, Path] = {}
    second = download_dir / "file_rename_manifest.csv"
    if second.exists():
        with second.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                moves[Path(row["old_path"])] = Path(row["new_path"])

    mapping: dict[str, Path] = {}
    with manifest.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            original = row.get("original_ms2_filename")
            target = row.get("new_ms2_local_path")
            if not original or not target:
                continue
            path = moves.get(Path(target), Path(target))
            if path.exists() and path.stat().st_size > 0:
                mapping[_squash(original.rsplit(".", 1)[0])] = path
    return mapping


def find_ms2_file(files: dict[str, Path], url: str | None,
                  renames: dict[str, Path] | None = None) -> Path | None:
    """Match an MS2 URL to a downloaded file, allowing for the rename pass.

    Matching is anchored: the normalised original stem must be the whole
    normalised filename or sit at its end behind an underscore boundary. A bare
    containment test would let `3_FeMeOHGa` match `53_FeMeOHGa`. Ambiguous hits
    are rejected rather than guessed at.
    """
    if not url:
        return None
    basename = url.rstrip("/").split("/")[-1]
    if basename in files:
        return files[basename]

    # authoritative: the rename manifest recorded by the download step
    hit = (renames or {}).get(_squash(basename.rsplit(".", 1)[0]))
    if hit is not None:
        return hit

    stem = _normalise_filename(basename.rsplit(".", 1)[0])
    if not stem:
        return None

    matches = [
        path for name, path in files.items()
        if (lambda n: n == stem or n.endswith("_" + stem))(
            _normalise_filename(name.rsplit(".", 1)[0]))
    ]
    if len(matches) == 1:
        return matches[0]
    return None


# --------------------------------------------------------------------------
# GNPS bridge for molecular-family links
# --------------------------------------------------------------------------
# Molecular-family links carry no MS2_URL, so the metabolomics_file -> genome_label
# bridge cannot fire for them. They do carry a GNPS task and componentindex, and
# the task's own result files close the gap:
#
#   componentindex -> cluster indices   (clusterinfosummary...withcomponentID)
#   cluster index  -> internal filename (clusterinfo)
#   internal name  -> original mzML/mzXML (params, upload_file_mapping)
#   original name  -> genome_label       (genome_metabolome_links)

GNPS_BLOCKS = {
    "params": "params/",
    "clusterinfo": "clusterinfo/",
    "summary": "clusterinfosummarygroup_attributes_withIDs_withcomponentID/",
}
GNPS_URL = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task={task}&block=main&file={block}"


def parse_network_nodes_url(url: str | None) -> tuple[str | None, str | None]:
    """Pull (task, componentindex) out of a network_nodes_URL."""
    if not url:
        return None, None
    task = re.search(r"task=([0-9a-f]{32})", url)
    component = re.search(r"componentindex=(-?\d+)", url)
    return (task.group(1) if task else None,
            component.group(1) if component else None)


def fetch_gnps_task(task: str, cache_dir: Path, timeout: int = 180) -> dict[str, Path] | None:
    """Download (and cache) the three result files a task needs. None if it fails."""
    import urllib.error
    import urllib.request

    task_dir = cache_dir / task
    task_dir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for name, block in GNPS_BLOCKS.items():
        path = task_dir / f"{name}.tsv"
        if not path.exists() or path.stat().st_size == 0:
            url = GNPS_URL.format(task=task, block=block)
            try:
                with urllib.request.urlopen(url, timeout=timeout) as response:
                    path.write_bytes(response.read())
            except (urllib.error.URLError, OSError, TimeoutError) as exc:
                print(f"    task {task[:8]} {name}: FAILED ({type(exc).__name__})")
                return None
        out[name] = path
    return out


def gnps_component_files(task: str, cache_dir: Path) -> dict[str, set[str]] | None:
    """Return {componentindex: {original file basenames}} for one GNPS task."""
    files = fetch_gnps_task(task, cache_dir)
    if not files:
        return None

    # internal name -> original basename
    mapping: dict[str, str] = {}
    for entry in re.findall(r"upload_file_mapping[^>]*>([^<]+)<",
                            files["params"].read_text(errors="replace")):
        if "|" in entry:
            internal, original = entry.split("|", 1)
            mapping[Path(internal).name] = Path(original).name

    # cluster index -> component
    cluster_to_component: dict[str, str] = {}
    with files["summary"].open(errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            idx, comp = row.get("cluster index"), row.get("componentindex")
            if idx and comp:
                cluster_to_component[idx] = comp

    # cluster index -> internal filenames, folded up into components
    component_files: dict[str, set[str]] = {}
    with files["clusterinfo"].open(errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            comp = cluster_to_component.get(row.get("#ClusterIdx") or "")
            internal = row.get("#Filename")
            if not comp or not internal:
                continue
            original = mapping.get(Path(internal).name)
            if original:
                component_files.setdefault(comp, set()).add(original)
    return component_files


class GnpsResolver:
    """Lazily fetches and caches GNPS component -> source file maps."""

    def __init__(self, cache_dir: Path):
        self.cache_dir = cache_dir
        self._tasks: dict[str, dict[str, set[str]] | None] = {}

    def files_for(self, url: str | None) -> set[str]:
        task, component = parse_network_nodes_url(url)
        if not task or component is None:
            return set()
        if task not in self._tasks:
            print(f"  resolving GNPS task {task[:8]} ...")
            self._tasks[task] = gnps_component_files(task, self.cache_dir)
        component_map = self._tasks[task]
        if not component_map:
            return set()
        return component_map.get(component, set())


# --------------------------------------------------------------------------
# row building
# --------------------------------------------------------------------------

def project_context(filename: str, payload: dict,
                    studies: dict[str, dict] | None = None) -> dict:
    """Project-level identifiers shared by every row of a project."""
    metabolomics = payload.get("metabolomics", {}) or {}
    project = metabolomics.get("project", {}) or {}
    personal = payload.get("personal", {}) or {}
    genomes = payload.get("genomes", []) or []
    podp_id = filename[:-5] if filename.endswith(".json") else filename
    study = (studies or {}).get(podp_id, {})
    return {
        # descriptor stem keeps the .N version suffix -- NPLinker's podp_id needs it
        "podp_id": podp_id,
        # folder on disk; equals podp_id unless this record was merged into another
        "study_id": study.get("study_id", podp_id),
        "study_members": study.get("study_members", podp_id),
        "study_role": study.get("study_role", "main"),
        "project_uuid": filename.split(".")[0],
        "gnps_massive_id": project.get("GNPSMassIVE_ID"),
        "gnps_task_id": project.get("molecular_network"),
        "gnps_massive_url": project.get("MaSSIVE_URL"),
        "pi_name": personal.get("PI_name"),
        "submitter_email": personal.get("submitter_email"),
        "n_genomes_declared": len(genomes),
        "metcalf_rankable": len(genomes) >= METCALF_MIN_GENOMES,
        "genomes": genomes,
    }


def genome_fields(genomes: list, label: str | None) -> dict:
    """Resolve a genome_label to its accessions. Reuses the locator's resolver."""
    out = {
        "genome_label": label or "",
        "genome_accession": "",
        "genome_accession_type": "",
        "biosample_accession": "",
        "genome_publications": "",
    }
    if not label:
        return out
    resolved = _genome_accession_for_label(genomes, label)
    if resolved:
        out["genome_accession_type"], out["genome_accession"] = resolved
    for genome in genomes:
        if isinstance(genome, dict) and genome.get("genome_label") == label:
            out["biosample_accession"] = genome.get("BioSample_accession") or ""
            out["genome_publications"] = str(genome.get("publications") or "")
            break
    return out


def build_link_rows(
    payloads: dict[str, dict],
    mibig: dict[str, dict],
    downloads: dict[str, dict[str, Path]],
    gnps: "GnpsResolver | None" = None,
    redirects: dict[str, str] | None = None,
    renames: dict[str, Path] | None = None,
    studies: dict[str, dict] | None = None,
) -> list[dict]:
    """One row per BGC_MS2_link across every descriptor."""
    rows: list[dict] = []

    for filename, payload in payloads.items():
        ctx = project_context(filename, payload, studies)
        genomes = ctx["genomes"]
        files = (downloads.get(ctx["study_id"])
                 or downloads.get(ctx["podp_id"])
                 or downloads.get(ctx["project_uuid"], {}))

        # the link -> genome bridge: MS2_URL -> metabolomics_file -> genome_label
        file_to_labels: dict[str, set[str]] = {}
        for gm in payload.get("genome_metabolome_links", []) or []:
            if not isinstance(gm, dict):
                continue
            mfile, label = gm.get("metabolomics_file"), gm.get("genome_label")
            if mfile and label:
                file_to_labels.setdefault(mfile, set()).add(label)

        # same index keyed by basename -- GNPS reports source files by name only
        base_to_labels: dict[str, set[str]] = {}
        for mfile, labels_ in file_to_labels.items():
            base_to_labels.setdefault(Path(mfile).name, set()).update(labels_)

        for link in payload.get("BGC_MS2_links", []) or []:
            bgc_id = link.get("BGC_ID") or {}
            exact = bgc_id.get("MIBiG_number")
            similar = bgc_id.get("similar_MIBiG_number")
            cited_accession = (exact or similar or "").strip()

            # Follow MIBiG's own duplicate redirect rather than editing the
            # source claim: mibig_id stays as PoDP cited it, mibig_id_resolved
            # is what was actually read.
            accession = (redirects or {}).get(cited_accession, cited_accession)
            entry = mibig.get(accession, {})
            redirect_note = ""
            if accession != cited_accession:
                redirect_note = (f"{cited_accession} retired as duplicate of "
                                 f"{accession}; resolved to the active entry")

            smiles_result = compare_structures(link.get("SMILES"), entry.get("compounds", []))
            tier, reason = assign_tier(link, smiles_result, entry)

            ms2_url = link.get("MS2_URL")
            labels = sorted(file_to_labels.get(ms2_url, set())) if ms2_url else []
            genome_source = "ms2_url_bridge" if labels else ""
            gnps_files: list[str] = []

            # molecular-family links have no MS2_URL; go via the GNPS component
            unresolved_reason = ""
            if not labels and gnps is not None:
                gnps_files = sorted(gnps.files_for(link.get("network_nodes_URL")))
                found: set[str] = set()
                for name in gnps_files:
                    found.update(base_to_labels.get(name, set()))
                if found:
                    labels = sorted(found)
                    genome_source = "gnps_component"
                elif gnps_files:
                    # GNPS answered, but none of the component's source files is
                    # one this project declares. Deliberately NOT fuzzy-matched:
                    # these names differ by naming convention or genuinely belong
                    # to another strain, and a near-miss join here would produce
                    # plausible-looking garbage. Candidates are kept in
                    # gnps_component_files for manual adjudication.
                    unresolved_reason = "gnps_files_match_no_declared_metabolomics_file"
                else:
                    unresolved_reason = "gnps_component_empty_or_fetch_failed"
            elif not labels:
                unresolved_reason = ("no_network_nodes_url"
                                     if not link.get("network_nodes_URL")
                                     else "gnps_resolution_not_run")

            gfields = genome_fields(genomes, labels[0] if labels else None)
            if len(labels) > 1:
                gfields["genome_label"] = "|".join(labels)

            genome_path = None
            if gfields["genome_accession"]:
                genome_path = files.get(f"{gfields['genome_accession']}.gbk")
            status_genome = classify_genome_file(genome_path)

            ms2_path = find_ms2_file(files, ms2_url, renames)
            mibig_path = entry.get("path")

            # the four file types. blank name == missing.
            present = {
                "bgc": mibig_path.name if mibig_path else "",
                "genome": genome_path.name if status_genome == "sequence_ok" else "",
                "ms2_run": ms2_path.name if ms2_path else "",
                # nothing extracts single spectra yet; MS2_scan is the pointer
                "ms2_spectrum": "",
            }
            paths = {
                "bgc": str(mibig_path) if mibig_path else "",
                "genome": str(genome_path) if status_genome == "sequence_ok" else "",
                "ms2_run": str(ms2_path) if ms2_path else "",
                "ms2_spectrum": "",
            }
            missing = [k for k in FILE_KINDS if not present[k]]

            first_compound = (entry.get("compounds") or [{}])[0]
            rows.append({
                "podp_id": ctx["podp_id"],
                "study_id": ctx["study_id"],
                "study_role": ctx["study_role"],
                "study_members": ctx["study_members"],
                "tier": tier,
                "tier_reason": reason,
                "data_status": summarise_status(missing),
                "missing_files": "|".join(missing),
                "file_bgc": present["bgc"],
                "file_genome": present["genome"],
                "file_ms2_run": present["ms2_run"],
                "file_ms2_spectrum": present["ms2_spectrum"],
                "path_bgc": paths["bgc"],
                "path_genome": paths["genome"],
                "path_ms2_run": paths["ms2_run"],
                "path_ms2_spectrum": paths["ms2_spectrum"],
                "status_genome": status_genome,
                "mibig_id": cited_accession,
                "mibig_id_resolved": accession,
                "mibig_redirect": redirect_note,
                "mibig_ref_type": "exact" if exact else ("similar" if similar else ""),
                "mibig_local": bool(entry),
                "mibig_version": entry.get("schema", ""),
                "mibig_status": entry.get("status") or "",
                "mibig_quality": entry.get("quality") or "",
                "mibig_retirement_reason": "; ".join(entry.get("retirement_reasons") or []),
                "gnps_massive_id": ctx["gnps_massive_id"] or "",
                "gnps_task_id": ctx["gnps_task_id"] or "",
                "gnps_massive_url": ctx["gnps_massive_url"] or "",
                "compound_name": smiles_result["compound"] or first_compound.get("name") or "",
                "compound_description": link.get("known_link") or "",
                "iupac": link.get("IUPAC") or "",
                "smiles_podp": link.get("SMILES") or "",
                "smiles_mibig": smiles_result["smiles"] or "",
                "smiles_match": smiles_result["match"],
                "inchikey_match": smiles_result["inchikey_match"],
                "inchikey_podp": smiles_result["inchikey_podp"],
                "inchikey_mibig": smiles_result["inchikey_mibig"],
                "smiles_similarity": smiles_result["similarity"],
                "molecular_formula": first_compound.get("formula") or "",
                "bgc_locus_accession": entry.get("locus_accession") or "",
                "bgc_start": entry.get("locus_start") or "",
                "bgc_end": entry.get("locus_end") or "",
                "bgc_organism": entry.get("organism") or "",
                **gfields,
                "genome_ambiguous": len(labels) > 1,
                "genome_source": genome_source,
                "genome_unresolved_reason": unresolved_reason,
                "gnps_component_files": "|".join(gnps_files),
                "ms2_run_url": ms2_url or "",
                "ms2_scan": link.get("MS2_scan") or "",
                "network_nodes_url": link.get("network_nodes_URL") or "",
                "link_type": link.get("link") or "",
                "verification_raw": "; ".join(_as_verification_list(link.get("verification"))),
                "n_genomes_declared": ctx["n_genomes_declared"],
                "metcalf_rankable": ctx["metcalf_rankable"],
                "pi_name": ctx["pi_name"] or "",
                "submitter_email": ctx["submitter_email"] or "",
                "notes": "",
            })
    return rows


def build_pairing_rows(
    payloads: dict[str, dict],
    downloads: dict[str, dict[str, Path]],
    renames: dict[str, Path] | None = None,
    studies: dict[str, dict] | None = None,
) -> list[dict]:
    """One row per genome_metabolome_link. Strain <-> file granularity.

    These are NOT ground truth and carry no tier. PoDP asserts only that the
    genome and the metabolomics file came from the same isolate -- there is no
    BGC <-> spectrum claim to grade. They are the pool a run draws its strain
    axis from, not rows to score against.
    """
    rows: list[dict] = []
    for filename, payload in payloads.items():
        ctx = project_context(filename, payload, studies)
        genomes = ctx["genomes"]
        files = (downloads.get(ctx["study_id"])
                 or downloads.get(ctx["podp_id"])
                 or downloads.get(ctx["project_uuid"], {}))

        for gm in payload.get("genome_metabolome_links", []) or []:
            if not isinstance(gm, dict):
                continue
            gfields = genome_fields(genomes, gm.get("genome_label"))
            ms2_url = gm.get("metabolomics_file")

            genome_path = None
            if gfields["genome_accession"]:
                genome_path = files.get(f"{gfields['genome_accession']}.gbk")
            status_genome = classify_genome_file(genome_path)
            ms2_path = find_ms2_file(files, ms2_url, renames)

            present = {
                "genome": genome_path.name if status_genome == "sequence_ok" else "",
                "ms2_run": ms2_path.name if ms2_path else "",
            }
            missing = [k for k, v in present.items() if not v]

            rows.append({
                "podp_id": ctx["podp_id"],
                "study_id": ctx["study_id"],
                "study_role": ctx["study_role"],
                "study_members": ctx["study_members"],
                "record_type": "strain_pairing",
                "data_status": "complete" if not missing else "incomplete",
                "missing_files": "|".join(missing),
                "file_genome": present["genome"],
                "path_genome": str(genome_path) if present["genome"] else "",
                "status_genome": status_genome,
                "file_ms2_run": present["ms2_run"],
                "path_ms2_run": str(ms2_path) if ms2_path else "",
                "ms2_run_url": ms2_url or "",
                **gfields,
                "gnps_massive_id": ctx["gnps_massive_id"] or "",
                "gnps_task_id": ctx["gnps_task_id"] or "",
                "sample_preparation_label": gm.get("sample_preparation_label") or "",
                "extraction_method_label": gm.get("extraction_method_label") or "",
                "instrumentation_method_label": gm.get("instrumentation_method_label") or "",
                "n_genomes_declared": ctx["n_genomes_declared"],
                "metcalf_rankable": ctx["metcalf_rankable"],
            })
    return rows


# --------------------------------------------------------------------------
# MIBiG 3.1 -> 4.0 migration report
# --------------------------------------------------------------------------

def mibig_migration_report(
    old_index: dict[str, dict],
    new_index: dict[str, dict],
    link_rows: list[dict],
) -> list[dict]:
    """One row per cited MIBiG accession, recording exactly what the move changed.

    Compares the two releases only over accessions the ground truth actually
    cites -- a full 2502-vs-3013 diff is noise for this purpose. Every row says
    what happened to the entry and, where it matters, what it does to the table.
    """
    cited: dict[str, int] = Counter(r["mibig_id"] for r in link_rows if r["mibig_id"])
    report: list[dict] = []

    for accession in sorted(cited):
        old = old_index.get(accession)
        new = new_index.get(accession)

        def has_structure(entry: dict | None) -> bool:
            return bool(entry and any(c.get("smiles") for c in entry["compounds"]))

        if old and not new:
            fate = "dropped_from_4.0"
        elif new and not old:
            fate = "added_in_4.0"
        elif not old and not new:
            fate = "absent_from_both"
        else:
            fate = "carried_over"

        structure_change = ""
        if has_structure(new) and not has_structure(old):
            structure_change = "gained_structure"
        elif has_structure(old) and not has_structure(new):
            structure_change = "lost_structure"

        report.append({
            "mibig_id": accession,
            "ground_truth_rows": cited[accession],
            "fate": fate,
            "status_3.1": (old or {}).get("status") or ("" if old else "ABSENT"),
            "status_4.0": (new or {}).get("status") or ("" if new else "ABSENT"),
            "quality_4.0": (new or {}).get("quality") or "",
            "retirement_reason_4.0": "; ".join((new or {}).get("retirement_reasons") or []),
            "structure_3.1": "yes" if has_structure(old) else "no",
            "structure_4.0": "yes" if has_structure(new) else "no",
            "structure_change": structure_change,
            "compound_3.1": "; ".join(
                filter(None, (c.get("name") for c in (old or {}).get("compounds", [])))),
            "compound_4.0": "; ".join(
                filter(None, (c.get("name") for c in (new or {}).get("compounds", [])))),
            "locus_3.1": (old or {}).get("locus_accession") or "",
            "locus_4.0": (new or {}).get("locus_accession") or "",
        })
    return report


def print_migration_summary(report: list[dict]) -> None:
    rows = sum(r["ground_truth_rows"] for r in report)
    print(f"\n--- MIBiG 3.1 -> 4.0, over the {len(report)} cited accessions "
          f"({rows} ground_truth rows) " + "-" * 6)

    for label, key in (("fate", "fate"), ("status in 4.0", "status_4.0")):
        print(f"  by {label}:")
        for value, count in Counter(r[key] for r in report).most_common():
            affected = sum(r["ground_truth_rows"] for r in report if r[key] == value)
            print(f"    {value:<22} {count:>3} accessions  ({affected} rows)")

    changed = [r for r in report if r["structure_change"]]
    print(f"  structure changes: {len(changed)} accessions")
    for r in changed:
        print(f"    {r['mibig_id']}  {r['structure_change']:<17} "
              f"{r['ground_truth_rows']} row(s)  {r['compound_4.0'] or r['compound_3.1']}")

    flagged = [r for r in report if r["status_4.0"] not in ("active", "ABSENT", "")]
    if flagged:
        print(f"  cited but NOT active in 4.0: {len(flagged)} accessions")
        for r in flagged:
            print(f"    {r['mibig_id']}  {r['status_4.0']:<9} {r['ground_truth_rows']} row(s)  "
                  f"{r['retirement_reason_4.0'] or '(no reason given)'}")

    missing = [r for r in report if r["status_4.0"] == "ABSENT"]
    if missing:
        print(f"  cited but absent from the 4.0 bundle: {len(missing)}")
        for r in missing:
            print(f"    {r['mibig_id']}  {r['ground_truth_rows']} row(s)  "
                  f"was: {r['compound_3.1'] or '(unnamed)'}")


# --------------------------------------------------------------------------
# funnel + output
# --------------------------------------------------------------------------

def build_funnel(
    payloads: dict[str, dict],
    skipped: list[str],
    link_rows: list[dict],
    pairing_rows: list[dict],
) -> list[dict]:
    tiers = Counter(r["tier"] for r in link_rows)
    with_links = sum(1 for p in payloads.values() if p.get("BGC_MS2_links"))
    rankable = sum(1 for r in link_rows if r["metcalf_rankable"])
    stages = [
        ("descriptors found", len(payloads) + len(skipped), "files in json_descriptors/"),
        ("descriptors parsed", len(payloads), f"skipped: {', '.join(skipped) or 'none'}"),
        ("projects with BGC_MS2_links", with_links, "optional PoDP curation field"),
        ("BGC-MS2 links total", len(link_rows), "ground_truth.csv rows"),
        ("tier 1 (experimental)", tiers[1], "verification startswith 'Experimentally validated'"),
        ("tier 2 (structural match)", tiers[2],
         "canonical InChIKey match (full key or connectivity layer)"),
        ("tier 3 (no evidence)", tiers[3], "retained, not dropped"),
        ("links in rankable projects", rankable, f">= {METCALF_MIN_GENOMES} declared genomes"),
        ("links data_status=complete", sum(1 for r in link_rows if r["data_status"] == "complete"),
         "all four file types on disk"),
        ("links spectrum_pending", sum(1 for r in link_rows if r["data_status"] == "spectrum_pending"),
         "only the un-extracted MS2 spectrum missing"),
        ("genome via MS2_URL bridge",
         sum(1 for r in link_rows if r["genome_source"] == "ms2_url_bridge"),
         "metabolomics_file -> genome_label"),
        ("genome via GNPS component",
         sum(1 for r in link_rows if r["genome_source"] == "gnps_component"),
         "molecular-family links, needs --resolve-gnps"),
        ("genome unresolved",
         sum(1 for r in link_rows if not r["genome_source"]), "no bridge fired"),
        ("USABLE: tier 1 + rankable + genome on disk",
         sum(1 for r in link_rows
             if r["tier"] == 1 and r["metcalf_rankable"]
             and r["status_genome"] == "sequence_ok"),
         "experimentally validated, scoreable, sequence present"),
        ("strain pairings (untiered)", len(pairing_rows),
         "strain_pairings.csv rows -- declared co-origin, not ground truth"),
    ]
    return [{"stage": s, "count": c, "filter": f} for s, c, f in stages]


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def report_screen_disagreement(rows: list[dict]) -> None:
    """Compare the (diagnostic) string screen against the authoritative InChIKey.

    Tiers no longer depend on the screen, so this is pure instrumentation -- it
    shows how far a chemistry-toolkit-free heuristic would have drifted.
    """
    comparable = [r for r in rows if r["inchikey_match"] != "uncomputable"]
    strong = lambda r: r["smiles_match"] in ("identical", "permutation")
    matched = lambda r: r["inchikey_match"] in ("full", "skeleton")

    agree = sum(1 for r in comparable if strong(r) == (r["inchikey_match"] == "full"))
    print(f"\n--- string screen vs InChIKey (diagnostic only) " + "-" * 23)
    print(f"  comparable rows              {len(comparable)}")
    print(f"  agree / disagree             {agree} / {len(comparable) - agree}")

    weak = [r for r in comparable if r["smiles_match"] in ("formula", "partial")]
    refuted = [r for r in weak if not matched(r)]
    print(f"  weak screen calls            {len(weak)}  "
          f"({len(refuted)} refuted by rdkit)")
    missed = [r for r in comparable if r["smiles_match"] == "none" and matched(r)]
    print(f"  screen false negatives       {len(missed)}")
    for verdict in ("identical", "permutation"):
        subset = [r for r in comparable if r["smiles_match"] == verdict]
        ok = sum(1 for r in subset if r["inchikey_match"] == "full")
        if subset:
            print(f"  screen '{verdict}'{'':<{13 - len(verdict)}}{ok}/{len(subset)} confirmed full")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="print the funnel without writing any CSV")
    parser.add_argument("--migration-report", action="store_true",
                        help="diff the MIBiG 3.1 and 4.0 releases over the cited "
                             "accessions and write mibig_migration_report.csv")
    parser.add_argument("--mibig-31-dir", type=Path, default=MIBIG_31_DIR,
                        help="MIBiG 3.1 tree, used only by --migration-report")
    parser.add_argument("--resolve-gnps", action="store_true",
                        help="resolve molecular-family links to genomes via the "
                             "GNPS task (network; cached under data/PoDP/gnps_cache)")
    parser.add_argument("--podp-dir", type=Path, default=PODP_DIR)
    parser.add_argument("--mibig-dir", type=Path, default=MIBIG_DIR)
    parser.add_argument("--download-dir", type=Path, default=DOWNLOAD_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    if not args.podp_dir.is_dir():
        print(f"error: {args.podp_dir} not found", file=sys.stderr)
        return 1

    try:
        import rdkit  # noqa: F401
    except ImportError:
        print("error: rdkit is required -- tier 2 is decided by canonical InChIKey.\n"
              "  conda install -n nplinker-application -c conda-forge rdkit\n"
              "  (or: pip install rdkit)", file=sys.stderr)
        return 1

    payloads, skipped = load_descriptors(args.podp_dir)
    mibig_dir = args.mibig_dir
    if not (mibig_dir.is_dir() and any(mibig_dir.glob("BGC*.json"))):
        mibig_dir = ensure_mibig(MIBIG_VERSION, MIBIG_CACHE)
    mibig = build_mibig_index(mibig_dir)
    downloads = index_downloads(args.download_dir)

    print(f"descriptors parsed : {len(payloads)}")
    schemas = Counter(e.get("schema") for e in mibig.values())
    print(f"MIBiG entries      : {len(mibig)}  schema {dict(schemas)}  ({mibig_dir})")
    statuses = Counter(e.get("status") for e in mibig.values())
    if statuses:
        print(f"MIBiG status       : {dict(statuses)}")
    print(f"download dirs      : {len(downloads)}")
    if skipped:
        print(f"skipped            : {', '.join(skipped)}")

    resolver = None
    if args.resolve_gnps:
        print("\n--- resolving molecular-family links via GNPS " + "-" * 24)
        resolver = GnpsResolver(args.out_dir / "gnps_cache")

    redirects = build_redirect_map(mibig)
    if redirects:
        print(f"MIBiG redirects    : {len(redirects)} retired-as-duplicate entries "
              f"with an active target")

    studies = resolve_studies(payloads)
    merged = {k: v for k, v in studies.items() if v["study_role"] == "merged"}
    if merged:
        print(f"study grouping     : {len(merged)} project(s) merged into a shared "
              f"study folder (identical MassIVE set)")
        for podp_id, info in sorted(merged.items()):
            print(f"    {podp_id}  ->  {info['study_id']}")

    renames = load_rename_manifest(args.download_dir)
    if renames:
        print(f"MS2 rename manifest: {len(renames)} original -> renamed files")

    link_rows = build_link_rows(payloads, mibig, downloads,
                                gnps=resolver, redirects=redirects, renames=renames,
                                studies=studies)
    pairing_rows = build_pairing_rows(payloads, downloads, renames=renames,
                                      studies=studies)
    funnel = build_funnel(payloads, skipped, link_rows, pairing_rows)

    print("\n--- funnel " + "-" * 58)
    for stage in funnel:
        print(f"  {stage['stage']:<32} {stage['count']:>6}   {stage['filter']}")

    print("\n--- genome file status, ground_truth rows " + "-" * 28)
    for status, count in Counter(r["status_genome"] for r in link_rows).most_common():
        print(f"  {status:<28} {count:>5}")

    print("\n--- SMILES match, all rows " + "-" * 43)
    for match, count in Counter(r["smiles_match"] for r in link_rows).most_common():
        print(f"  {match:<28} {count:>5}")

    tier1_nomatch = sum(1 for r in link_rows if r["tier"] == 1 and r["smiles_match"] == "none")
    print(f"\nQC: tier 1 rows whose SMILES does not match their cited MIBiG entry: {tier1_nomatch}")

    report: list[dict] = []
    if args.migration_report:
        old_index = build_mibig_index(args.mibig_31_dir)
        if not old_index:
            print(f"\nmigration report: {args.mibig_31_dir} not found -- skipped")
        else:
            report = mibig_migration_report(old_index, mibig, link_rows)
            print_migration_summary(report)

    report_screen_disagreement(link_rows)

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return 0

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "ground_truth.csv", link_rows)
    write_csv(args.out_dir / "strain_pairings.csv", pairing_rows)
    if report:
        migration_path = args.out_dir / "reports" / "mibig_migration_report.csv"
        migration_path.parent.mkdir(parents=True, exist_ok=True)
        write_csv(migration_path, report)
        print(f"wrote mibig_migration_report.csv ({len(report)} accessions)")
    print(f"\nwrote ground_truth.csv    ({len(link_rows)} rows)")
    print(f"wrote strain_pairings.csv ({len(pairing_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
