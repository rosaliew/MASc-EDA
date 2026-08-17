import json
from pathlib import Path

from scripts.podp_online_file_locator import (
    build_download_plan,
    collect_verified_bgc_ms2_genome_links,
    collect_verified_bgc_ms2_links,
    extract_project_links,
)


def test_extract_project_links_from_podp_json():
    payload = {
        "metabolomics": {
            "project": {
                "GNPSMassIVE_ID": "MSV000084950",
                "MaSSIVE_URL": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=07d6a1a51961435bb3c1126674857790",
            }
        },
        "genomes": [
            {
                "genome_label": "Test strain",
                "genome_ID": {"GenBank_accession": "GCF_000203835.1"},
            }
        ],
    }

    result = extract_project_links(payload)

    assert result["massive_id"] == "MSV000084950"
    assert result["massive_url"] == "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=07d6a1a51961435bb3c1126674857790"
    assert result["genomes"][0]["genome_accession"] == "GCF_000203835.1"
    assert result["genomes"][0]["ncbi_assembly_url"] == "https://www.ncbi.nlm.nih.gov/assembly/?term=GCF_000203835.1"
    assert result["genomes"][0]["ncbi_nuccore_url"] == "https://www.ncbi.nlm.nih.gov/nuccore/GCF_000203835.1"


def test_build_download_plan_includes_fasta_gbk_and_ms2():
    payload = {
        "metabolomics": {
            "project": {"GNPSMassIVE_ID": "MSV000084950"},
            "files": [{"metabolomics_file": "ftp://massive.ucsd.edu/MSV000084950/raw/sample.mzML"}],
        },
        "genomes": [{
            "genome_label": "Test strain",
            "genome_ID": {"GenBank_accession": "GCF_000203835.1"},
        }],
    }

    plan = build_download_plan(payload, project_id="demo_project", output_dir="/tmp/demo_output")

    assert plan["project_id"] == "demo_project"
    assert any(request["kind"] == "fasta" for request in plan["genome_requests"])
    assert any(request["kind"] == "gbk" for request in plan["genome_requests"])
    assert any(request["kind"] == "ms2" for request in plan["ms2_requests"])
    assert plan["ms2_requests"][0]["url"].startswith("ftp://massive.ucsd.edu/")


def test_collect_verified_bgc_ms2_links_filters_verified_massive_entries():
    payload = {
        "_project_id": "demo_project",
        "BGC_MS2_links": [
            {
                "link": "single molecule",
                "MS2_URL": "ftp://massive.ucsd.edu/MSV000084950/peak/040220001.mzML",
                "MS2_scan": "2758",
                "verification": ["Experimentally validated with NMR and/or detailed MS/MS analysis"],
            },
            {
                "link": "single molecule",
                "MS2_URL": "ftp://massive.ucsd.edu/MSV000084950/peak/040220002.mzML",
                "verification": [],
            },
            {
                "link": "GNPS molecular family",
                "network_nodes_URL": "https://gnps.ucsd.edu/ProteoSAFe/result.jsp",
                "verification": ["Experimentally validated with NMR and/or detailed MS/MS analysis"],
            },
        ],
    }

    result = collect_verified_bgc_ms2_links(payload)

    assert len(result) == 1
    assert result[0]["MS2_URL"] == "ftp://massive.ucsd.edu/MSV000084950/peak/040220001.mzML"
    assert result[0]["project_id"] == "demo_project"


def test_collect_verified_bgc_ms2_genome_links_maps_ms2_url_to_genome_accession():
    payload = {
        "_project_id": "demo_project",
        "genomes": [
            {
                "genome_label": "Test strain A",
                "genome_ID": {"GenBank_accession": "GCF_000203835.1"},
            },
            {
                "genome_label": "Test strain B",
                "genome_ID": {"GenBank_accession": "GCF_000000002.1"},
            },
        ],
        "genome_metabolome_links": [
            {
                "genome_label": "Test strain A",
                "metabolomics_file": "ftp://massive.ucsd.edu/MSV000084950/peak/040220001.mzML",
            },
            {
                "genome_label": "Test strain B",
                "metabolomics_file": "ftp://massive.ucsd.edu/MSV000084950/peak/unlinked.mzML",
            },
        ],
        "BGC_MS2_links": [
            {
                "link": "single molecule",
                "MS2_URL": "ftp://massive.ucsd.edu/MSV000084950/peak/040220001.mzML",
                "verification": ["Experimentally validated with NMR and/or detailed MS/MS analysis"],
            },
        ],
    }

    result = collect_verified_bgc_ms2_genome_links(payload)

    assert len(result) == 1
    assert result[0]["genome_label"] == "Test strain A"
    assert result[0]["genome_accession"] == "GCF_000203835.1"
    assert result[0]["genome_accession_type"] == "GenBank_accession"
