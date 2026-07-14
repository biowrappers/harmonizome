"""Focused tests for the Harmonizome click interface."""

import json
from pathlib import Path

from click.testing import CliRunner
import pytest

import harmonizome.harmonizome_cli as cli_module
from harmonizome.harmonizome_cli import main


def test_cli_without_subcommand_shows_help():
    runner = CliRunner()

    result = runner.invoke(main)

    assert result.exit_code == 1
    assert "Welcome to Harmonizome CLI" in result.output
    assert "Commands:" in result.output


def test_cli_list_datasets_runs(monkeypatch):
    runner = CliRunner()
    monkeypatch.setattr(
        "harmonizome.harmonizome_cli.Harmonizome.DATASETS",
        ["ENCODE"],
    )
    monkeypatch.setattr(
        "harmonizome.harmonizome._get_dataset_to_path",
        lambda: {"ENCODE": "encode/path"},
        raising=False,
    )

    result = runner.invoke(main, ["list-datasets"])

    assert result.exit_code == 0
    assert "Available Harmonizome Datasets" in result.output
    assert "ENCODE" in result.output


def test_cli_get_entity_runs(monkeypatch):
    runner = CliRunner()
    monkeypatch.setattr(
        "harmonizome.harmonizome_cli.Harmonizome.get",
        lambda entity_type, name: {
            "symbol": name,
            "description": "Test description",
        },
        raising=False,
    )

    result = runner.invoke(main, ["get-entity", "gene", "STAT3"])

    assert result.exit_code == 0
    assert "Gene: STAT3" in result.output
    assert "STAT3" in result.output


def test_display_key_name_uses_overrides():
    assert cli_module.display_key_name("uniprotId") == "UniProt ID"
    assert cli_module.display_key_name("gene_symbol") == "Gene Symbol"


def test_format_value_covers_supported_shapes():
    long_text = "A" * 130

    assert cli_module.format_value("description", long_text) == long_text
    assert cli_module.format_value("name", long_text) == f"{'A' * 120}..."
    assert cli_module.format_value("synonyms", ["A", "B"]) == "A, B"
    assert cli_module.format_value(
        "proteins",
        [{"symbol": "STAT3"}, {}],
    ) == "  • STAT3\n  • Unknown"
    assert cli_module.format_value(
        "hgncRootFamilies",
        [{"name": "Kinase"}, {}],
    ) == "  • Kinase\n  • Unknown"
    assert cli_module.format_value(
        "geneSets",
        [{"name": f"Set{i}"} for i in range(6)],
    ) == (
        "  • Set0\n  • Set1\n  • Set2\n  • Set3\n  • Set4\n  ... and 1 more"
    )
    assert cli_module.format_value("other", [1, 2, 3]) == "1, 2, 3"
    assert cli_module.format_value("other", [1, 2, 3, 4]) == "1, 2, 3 ... and 1 more"
    assert cli_module.format_value("other", {"href": "/x", "name": "ENCODE", "id": 5}) == (
        "ENCODE (ID: 5)"
    )
    assert cli_module.format_value("other", {"a": 1, "b": 2}) == "<2 items>"
    assert cli_module.format_value("other", 1.25) == "1.25"


def test_find_dataset_by_partial_name_matches_case_insensitively(monkeypatch):
    monkeypatch.setattr(
        cli_module,
        "_get_dataset_to_path",
        lambda: {
            "ENCODE TF": "encode/tf",
            "GTEx Tissue": "gtex/tissue",
            "Other": "other",
        },
    )

    matches = cli_module.find_dataset_by_partial_name("encode")

    assert matches == ["ENCODE TF"]


def test_cli_get_entity_shows_extra_fields_and_multiline_values(monkeypatch):
    runner = CliRunner()
    monkeypatch.setattr(
        "harmonizome.harmonizome_cli.Harmonizome.get",
        lambda entity_type, name: {
            "symbol": name,
            "description": "Test description",
            "proteins": [{"symbol": "P1"}, {}],
            "customField": "custom value",
        },
        raising=False,
    )

    result = runner.invoke(main, ["get-entity", "gene", "STAT3"])

    assert result.exit_code == 0
    assert "Proteins:" in result.output
    assert "  • P1" in result.output
    assert "  • Unknown" in result.output
    assert "Customfield:" in result.output
    assert "custom value" in result.output


def test_cli_get_entity_handles_errors(monkeypatch):
    runner = CliRunner()

    def raise_error(entity_type: str, name: str) -> dict[str, str]:
        raise RuntimeError("broken lookup")

    monkeypatch.setattr(
        "harmonizome.harmonizome_cli.Harmonizome.get",
        raise_error,
        raising=False,
    )

    result = runner.invoke(main, ["get-entity", "gene", "STAT3"])

    assert result.exit_code == 1
    assert "Error getting gene 'STAT3': broken lookup" in result.output


def test_cli_download_expands_short_codes_and_restores_working_directory(
    monkeypatch, tmp_path: Path
):
    runner = CliRunner()
    original_working_directory = Path.cwd()
    output_dir = tmp_path / "downloads"
    captured: dict[str, list[str] | str] = {}

    def fake_find_dataset_by_partial_name(partial_name: str) -> list[str]:
        if partial_name == "ENCODE":
            return ["ENCODE TF", "ENCODE Histone"]
        return []

    def fake_download(datasets: list[str]) -> list[str]:
        captured["datasets"] = datasets
        captured["cwd"] = str(Path.cwd())
        return ["first.tsv", "second.tsv"]

    monkeypatch.setattr(cli_module, "find_dataset_by_partial_name", fake_find_dataset_by_partial_name)
    monkeypatch.setattr(cli_module.Harmonizome, "download", fake_download)

    result = runner.invoke(
        main,
        ["download", "ENCODE", "Other Dataset", "--output-dir", str(output_dir)],
    )

    assert result.exit_code == 0
    assert captured["datasets"] == ["ENCODE TF", "ENCODE Histone", "Other Dataset"]
    assert captured["cwd"] == str(output_dir)
    assert Path.cwd() == original_working_directory
    assert output_dir.exists()
    assert "Found 2 datasets matching 'ENCODE':" in result.output
    assert "Downloaded: first.tsv" in result.output
    assert "Downloaded: second.tsv" in result.output


def test_cli_download_handles_keyboard_interrupt(monkeypatch):
    runner = CliRunner()

    def interrupt_download(datasets: list[str]) -> list[str]:
        raise KeyboardInterrupt

    monkeypatch.setattr(cli_module.Harmonizome, "download", interrupt_download)

    result = runner.invoke(main, ["download", "ENCODE"])

    assert result.exit_code == 1
    assert "Download interrupted by user." in result.output


def test_cli_download_keeps_short_code_when_no_dataset_matches(monkeypatch):
    runner = CliRunner()
    captured: dict[str, list[str]] = {}

    monkeypatch.setattr(cli_module, "find_dataset_by_partial_name", lambda partial_name: [])

    def capture_download(datasets: list[str]) -> list[str]:
        captured["datasets"] = datasets
        return []

    monkeypatch.setattr(cli_module.Harmonizome, "download", capture_download)

    result = runner.invoke(main, ["download", "GTEX"])

    assert result.exit_code == 0
    assert captured["datasets"] == ["GTEX"]


def test_cli_download_handles_generic_errors(monkeypatch):
    runner = CliRunner()

    def fail_download(datasets: list[str]) -> list[str]:
        raise RuntimeError("network failure")

    monkeypatch.setattr(cli_module.Harmonizome, "download", fail_download)

    result = runner.invoke(main, ["download", "ENCODE"])

    assert result.exit_code == 1
    assert "Error downloading datasets: network failure" in result.output


def test_cli_functional_associations_writes_output_file(
    monkeypatch, tmp_path: Path
):
    runner = CliRunner()
    output_file = tmp_path / "functional-associations.json"

    datasets = []
    for index in range(6):
        datasets.append(
            {
                "dataset": f"Dataset {index}",
                "summary": "S" * 140,
                "associations": [
                    {
                        "description": "Increased",
                        "items": [
                            {"name": f"Gene {item_index}", "score": float(item_index)}
                            for item_index in range(4)
                        ],
                    }
                ],
            }
        )

    expected_results = {
        "gene_info": {
            "symbol": "STAT3",
            "name": "Signal transducer",
            "ncbi_id": "6774",
            "description": "D" * 240,
        },
        "functional_associations": {
            "total_datasets": 6,
            "total_associations": 24,
            "total_increased": 12,
            "total_decreased": 12,
            "datasets": datasets,
        },
    }

    monkeypatch.setattr(
        cli_module.Harmonizome,
        "get_gene_functional_annotations",
        lambda gene_symbol, requested_datasets: expected_results,
    )

    result = runner.invoke(
        main,
        [
            "functional-associations",
            "STAT3",
            "--datasets",
            "ENCODE TF",
            "--output-file",
            str(output_file),
        ],
    )

    assert result.exit_code == 0
    assert "Datasets: ENCODE TF" in result.output
    assert "Gene Information:" in result.output
    assert "Functional Associations Summary:" in result.output
    assert "Dataset Details:" in result.output
    assert "... and 1 more datasets" in result.output
    assert "... and 1 more" in result.output
    assert "Results saved to:" in result.output
    assert json.loads(output_file.read_text(encoding="utf-8")) == expected_results


def test_cli_functional_associations_handles_errors(monkeypatch):
    runner = CliRunner()

    def raise_error(gene_symbol: str, datasets: list[str] | None) -> dict[str, str]:
        raise RuntimeError("annotation lookup failed")

    monkeypatch.setattr(
        cli_module.Harmonizome,
        "get_gene_functional_annotations",
        raise_error,
    )

    result = runner.invoke(main, ["functional-associations", "STAT3"])

    assert result.exit_code == 1
    assert "Error getting functional associations: annotation lookup failed" in result.output
