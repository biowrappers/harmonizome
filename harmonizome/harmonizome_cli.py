#!/usr/bin/env python3
"""Command-line interface for the Harmonizome package."""

import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Optional

import click

from harmonizome import Harmonizome
from harmonizome.harmonizome import _get_dataset_to_path

FIELD_ORDER: dict[str, list[str]] = {
    "gene": [
        "symbol",
        "name",
        "synonyms",
        "description",
        "ncbiEntrezGeneId",
        "ncbiEntrezGeneUrl",
        "proteins",
        "hgncRootFamilies",
    ],
    "protein": [
        "symbol",
        "name",
        "description",
        "uniprotId",
        "uniprotUrl",
        "genes",
    ],
    "dataset": ["name", "description", "version", "resource"],
    "attribute": ["name", "description", "dataset", "resource"],
    "gene_set": ["name", "description", "dataset", "genes"],
    "resource": ["name", "description", "url", "version"],
}

DISPLAY_KEY_OVERRIDES = {
    "ncbiEntrezGeneId": "NCBI Entrez Gene ID",
    "ncbiEntrezGeneUrl": "NCBI Entrez Gene URL",
    "uniprotId": "UniProt ID",
    "uniprotUrl": "UniProt URL",
}


def write_header(text: str) -> None:
    """Render a section header using click styling."""
    click.secho(f"\n{'=' * 60}", fg="magenta", bold=True)
    click.secho(f"  {text}", fg="magenta", bold=True)
    click.secho(f"{'=' * 60}\n", fg="magenta", bold=True)


def write_status(symbol: str, text: str, color: str) -> None:
    """Render a compact status line with consistent coloring."""
    click.secho(f"{symbol} {text}", fg=color)


def write_dataset_item(index: int, dataset: str, short_code: str) -> None:
    """Render one dataset list entry."""
    click.secho(f"{index:3d}.", fg="cyan", nl=False)
    click.secho(f" {dataset}", bold=True)
    click.secho("       Short code:", fg="yellow", nl=False)
    click.echo(f" {short_code}\n")


def display_key_name(key: str) -> str:
    """Convert raw response keys into human-readable labels."""
    return DISPLAY_KEY_OVERRIDES.get(key, key.replace("_", " ").title())


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def list_datasets() -> None:
    """List all available datasets."""
    dataset_to_path = _get_dataset_to_path()
    write_header("Available Harmonizome Datasets")
    write_status("i", f"Total datasets: {len(Harmonizome.DATASETS)}", "blue")
    write_status(
        "!",
        "Use the full dataset name for downloading (e.g., 'ENCODE Transcription Factor Binding Site Profiles')",
        "yellow",
    )
    click.echo()

    for i, dataset in enumerate(Harmonizome.DATASETS, 1):
        short_code = dataset_to_path.get(dataset, "N/A")
        write_dataset_item(i, dataset, short_code)


def format_value(key: str, value: Any) -> str:
    """Format a value for display based on its type and key."""
    if isinstance(value, str):
        # Don't truncate descriptions
        if key == "description":
            return value
        # Truncate other long strings
        elif len(value) > 120:
            return value[:120] + "..."
        return value
    elif isinstance(value, list):
        if key == "synonyms":
            # Show all synonyms
            return ", ".join(value)
        elif key == "proteins":
            return "\n".join([f"  • {p.get('symbol', 'Unknown')}" for p in value])
        elif key == "hgncRootFamilies":
            return "\n".join([f"  • {f.get('name', 'Unknown')}" for f in value])
        elif key == "geneSets":
            return "\n".join([f"  • {g.get('name', 'Unknown')}" for g in value[:5]]) + (
                f"\n  ... and {len(value) - 5} more" if len(value) > 5 else ""
            )
        else:
            if len(value) > 3:
                return (
                    ", ".join(str(v) for v in value[:3])
                    + f" ... and {len(value) - 3} more"
                )
            return ", ".join(str(v) for v in value)
    elif isinstance(value, dict):
        if "href" in value:
            return f"{value.get('name', 'Unknown')} (ID: {value.get('id', 'N/A')})"
        else:
            return f"<{len(value)} items>"
    else:
        return str(value)


def get_entity_info(entity_type: str, name: str) -> None:
    """Get information about a specific entity."""
    try:
        info = Harmonizome.get(entity_type, name)
        write_header(f"{entity_type.title()}: {name}")

        ordered_fields = FIELD_ORDER.get(entity_type, list(info.keys()))

        # Add any missing fields to the end
        for key in info.keys():
            if key not in ordered_fields:
                ordered_fields.append(key)

        for key in ordered_fields:
            if key in info:
                value = info[key]
                formatted_value = format_value(key, value)
                click.secho(f"{display_key_name(key)}:", bold=True)
                if "\n" in formatted_value:
                    click.echo(formatted_value)
                else:
                    click.echo(f"  {formatted_value}")
                click.echo()

    except Exception as e:
        write_status("x", f"Error getting {entity_type} '{name}': {e}", "red")
        sys.exit(1)


def find_dataset_by_partial_name(partial_name: str) -> list[str]:
    """Find datasets that contain the partial name."""
    dataset_to_path = _get_dataset_to_path()
    matching_datasets = []
    for dataset_name in dataset_to_path:
        if partial_name.upper() in dataset_name.upper():
            matching_datasets.append(dataset_name)
    return matching_datasets


def download_datasets(datasets: list[str], output_dir: Optional[str] = None) -> None:
    """Download specified datasets."""
    write_header("Harmonizome Dataset Download")

    if output_dir:
        write_status("i", f"Output directory: {output_dir}", "blue")
        Path(output_dir).mkdir(exist_ok=True)
        original_working_directory = Path.cwd()
        os.chdir(output_dir)
    else:
        write_status("i", "Output directory: current directory", "blue")
        original_working_directory = None

    # Check if any datasets need to be expanded (e.g., "ENCODE" -> all ENCODE datasets)
    expanded_datasets = []
    for dataset in datasets:
        if dataset.upper() in ["ENCODE", "GTEX", "MSIGDB"]:
            # Find all datasets containing this name
            matching = find_dataset_by_partial_name(dataset)
            if matching:
                write_status(
                    "+",
                    f"Found {len(matching)} datasets matching '{dataset}':",
                    "green",
                )
                for match in matching:
                    click.secho("  •", fg="cyan", nl=False)
                    click.echo(f" {match}")
                expanded_datasets.extend(matching)
            else:
                expanded_datasets.append(dataset)
        else:
            expanded_datasets.append(dataset)

    write_status("i", f"Total datasets to download: {len(expanded_datasets)}", "blue")
    click.echo()

    try:
        for filename in Harmonizome.download(expanded_datasets):
            write_status("+", f"Downloaded: {filename}", "green")
    except KeyboardInterrupt:
        write_status("!", "Download interrupted by user.", "yellow")
        sys.exit(1)
    except Exception as e:
        write_status("x", f"Error downloading datasets: {e}", "red")
        sys.exit(1)
    finally:
        if original_working_directory is not None:
            os.chdir(original_working_directory)


def get_functional_associations(
    gene_symbol: str,
    datasets: Optional[list[str]] = None,
    output_file: Optional[str] = None,
) -> None:
    """Get functional associations for a gene."""
    write_header(f"Functional Associations for {gene_symbol}")

    if datasets:
        write_status("i", f"Datasets: {', '.join(datasets)}", "blue")
    else:
        write_status("i", "Datasets: All available", "blue")

    try:
        # Get functional associations
        results = Harmonizome.get_gene_functional_annotations(gene_symbol, datasets)

        # Display results
        gene_info = results["gene_info"]
        func_assoc = results["functional_associations"]

        click.secho("\nGene Information:", bold=True)
        click.echo(f"  Symbol: {gene_info['symbol']}")
        click.echo(f"  Name: {gene_info['name']}")
        click.echo(f"  NCBI ID: {gene_info['ncbi_id']}")
        click.echo(f"  Description: {gene_info['description'][:200]}...")

        click.secho("\nFunctional Associations Summary:", bold=True)
        click.echo(f"  Total datasets: {func_assoc['total_datasets']}")
        click.echo(f"  Total associations: {func_assoc['total_associations']}")
        click.echo(f"  Increased associations: {func_assoc['total_increased']}")
        click.echo(f"  Decreased associations: {func_assoc['total_decreased']}")

        # Show dataset details
        click.secho("\nDataset Details:", bold=True)
        for dataset in func_assoc["datasets"][:5]:  # Show first 5 datasets
            click.secho(f"\nDataset: {dataset['dataset']}", fg="cyan")
            click.echo(f"  Summary: {dataset['summary'][:100]}...")

            for assoc_group in dataset["associations"]:
                click.echo(f"  {assoc_group['description']}:")
                for item in assoc_group["items"][:3]:  # Show first 3
                    click.echo(f"    {item['name']} [{item['score']:.5f}]")
                if len(assoc_group["items"]) > 3:
                    click.echo(f"    ... and {len(assoc_group['items']) - 3} more")

        if len(func_assoc["datasets"]) > 5:
            click.secho(
                f"\n... and {len(func_assoc['datasets']) - 5} more datasets",
                fg="yellow",
            )

        # Save to file if requested
        if output_file:
            with Path(output_file).open("w", encoding="utf-8") as file_handle:
                json.dump(results, file_handle, indent=2)
            write_status("+", f"Results saved to: {output_file}", "green")

    except Exception as e:
        write_status("x", f"Error getting functional associations: {e}", "red")
        sys.exit(1)


ENTITY_TYPE_CHOICES = click.Choice(
    ["gene", "gene_set", "attribute", "dataset", "protein", "resource"],
    case_sensitive=False,
)


@click.group(
    invoke_without_command=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose logging.")
@click.pass_context
def main(ctx: click.Context, verbose: bool) -> None:
    """Access harmonized datasets of genes and proteins."""
    setup_logging(verbose)

    if ctx.invoked_subcommand is None:
        write_header("Welcome to Harmonizome CLI")
        write_status("i", "Use --help to see available commands", "blue")
        click.echo()
        click.echo(ctx.get_help())
        ctx.exit(1)


@main.command("list-datasets")
def list_datasets_command() -> None:
    """List all available datasets."""
    list_datasets()


@main.command("get-entity")
@click.argument("entity_type", type=ENTITY_TYPE_CHOICES)
@click.argument("name")
def get_entity_command(entity_type: str, name: str) -> None:
    """Get information about a specific entity."""
    get_entity_info(entity_type.lower(), name)


@main.command("download")
@click.argument("datasets", nargs=-1, required=True)
@click.option(
    "--output-dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=str),
    help="Output directory for downloads. Defaults to the current directory.",
)
def download_command(datasets: tuple[str, ...], output_dir: Optional[str]) -> None:
    """Download datasets."""
    download_datasets(list(datasets), output_dir)


@main.command("functional-associations")
@click.argument("gene_symbol")
@click.option(
    "--datasets",
    multiple=True,
    help="Specific datasets to search. Defaults to all datasets.",
)
@click.option(
    "--output-file",
    type=click.Path(dir_okay=False, path_type=str),
    help="Save results to a JSON file.",
)
def functional_associations_command(
    gene_symbol: str,
    datasets: tuple[str, ...],
    output_file: Optional[str],
) -> None:
    """Get functional associations for a gene."""
    requested_datasets = list(datasets) if datasets else None
    get_functional_associations(
        gene_symbol,
        datasets=requested_datasets,
        output_file=output_file,
    )


if __name__ == "__main__":
    main()
