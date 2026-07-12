#!/usr/bin/env python3
"""Command-line interface for the Harmonizome package."""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any, Optional

from harmonizome import Harmonizome


# ANSI color codes for terminal output
class Colors:
    HEADER = "\033[95m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


def print_header(text: str) -> None:
    """Print a formatted header."""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}")
    print(f"  {text}")
    print(f"{'='*60}{Colors.END}\n")


def print_success(text: str) -> None:
    """Print success message in green."""
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")


def print_warning(text: str) -> None:
    """Print warning message in yellow."""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")


def print_error(text: str) -> None:
    """Print error message in red."""
    print(f"{Colors.RED}✗ {text}{Colors.END}")


def print_info(text: str) -> None:
    """Print info message in blue."""
    print(f"{Colors.BLUE}ℹ {text}{Colors.END}")


def print_dataset_item(index: int, dataset: str, short_code: str) -> None:
    """Print a formatted dataset item."""
    print(f"{Colors.CYAN}{index:3d}.{Colors.END} {Colors.BOLD}{dataset}{Colors.END}")
    print(f"       {Colors.YELLOW}Short code:{Colors.END} {short_code}")
    print()


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def list_datasets() -> None:
    """List all available datasets."""
    from harmonizome.harmonizome import DATASET_TO_PATH

    print_header("Available Harmonizome Datasets")
    print_info(f"Total datasets: {len(Harmonizome.DATASETS)}")
    print_warning(
        "Use the full dataset name for downloading (e.g., 'ENCODE Transcription Factor Binding Site Profiles')"
    )
    print()

    for i, dataset in enumerate(Harmonizome.DATASETS, 1):
        short_code = DATASET_TO_PATH.get(dataset, "N/A")
        print_dataset_item(i, dataset, short_code)


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
        print_header(f"{entity_type.title()}: {name}")

        # Define the order and labels for common fields
        field_order = {
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

        # Get the field order for this entity type, or use all keys
        ordered_fields = field_order.get(entity_type, list(info.keys()))

        # Add any missing fields to the end
        for key in info.keys():
            if key not in ordered_fields:
                ordered_fields.append(key)

        for key in ordered_fields:
            if key in info:
                value = info[key]
                formatted_value = format_value(key, value)

                # Format the key name for display
                display_key = key.replace("_", " ").title()
                if key == "ncbiEntrezGeneId":
                    display_key = "NCBI Entrez Gene ID"
                elif key == "ncbiEntrezGeneUrl":
                    display_key = "NCBI Entrez Gene URL"
                elif key == "uniprotId":
                    display_key = "UniProt ID"
                elif key == "uniprotUrl":
                    display_key = "UniProt URL"

                print(f"{Colors.BOLD}{display_key}:{Colors.END}")
                if "\n" in formatted_value:
                    print(f"{formatted_value}")
                else:
                    print(f"  {formatted_value}")
                print()

    except Exception as e:
        print_error(f"Error getting {entity_type} '{name}': {e}")
        sys.exit(1)


def find_dataset_by_partial_name(partial_name: str) -> list[str]:
    """Find datasets that contain the partial name."""
    from harmonizome.harmonizome import DATASET_TO_PATH

    matching_datasets = []
    for dataset_name in DATASET_TO_PATH.keys():
        if partial_name.upper() in dataset_name.upper():
            matching_datasets.append(dataset_name)
    return matching_datasets


def download_datasets(
    datasets: list[str], output_dir: Optional[str] = None
) -> None:
    """Download specified datasets."""
    print_header("Harmonizome Dataset Download")

    if output_dir:
        print_info(f"Output directory: {output_dir}")
        Path(output_dir).mkdir(exist_ok=True)
        original_working_directory = Path.cwd()
        os.chdir(output_dir)
    else:
        print_info("Output directory: current directory")
        original_working_directory = None

    # Check if any datasets need to be expanded (e.g., "ENCODE" -> all ENCODE datasets)
    expanded_datasets = []
    for dataset in datasets:
        if dataset.upper() in ["ENCODE", "GTEX", "MSIGDB"]:
            # Find all datasets containing this name
            matching = find_dataset_by_partial_name(dataset)
            if matching:
                print_success(f"Found {len(matching)} datasets matching '{dataset}':")
                for match in matching:
                    print(f"  {Colors.CYAN}•{Colors.END} {match}")
                expanded_datasets.extend(matching)
            else:
                expanded_datasets.append(dataset)
        else:
            expanded_datasets.append(dataset)

    print_info(f"Total datasets to download: {len(expanded_datasets)}")
    print()

    try:
        for filename in Harmonizome.download(expanded_datasets):
            print_success(f"Downloaded: {filename}")
    except KeyboardInterrupt:
        print_warning("Download interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print_error(f"Error downloading datasets: {e}")
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
    print_header(f"Functional Associations for {gene_symbol}")

    if datasets:
        print_info(f"Datasets: {', '.join(datasets)}")
    else:
        print_info("Datasets: All available")

    try:
        # Get functional associations
        results = Harmonizome.get_gene_functional_annotations(gene_symbol, datasets)

        # Display results
        gene_info = results["gene_info"]
        func_assoc = results["functional_associations"]

        print(f"\n{Colors.BOLD}Gene Information:{Colors.END}")
        print(f"  Symbol: {gene_info['symbol']}")
        print(f"  Name: {gene_info['name']}")
        print(f"  NCBI ID: {gene_info['ncbi_id']}")
        print(f"  Description: {gene_info['description'][:200]}...")

        print(f"\n{Colors.BOLD}Functional Associations Summary:{Colors.END}")
        print(f"  Total datasets: {func_assoc['total_datasets']}")
        print(f"  Total associations: {func_assoc['total_associations']}")
        print(f"  Increased associations: {func_assoc['total_increased']}")
        print(f"  Decreased associations: {func_assoc['total_decreased']}")

        # Show dataset details
        print(f"\n{Colors.BOLD}Dataset Details:{Colors.END}")
        for dataset in func_assoc["datasets"][:5]:  # Show first 5 datasets
            print(f"\n{Colors.CYAN}Dataset: {dataset['dataset']}{Colors.END}")
            print(f"  Summary: {dataset['summary'][:100]}...")

            for assoc_group in dataset["associations"]:
                print(f"  {assoc_group['description']}:")
                for item in assoc_group["items"][:3]:  # Show first 3
                    print(f"    {item['name']} [{item['score']:.5f}]")
                if len(assoc_group["items"]) > 3:
                    print(f"    ... and {len(assoc_group['items']) - 3} more")

        if len(func_assoc["datasets"]) > 5:
            print(
                f"\n{Colors.YELLOW}... and {len(func_assoc['datasets']) - 5} more datasets{Colors.END}"
            )

        # Save to file if requested
        if output_file:
            import json

            with open(output_file, "w") as f:
                json.dump(results, f, indent=2)
            print_success(f"Results saved to: {output_file}")

    except Exception as e:
        print_error(f"Error getting functional associations: {e}")
        sys.exit(1)


def main() -> None:
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description=f"{Colors.BOLD}Harmonizome CLI{Colors.END} - Access to harmonized datasets of genes and proteins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
{Colors.BOLD}Examples:{Colors.END}
  %(prog)s list-datasets
  %(prog)s get-entity gene BRCA1
  %(prog)s download ENCODE GTEx
  %(prog)s download --output-dir ./data ENCODE
  %(prog)s functional-associations STAT3
  %(prog)s functional-associations BRCA1 --datasets ENCODE GTEx
  %(prog)s functional-associations STAT3 --use-download --output-file results.json
        """,
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # List datasets command
    subparsers.add_parser("list-datasets", help="List all available datasets")

    # Get entity command
    entity_parser = subparsers.add_parser(
        "get-entity", help="Get information about a specific entity"
    )
    entity_parser.add_argument(
        "entity_type",
        choices=["gene", "gene_set", "attribute", "dataset", "protein", "resource"],
        help="Type of entity to query",
    )
    entity_parser.add_argument("name", help="Name of the entity")

    # Download command
    download_parser = subparsers.add_parser("download", help="Download datasets")
    download_parser.add_argument(
        "datasets", nargs="+", help="Names of datasets to download"
    )
    download_parser.add_argument(
        "--output-dir",
        help="Output directory for downloads (default: current directory)",
    )

    # Functional associations command
    assoc_parser = subparsers.add_parser(
        "functional-associations", help="Get functional associations for a gene"
    )
    assoc_parser.add_argument("gene_symbol", help="Gene symbol (e.g., STAT3, BRCA1)")
    assoc_parser.add_argument(
        "--datasets",
        nargs="+",
        help="Specific datasets to search (default: all datasets)",
    )

    assoc_parser.add_argument("--output-file", help="Save results to JSON file")

    args = parser.parse_args()

    if not args.command:
        print_header("Welcome to Harmonizome CLI")
        print_info("Use --help to see available commands")
        print()
        parser.print_help()
        sys.exit(1)

    setup_logging(args.verbose)

    if args.command == "list-datasets":
        list_datasets()
    elif args.command == "get-entity":
        get_entity_info(args.entity_type, args.name)
    elif args.command == "download":
        download_datasets(args.datasets, args.output_dir)
    elif args.command == "functional-associations":
        get_functional_associations(
            args.gene_symbol, datasets=args.datasets, output_file=args.output_file
        )


if __name__ == "__main__":
    main()
