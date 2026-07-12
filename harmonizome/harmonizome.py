"""Class for reading, parsing, and downloading data from the Harmonizome API."""

import gzip
import json
import logging
import ssl
from collections.abc import Iterator, Set
from io import BytesIO
from pathlib import Path
from typing import Any, BinaryIO, Optional, Union
from urllib.error import HTTPError, URLError
from urllib.parse import quote_plus
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix

from .utils import cache_to_file

logger = logging.getLogger(__name__)

GENE_MATRIX_DOWNLOADS = [
    "gene_attribute_matrix.txt.gz",
    "attribute_list_entries.txt.gz",
]
DEFAULT_GENE_INFO = {
    "name": "Unknown",
    "description": "No description available",
}
DEFAULT_CONFIG = {
    "downloads": [
        "gene_attribute_matrix.txt.gz",
        "gene_list_terms.txt.gz",
        "attribute_list_entries.txt.gz",
    ],
    "datasets": {"ENCODE": "encode", "GTEx": "gtex"},
}
_CONFIG: Optional[dict[str, Any]] = None
DOWNLOADS: list[str] = []
DATASET_TO_PATH: dict[str, str] = {}


class Enum(set):
    """Simple Enum shim used by the historical public API."""

    def __getattr__(self, name: str) -> str:
        if name in self:
            return name
        raise AttributeError


# The entity types supported by the Harmonizome API.
class Entity(Enum):

    DATASET = "dataset"
    GENE = "gene"
    GENE_SET = "gene_set"
    ATTRIBUTE = "attribute"
    GENE_FAMILY = "gene_family"
    NAMING_AUTHORITY = "naming_authority"
    PROTEIN = "protein"
    RESOURCE = "resource"


def json_from_url(url: str) -> dict[str, Any]:
    """Returns API response after decoding and loading JSON.

    Note: Uses unverified SSL context due to certificate verification issues
    with the Harmonizome API. This is a pragmatic solution for accessing
    publicly available genomic data. For production use with sensitive data,
    consider implementing proper certificate verification.
    """
    # Use unverified context to handle SSL certificate issues
    # This is necessary due to the API's certificate configuration
    context = ssl._create_unverified_context()

    try:
        response = urlopen(url, context=context)
        data = response.read()

        # Try UTF-8 first, then fallback to latin-1 for problematic responses
        try:
            decoded_data = data.decode("utf-8")
        except UnicodeDecodeError:
            logger.warning(f"UTF-8 decode failed for {url}, trying latin-1")
            decoded_data = data.decode("latin-1")

        return json.loads(decoded_data)
    except (HTTPError, URLError, json.JSONDecodeError) as e:
        logger.error(f"Failed to fetch data from {url}: {e}")
        raise


def download_from_url(url: str) -> BinaryIO:
    """Downloads a file from URL with SSL workaround.

    Note: Uses unverified SSL context due to certificate verification issues
    with the Harmonizome API. This is a pragmatic solution for accessing
    publicly available genomic data.
    """
    # Use unverified context to handle SSL certificate issues
    context = ssl._create_unverified_context()

    try:
        response = urlopen(url, context=context)
        return response
    except (HTTPError, URLError) as e:
        logger.error(f"Failed to download from {url}: {e}")
        raise


VERSION = "1.0.1"
API_URL = "https://maayanlab.cloud/Harmonizome/api"
DOWNLOAD_URL = "https://maayanlab.cloud/static/hdfs/harmonizome/data"


def _load_config() -> dict[str, Any]:
    """Load configuration from API with SSL fallback."""
    try:
        return json_from_url(f"{API_URL}/dark/script_config")
    except (HTTPError, URLError, json.JSONDecodeError) as e:
        logger.error(f"Failed to load configuration from API: {e}")
        return DEFAULT_CONFIG.copy()


def _ensure_config_loaded() -> dict[str, Any]:
    """Populate module-level config only when a caller first needs it."""
    global _CONFIG, DOWNLOADS, DATASET_TO_PATH

    if _CONFIG is None:
        _CONFIG = _load_config()

    if not DOWNLOADS:
        DOWNLOADS = list(_CONFIG.get("downloads", []))
    if not DATASET_TO_PATH:
        DATASET_TO_PATH = dict(_CONFIG.get("datasets", {}))

    return _CONFIG


def _get_downloads() -> list[str]:
    """Return the known dataset download targets."""
    _ensure_config_loaded()
    return DOWNLOADS


def _get_dataset_to_path() -> dict[str, str]:
    """Return the API dataset-name to path mapping."""
    _ensure_config_loaded()
    return DATASET_TO_PATH


class _LazyDatasetNames(Set[str]):
    """Resolve the dataset-name collection only when it is first iterated."""

    def __contains__(self, value: object) -> bool:
        return value in _get_dataset_to_path()

    def __iter__(self) -> Iterator[str]:
        return iter(_get_dataset_to_path())

    def __len__(self) -> int:
        return len(_get_dataset_to_path())


class GeneData:
    """Container for one gene record and its Harmonizome associations."""

    def __init__(self, gene_info: dict[str, Any]) -> None:
        self.gene_info = gene_info
        self.associations = gene_info.get("associations", [])

    @staticmethod
    def _parse_gene_set(assoc: dict[str, Any]) -> tuple[str, str]:
        """Extract gene-set and dataset names from either API payload shape."""
        gene_set = assoc.get("geneSet", {})
        gene_set_full_name = gene_set.get("name", "")
        dataset_name = gene_set.get("dataset", {}).get("name")

        if "/" in gene_set_full_name:
            gene_set_name, parsed_dataset_name = gene_set_full_name.split("/", 1)
        else:
            gene_set_name = gene_set_full_name
            parsed_dataset_name = ""

        return gene_set_name, dataset_name or parsed_dataset_name

    @staticmethod
    def _association_to_row(assoc: dict[str, Any]) -> dict[str, Any]:
        """Normalize one API association into a tabular row."""
        gene_set_name, dataset_name = GeneData._parse_gene_set(assoc)

        return {
            "gene_set": gene_set_name,
            "dataset": dataset_name,
            "thresholdValue": assoc.get("thresholdValue"),
            "standardizedValue": assoc.get("standardizedValue"),
        }

    def get_associations(self, dataset: Optional[str] = None) -> list[dict[str, Any]]:
        if dataset is None:
            return self.associations
        return [
            assoc
            for assoc in self.associations
            if self._parse_gene_set(assoc)[1] == dataset
        ]

    def save(
        self, path: str, format: str = "json", dataset: Optional[str] = None
    ) -> None:
        output_path = Path(path)
        rows = [
            self._association_to_row(assoc) for assoc in self.get_associations(dataset)
        ]
        if format == "json":
            with output_path.open("w", encoding="utf-8") as file_handle:
                json.dump(rows, file_handle, indent=2)
        elif format == "csv":
            pd.DataFrame(rows).to_csv(output_path, index=False)
        else:
            raise ValueError("format must be 'json' or 'csv'")

    def to_dataframe(self, dataset: Optional[str] = None) -> pd.DataFrame:
        """
        Return associations as a pandas DataFrame, with columns:
        'gene_set', 'dataset', 'thresholdValue', 'standardizedValue'.
        Optionally filter by dataset name.
        """
        rows = [
            self._association_to_row(assoc) for assoc in self.get_associations(dataset)
        ]
        return pd.DataFrame(rows)


class Harmonizome:

    __version__ = VERSION
    DATASETS = _LazyDatasetNames()

    @staticmethod
    def _build_association_group(
        associations: list[dict[str, Any]], association_type: str
    ) -> dict[str, Any]:
        """Format one directional association group for output."""
        return {
            "type": association_type,
            "count": len(associations),
            "description": f"{len(associations)} {association_type} fitness associations",
            "items": [
                {"name": assoc["name"], "score": assoc["score"]}
                for assoc in associations
            ],
        }

    @classmethod
    def _format_functional_associations_from_grouped_data(
        cls, grouped_data: dict[str, Any]
    ) -> dict[str, Any]:
        """Convert grouped dataset associations into the public response format."""
        functional_associations = []
        total_associations = 0
        total_increased = 0
        total_decreased = 0

        for dataset_name, dataset_data in grouped_data.items():
            dataset_entry = {
                "dataset": dataset_name,
                "summary": dataset_data["summary"],
                "associations": [],
            }

            increased_associations = dataset_data.get(
                "increased", dataset_data.get("increased_associations", [])
            )
            decreased_associations = dataset_data.get(
                "decreased", dataset_data.get("decreased_associations", [])
            )

            if increased_associations:
                dataset_entry["associations"].append(
                    cls._build_association_group(increased_associations, "increased")
                )
                total_increased += len(increased_associations)

            if decreased_associations:
                dataset_entry["associations"].append(
                    cls._build_association_group(decreased_associations, "decreased")
                )
                total_decreased += len(decreased_associations)

            if dataset_entry["associations"]:
                functional_associations.append(dataset_entry)
                total_associations += len(increased_associations) + len(
                    decreased_associations
                )

        return {
            "total_datasets": len(functional_associations),
            "total_associations": total_associations,
            "total_increased": total_increased,
            "total_decreased": total_decreased,
            "datasets": functional_associations,
        }

    @staticmethod
    def _format_gene_info(
        gene_info: dict[str, Any], gene_symbol: str
    ) -> dict[str, Any]:
        """Normalize gene metadata into the public response shape."""
        return {
            "symbol": gene_info.get("symbol", gene_symbol),
            "name": gene_info.get("name", DEFAULT_GENE_INFO["name"]),
            "description": gene_info.get(
                "description", DEFAULT_GENE_INFO["description"]
            ),
            "ncbi_id": gene_info.get("ncbiEntrezGeneId", "Unknown"),
        }

    @staticmethod
    def _build_fallback_gene_info(gene_symbol: str) -> dict[str, Any]:
        """Return a stable fallback payload when gene metadata lookup fails."""
        return {"symbol": gene_symbol, **DEFAULT_GENE_INFO}

    @classmethod
    def _build_functional_associations_response(
        cls,
        gene_info: dict[str, Any],
        gene_symbol: str,
        grouped_annotations: dict[str, Any],
    ) -> dict[str, Any]:
        """Assemble the public functional-associations response."""
        return {
            "gene_info": cls._format_gene_info(gene_info, gene_symbol),
            "functional_associations": cls._format_functional_associations_from_grouped_data(
                grouped_annotations
            ),
        }

    @staticmethod
    def _group_api_associations_by_dataset(
        gene_symbol: str,
        associations: list[dict[str, Any]],
        datasets: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """Group API association records by dataset and direction."""
        dataset_associations: dict[str, Any] = {}

        for assoc in associations:
            gene_set = assoc.get("geneSet", {})
            gene_set_name, dataset_name = GeneData._parse_gene_set(assoc)
            dataset_name = dataset_name or "Unknown Dataset"

            if datasets and dataset_name not in datasets:
                continue

            if dataset_name not in dataset_associations:
                dataset_associations[dataset_name] = {
                    "increased": [],
                    "decreased": [],
                    "summary": f"Associations for {gene_symbol} in {dataset_name}",
                }

            threshold_value = assoc.get("thresholdValue", 0)
            standardized_value = assoc.get("standardizedValue", 0)
            score = standardized_value if standardized_value != 0 else threshold_value

            association_item = {
                "name": gene_set_name or gene_set.get("name", "Unknown"),
                "score": score,
                "gene_set_id": gene_set.get("id", ""),
                "dataset": dataset_name,
            }

            if score > 0:
                dataset_associations[dataset_name]["increased"].append(association_item)
            elif score < 0:
                dataset_associations[dataset_name]["decreased"].append(association_item)

        return dataset_associations

    @classmethod
    def get(
        cls, entity: str, name: Optional[str] = None, start_at: Optional[int] = None
    ) -> dict[str, Any]:
        """Returns a single entity or a list, depending on if a name is
        provided. If no name is provided and start_at is specified, returns a
        list starting at that cursor position.
        """
        if name:
            name = quote_plus(name)
            return _get_by_name(entity, name)
        if isinstance(start_at, int):
            return _get_with_cursor(entity, start_at)
        return json_from_url(f"{API_URL}/{VERSION}/{entity}")

    @classmethod
    def next(cls, response: dict[str, Any]) -> dict[str, Any]:
        """Returns the next set of entities based on a previous API response."""
        start_at = _get_next(response)
        entity = _get_entity(response)
        return cls.get(entity=entity, start_at=start_at)

    @classmethod
    def download(
        cls,
        datasets: Optional[list[str]] = None,
        what: Optional[list[str]] = None,
    ) -> Iterator[str]:
        """For each dataset, creates a directory and downloads files into it."""
        # Why not check `if not datasets`? Because in principle, a user could
        # call `download([])`, which should download nothing, not everything.
        # Why might they do this? Imagine that the list of datasets is
        # dynamically generated in another user script.
        if datasets is None:
            datasets = cls.DATASETS
            warning = (
                "Warning: You are going to download all Harmonizome "
                "data. This is roughly 30GB. Do you accept?\n(Y/N) "
            )
            resp = input(warning)
            if resp.lower() != "y":
                return

        dataset_to_path = _get_dataset_to_path()
        download_targets = what if what is not None else _get_downloads()

        for dataset in datasets:
            if dataset not in cls.DATASETS:
                msg = (
                    f'"{dataset}" is not a valid dataset name. Check the `DATASETS` '
                    "property for a complete list of names."
                )
                raise AttributeError(msg)
            dataset_dir = Path(dataset)
            dataset_dir.mkdir(exist_ok=True)

            for dl in download_targets:
                path = dataset_to_path[dataset]
                url = f"{DOWNLOAD_URL}/{path}/{dl}"

                try:
                    response = download_from_url(url)
                except HTTPError as e:
                    # Not every dataset has all downloads.
                    logger.warning("Skipping %s: HTTP Error %s", dl, e.code)
                    continue
                except URLError as e:
                    logger.warning("Skipping %s: %s", dl, e)
                    continue

                file_path = dataset_dir / dl.replace(".gz", "")

                if response.code != 200:
                    logger.warning("Skipping %s: HTTP status %s", dl, response.code)
                    continue

                if file_path.is_file():
                    logger.info("Using cached `%s`", file_path)
                else:
                    logger.info("Downloading `%s`", file_path)
                    _download_and_decompress_file(response, file_path)

                yield str(file_path)

    @classmethod
    def download_df(
        cls,
        datasets: Optional[list[str]] = None,
        what: Optional[list[str]] = None,
        sparse: bool = False,
        **kwargs: Any,
    ) -> Iterator[pd.DataFrame]:
        for file in cls.download(datasets, what):
            if sparse:
                yield _read_as_sparse_dataframe(file, **kwargs)
            else:
                yield _read_as_dataframe(file, **kwargs)

    @classmethod
    def get_gene_functional_annotations(
        cls, gene_symbol: str, datasets: Optional[list[str]] = None
    ) -> dict[str, Any]:
        """Get functional annotations for a gene using the Harmonizome API.

        This method uses the API directly without downloading any files.
        It leverages the showAssociations=true parameter to get functional
        associations for a gene.

        Args:
            gene_symbol: Gene symbol (e.g., 'BRCA1', 'STAT3')
            datasets: List of dataset names to search. If None, searches all datasets.

        Returns:
            Dictionary formatted like Harmonizome web interface with:
            - gene_info: Basic gene information
            - functional_associations: List of datasets with associations
        """
        # Get gene info with associations
        try:
            gene_info = cls.get_gene_with_associations(gene_symbol)
        except (HTTPError, URLError, json.JSONDecodeError) as e:
            logger.warning(f"Could not get gene info for {gene_symbol}: {e}")
            gene_info = cls._build_fallback_gene_info(gene_symbol)

        dataset_associations = cls._group_api_associations_by_dataset(
            gene_symbol=gene_symbol,
            associations=gene_info.get("associations", []),
            datasets=datasets,
        )
        return cls._build_functional_associations_response(
            gene_info=gene_info,
            gene_symbol=gene_symbol,
            grouped_annotations=dataset_associations,
        )

    @classmethod
    def get_gene_with_associations(cls, gene_symbol: str) -> dict[str, Any]:
        """Get gene information with associations using the API.

        This uses the showAssociations=true parameter to get functional
        associations directly from the API without downloading files.

        Args:
            gene_symbol: Gene symbol

        Returns:
            Dictionary with gene info and associations
        """
        name = quote_plus(gene_symbol)
        return json_from_url(f"{API_URL}/{VERSION}/gene/{name}?showAssociations=true")

    @classmethod
    def download_gene_functional_annotations(
        cls, gene_symbol: str, datasets: Optional[list[str]] = None
    ) -> dict[str, Any]:
        """Download and get all functional annotations for a gene across specified datasets.

        This method downloads the necessary data files and extracts associations
        for the specified gene from the gene-attribute matrices.

        Args:
            gene_symbol: Gene symbol (e.g., 'BRCA1', 'STAT3')
            datasets: List of dataset names to search. If None, searches all datasets.

        Returns:
            Dictionary with dataset names as keys and annotation data as values.
        """
        if datasets is None:
            datasets = list(cls.DATASETS)

        results = {}

        for dataset in datasets:
            dataset_annotations = cls._get_gene_dataset_annotations_from_files(
                gene_symbol, dataset
            )
            if dataset_annotations:
                results[dataset] = dataset_annotations

        return results

    @classmethod
    def _download_dataset_annotation_files(
        cls, dataset: str
    ) -> tuple[Optional[str], Optional[str]]:
        """Download the files needed to recover one dataset's annotations."""
        gene_matrix_file = None
        attribute_list_file = None

        for filename in cls.download([dataset], GENE_MATRIX_DOWNLOADS):
            if "gene_attribute_matrix.txt" in filename:
                gene_matrix_file = filename
            elif "attribute_list_entries.txt" in filename:
                attribute_list_file = filename

        return gene_matrix_file, attribute_list_file

    @staticmethod
    def _get_gene_row_from_matrix(
        matrix_df: pd.DataFrame, gene_symbol: str, dataset: str
    ) -> Optional[pd.Series]:
        """Locate the matrix row that corresponds to the requested gene."""
        for index_value in matrix_df.index:
            try:
                gene_data = json.loads(index_value)
            except (json.JSONDecodeError, TypeError):
                continue

            if (
                isinstance(gene_data, list)
                and gene_data
                and gene_data[0] == gene_symbol
            ):
                return matrix_df.loc[index_value]

        logger.debug(f"Gene '{gene_symbol}' not found in dataset '{dataset}'")
        return None

    @staticmethod
    def _load_attribute_names(attribute_list_file: Optional[str]) -> dict[str, str]:
        """Map attribute identifiers to readable names when the sidecar file exists."""
        if attribute_list_file is None:
            return {}

        try:
            attribute_frame = pd.read_csv(
                attribute_list_file, sep="\t", encoding="latin-1"
            )
        except (
            FileNotFoundError,
            pd.errors.EmptyDataError,
            pd.errors.ParserError,
            UnicodeDecodeError,
        ) as e:
            logger.debug(f"Could not read attribute list '{attribute_list_file}': {e}")
            return {}

        if len(attribute_frame.columns) < 2:
            return {}

        attribute_names: dict[str, str] = {}
        for _, row in attribute_frame.iterrows():
            attribute_id = str(row.iloc[0])
            try:
                attribute_data = json.loads(row.iloc[0])
            except (json.JSONDecodeError, TypeError):
                attribute_names[attribute_id] = str(row.iloc[1])
                continue

            attribute_name = (
                attribute_data[0]
                if isinstance(attribute_data, list) and attribute_data
                else attribute_id
            )
            attribute_names[attribute_id] = str(attribute_name)

        return attribute_names

    @staticmethod
    def _build_annotation_summary(
        summary: str, gene_row: pd.Series, attribute_names: dict[str, str]
    ) -> Optional[dict[str, Any]]:
        """Convert one matrix row into grouped positive and negative associations."""
        associations: list[dict[str, Any]] = []
        increased_associations: list[dict[str, Any]] = []
        decreased_associations: list[dict[str, Any]] = []

        for attr_id, score in gene_row.items():
            if pd.isna(score) or score == 0:
                continue

            try:
                score_float = float(score)
            except (ValueError, TypeError):
                continue

            association = {
                "name": attribute_names.get(attr_id, attr_id),
                "score": score_float,
                "attribute_id": attr_id,
            }
            associations.append(association)

            if score_float > 0:
                increased_associations.append(association)
            elif score_float < 0:
                decreased_associations.append(association)

        if not associations:
            return None

        return {
            "summary": summary,
            "associations": associations,
            "increased_associations": increased_associations,
            "decreased_associations": decreased_associations,
            "total_associations": len(associations),
            "increased_count": len(increased_associations),
            "decreased_count": len(decreased_associations),
        }

    @classmethod
    def _get_gene_dataset_annotations_from_files(
        cls, gene_symbol: str, dataset: str
    ) -> Optional[dict[str, Any]]:
        """Get functional annotations for a gene in a specific dataset using downloaded files.

        Args:
            gene_symbol: Gene symbol
            dataset: Dataset name

        Returns:
            Dictionary with annotation data for the dataset
        """
        dataset_info = cls.get("dataset", dataset)
        summary = dataset_info.get(
            "description", f"Associations for {gene_symbol} in {dataset}"
        )

        gene_matrix_file, attribute_list_file = cls._download_dataset_annotation_files(
            dataset
        )
        if not gene_matrix_file:
            logger.debug(
                f"Could not download gene-attribute matrix for dataset '{dataset}'"
            )
            return None

        matrix_df = _read_as_dataframe(gene_matrix_file)
        gene_row = cls._get_gene_row_from_matrix(matrix_df, gene_symbol, dataset)
        if gene_row is None:
            return None

        attribute_names = cls._load_attribute_names(attribute_list_file)
        return cls._build_annotation_summary(summary, gene_row, attribute_names)

    @classmethod
    def get_gene_associations_summary(
        cls,
        gene_symbol: str,
        datasets: Optional[list[str]] = None,
        use_download: bool = False,
    ) -> dict[str, Any]:
        """Get a summary of all functional associations for a gene.

        Args:
            gene_symbol: Gene symbol
            datasets: List of dataset names to search. If None, searches all datasets.
            use_download: If True, downloads files to get associations. If False, uses API.

        Returns:
            Dictionary with summary statistics and dataset breakdown
        """
        if use_download:
            annotations = cls.download_gene_functional_annotations(
                gene_symbol, datasets
            )
        else:
            annotations = cls.get_gene_functional_annotations(gene_symbol, datasets)

        total_datasets = len(annotations)
        total_associations = sum(d["total_associations"] for d in annotations.values())
        total_increased = sum(d["increased_count"] for d in annotations.values())
        total_decreased = sum(d["decreased_count"] for d in annotations.values())

        return {
            "gene_symbol": gene_symbol,
            "total_datasets": total_datasets,
            "total_associations": total_associations,
            "total_increased_associations": total_increased,
            "total_decreased_associations": total_decreased,
            "datasets": annotations,
        }

    @classmethod
    def get_gene_functional_associations_formatted(
        cls,
        gene_symbol: str,
        datasets: Optional[list[str]] = None,
        use_download: bool = False,
    ) -> dict[str, Any]:
        """Get functional associations for a gene in Harmonizome web interface format.

        This method returns data structured exactly like the Harmonizome web interface,
        with datasets, summaries, and categorized associations with scores.

        Args:
            gene_symbol: Gene symbol (e.g., 'STAT3', 'BRCA1')
            datasets: List of dataset names to search. If None, searches all datasets.
            use_download: If True, downloads files to get associations. If False, uses API.

        Returns:
            Dictionary formatted like Harmonizome web interface with:
            - gene_info: Basic gene information
            - functional_associations: List of datasets with associations
        """
        try:
            gene_info = cls.get("gene", gene_symbol)
        except (HTTPError, URLError, json.JSONDecodeError) as e:
            logger.warning(f"Could not get gene info for {gene_symbol}: {e}")
            gene_info = cls._build_fallback_gene_info(gene_symbol)

        if use_download:
            annotations = cls.download_gene_functional_annotations(
                gene_symbol, datasets
            )
        else:
            return cls.get_gene_functional_annotations(gene_symbol, datasets)
        return cls._build_functional_associations_response(
            gene_info=gene_info,
            gene_symbol=gene_symbol,
            grouped_annotations=annotations,
        )

    @classmethod
    def get_gene_data(cls, gene_symbol: str, use_cache: bool = False) -> GeneData:
        if use_cache:
            return GeneData(cls._get_gene_with_associations_cached(gene_symbol))
        else:
            return GeneData(cls.get_gene_with_associations(gene_symbol))

    @staticmethod
    @cache_to_file
    def _get_gene_with_associations_cached(
        gene_symbol: str,
    ) -> dict[str, Any]:
        return Harmonizome.get_gene_with_associations(gene_symbol)


def _get_with_cursor(entity: str, start_at: int) -> dict[str, Any]:
    """Returns a list of entities based on cursor position."""
    return json_from_url(f"{API_URL}/{VERSION}/{entity}?cursor={start_at}")


def _get_by_name(entity: str, name: str) -> dict[str, Any]:
    """Returns a single entity based on name."""
    return json_from_url(f"{API_URL}/{VERSION}/{entity}/{name}")


def _get_entity(response: dict[str, Any]) -> str:
    """Returns the entity from an API response."""
    path = response["next"].split("?")[0]
    return path.split("/")[3]


def _get_next(response: dict[str, Any]) -> Optional[int]:
    """Returns the next property from an API response."""
    if response["next"]:
        return int(response["next"].split("=")[1])
    return None


def _download_and_decompress_file(
    response: BinaryIO, filename: Union[Path, str]
) -> None:
    """Write one gzip-compressed API download to its decompressed on-disk form."""
    compressed_file = BytesIO(response.read())
    decompressed_file = gzip.GzipFile(fileobj=compressed_file)
    output_path = Path(filename)
    with output_path.open("wb") as outfile:
        outfile.write(decompressed_file.read())


def _getfshape(
    fn: str,
    row_sep: str = "\n",
    col_sep: str = "\t",
    open_args: Optional[dict[str, Any]] = None,
) -> tuple[int, int]:
    """Fast and efficient way of finding row/col height of file"""
    open_kwargs = {} if open_args is None else dict(open_args)
    with Path(fn).open("r", newline=row_sep, **open_kwargs) as file_handle:
        col_size = file_handle.readline().count(col_sep) + 1
        row_size = sum(1 for line in file_handle) + 1
        return row_size, col_size


def _parse(
    fn: str,
    column_size: int = 3,
    index_size: int = 3,
    shape: Optional[tuple[int, int]] = None,
    index_fmt: Any = None,
    data_fmt: Any = None,
    index_dtype: Any = None,
    data_dtype: Any = None,
    col_sep: str = "\t",
    row_sep: str = "\n",
    open_args: Optional[dict[str, Any]] = None,
) -> tuple[Any, Any, Any, Any, Any]:
    """
    Smart(er) parser for processing matrix formats. Evaluate size and construct
     ndframes with the right size before parsing, this allows for more efficient
     loading of sparse dataframes as well. To obtain a sparse representation use:
         data_fmt=scipy.lil_matrix
    This only works if all of the data is of the same type, if it isn't a float
     use:
         data_dtype=np.float64

    Returns:
        (column_names, columns, index_names, index, data)
    """
    if index_fmt is None:
        index_fmt = np.ndarray
    if data_fmt is None:
        data_fmt = np.ndarray
    if index_dtype is None:
        index_dtype = object
    if data_dtype is None:
        data_dtype = np.float64

    if shape is not None:
        rows, cols = shape
    else:
        rows, cols = _getfshape(
            fn, row_sep=row_sep, col_sep=col_sep, open_args=open_args
        )

    columns = index_fmt((column_size, cols - index_size), dtype=index_dtype)
    index = index_fmt((rows - column_size, index_size), dtype=index_dtype)
    data = data_fmt((rows - column_size, cols - index_size), dtype=data_dtype)

    open_kwargs = {} if open_args is None else dict(open_args)

    with Path(fn).open("r", newline=row_sep, **open_kwargs) as file_handle:
        header = np.array(
            [next(file_handle).strip().split(col_sep) for _ in range(column_size)]
        )

        column_names = header[:column_size, index_size - 1]
        index_names = header[column_size - 1, :index_size]

        columns[:, :] = header[:column_size, index_size:]

        for ind, line in enumerate(file_handle):
            lh = line.strip().split(col_sep)
            index[ind, :] = lh[:index_size]
            data[ind, :] = lh[index_size:]

        return column_names, columns, index_names, index, data


def _parse_df(
    fn: str,
    sparse: bool = False,
    default_fill_value: Any = None,
    column_apply: Any = None,
    index_apply: Any = None,
    df_args: Optional[dict[str, Any]] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    data_fmt = lil_matrix if sparse else np.ndarray
    dataframe_kwargs = {} if df_args is None else dict(df_args)
    (
        column_names,
        columns,
        index_names,
        index,
        data,
    ) = _parse(fn, data_fmt=data_fmt, **kwargs)

    if column_apply is not None:
        column_names, columns = column_apply(column_names.T, columns.T)
    else:
        column_names, columns = (column_names.T, columns.T)

    if index_apply is not None:
        index_names, index = index_apply(index_names, index)

    index_values = pd.Index(data=index, name=str(index_names), dtype=object)
    column_values = pd.Index(data=columns, name=str(column_names), dtype=object)

    if sparse:
        fill_value = dataframe_kwargs.pop("default_fill_value", default_fill_value)
        dataframe = pd.DataFrame.sparse.from_spmatrix(
            data.tocsr(),
            index=index_values,
            columns=column_values,
        )
        if fill_value is not None:
            sparse_dtype = pd.SparseDtype(data.dtype, fill_value=fill_value)
            for column_name in dataframe.columns:
                dataframe[column_name] = dataframe[column_name].astype(sparse_dtype)
        if dataframe_kwargs:
            raise TypeError(
                "Unexpected DataFrame arguments for sparse parsing: "
                f"{sorted(dataframe_kwargs.keys())}"
            )
        return dataframe

    return pd.DataFrame(
        data=data,
        index=index_values,
        columns=column_values,
        **dataframe_kwargs,
    )


def _df_column_uniquify(df: pd.DataFrame) -> pd.DataFrame:
    """Append numeric suffixes to duplicate column names."""
    df_columns = df.columns
    new_columns = []
    for item in df_columns:
        counter = 0
        newitem = item
        while newitem in new_columns:
            counter += 1
            newitem = f"{item}_{counter}"
        new_columns.append(newitem)
    df.columns = new_columns
    return df


def _json_ind_no_slash(ind_names: Any, ind: Any) -> tuple[str, list[str]]:
    return (
        json.dumps([ind_name.replace("/", "|") for ind_name in ind_names]),
        [json.dumps([ii.replace("/", "|") for ii in i]) for i in ind],
    )


def _read_as_dataframe(fn: str) -> pd.DataFrame:
    """Load one Harmonizome text artifact into a pandas DataFrame."""
    if fn.endswith("gene_attribute_matrix.txt"):
        return _df_column_uniquify(
            _parse_df(
                fn,
                sparse=False,
                index_apply=_json_ind_no_slash,
                column_apply=_json_ind_no_slash,
                open_args=dict(encoding="latin-1"),
            )
        )
    elif fn.endswith("gene_list_terms.txt") or fn.endswith(
        "attribute_list_entries.txt"
    ):
        return pd.read_table(fn, encoding="latin-1", index_col=None)
    else:
        raise ValueError("Unable to parse this file into a dataframe.")


def _read_as_sparse_dataframe(fn: str, fill_value: int = 0) -> pd.DataFrame:
    """Load the gene-attribute matrix as a sparse pandas DataFrame."""
    if fn.endswith("gene_attribute_matrix.txt"):
        return _df_column_uniquify(
            _parse_df(
                fn,
                sparse=True,
                index_apply=_json_ind_no_slash,
                column_apply=_json_ind_no_slash,
                df_args=dict(default_fill_value=0),
                open_args=dict(encoding="latin-1"),
            )
        )
    else:
        raise ValueError("Unable to parse this file into a dataframe.")
