"""Class for reading, parsing, and downloading data from the Harmonizome API."""

import gzip
import json
import logging
import os
import ssl

# Support for both Python2.X and 3.X.
# -----------------------------------------------------------------------------
try:
    from io import BytesIO
    from urllib.error import HTTPError
    from urllib.parse import quote_plus
    from urllib.request import urlopen
except ImportError:
    from urllib import quote_plus

    from StringIO import StringIO as BytesIO
    from urllib2 import HTTPError, urlopen

try:
    input_shim = raw_input
except NameError:
    # If `raw_input` throws a `NameError`, the user is using Python 2.X.
    input_shim = input

import pandas as pd

from .utils import cache_to_file

# Enumerables and constants
# -----------------------------------------------------------------------------


class Enum(set):
    """Simple Enum shim since Python 2.X does not have them."""

    def __getattr__(self, name):
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


def json_from_url(url):
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
            logging.warning(f"UTF-8 decode failed for {url}, trying latin-1")
            decoded_data = data.decode("latin-1")

        return json.loads(decoded_data)
    except Exception as e:
        logging.error(f"Failed to fetch data from {url}: {e}")
        raise


def download_from_url(url):
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
    except Exception as e:
        logging.error(f"Failed to download from {url}: {e}")
        raise


VERSION = "1.0"
API_URL = "https://maayanlab.cloud/Harmonizome/api"
DOWNLOAD_URL = "https://maayanlab.cloud/static/hdfs/harmonizome/data"


# This config objects pulls the names of the datasets, their directories, and
# the possible downloads from the API. This allows us to add new datasets and
# downloads without breaking this file.
def _load_config():
    """Load configuration from API with SSL fallback."""
    try:
        config = json_from_url(API_URL + "/dark/script_config")
        return config
    except Exception as e:
        logging.error(f"Failed to load configuration from API: {e}")
        # Return minimal config to prevent import errors
        return {
            "downloads": [
                "gene_attribute_matrix.txt.gz",
                "gene_list_terms.txt.gz",
                "attribute_list_entries.txt.gz",
            ],
            "datasets": {"ENCODE": "encode", "GTEx": "gtex"},  # Minimal fallback
        }


config = _load_config()
DOWNLOADS = [x for x in config.get("downloads", [])]
DATASET_TO_PATH = config.get("datasets", {})


# Harmonizome class
# -----------------------------------------------------------------------------


class GeneData:
    def __init__(self, gene_info: dict):
        self.gene_info = gene_info
        self.associations = gene_info.get("associations", [])

    @staticmethod
    def _parse_gene_set(assoc: dict):
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
    def _association_to_row(assoc: dict) -> dict:
        """Normalize one API association into a tabular row."""
        gene_set_name, dataset_name = GeneData._parse_gene_set(assoc)

        return {
            "gene_set": gene_set_name,
            "dataset": dataset_name,
            "thresholdValue": assoc.get("thresholdValue"),
            "standardizedValue": assoc.get("standardizedValue"),
        }

    def get_associations(self, dataset: str = None):
        if dataset is None:
            return self.associations
        return [
            assoc
            for assoc in self.associations
            if self._parse_gene_set(assoc)[1] == dataset
        ]

    def save(self, path: str, format: str = "json", dataset: str = None):
        assocs = self.get_associations(dataset)
        rows = [self._association_to_row(assoc) for assoc in assocs]
        if format == "json":
            with open(path, "w") as f:
                json.dump(rows, f, indent=2)
        elif format == "csv":
            import pandas as pd

            pd.DataFrame(rows).to_csv(path, index=False)
        else:
            raise ValueError("format must be 'json' or 'csv'")

    def to_dataframe(self, dataset: str = None):
        """
        Return associations as a pandas DataFrame, with columns:
        'gene_set', 'dataset', 'thresholdValue', 'standardizedValue'.
        Optionally filter by dataset name.
        """
        assocs = self.get_associations(dataset)
        rows = [self._association_to_row(assoc) for assoc in assocs]
        import pandas as pd

        return pd.DataFrame(rows)


class Harmonizome(object):

    __version__ = VERSION
    DATASETS = DATASET_TO_PATH.keys()

    @staticmethod
    def _build_association_group(associations: list, association_type: str) -> dict:
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
        cls, grouped_data: dict
    ) -> dict:
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

    @classmethod
    def get(cls, entity, name=None, start_at=None):
        """Returns a single entity or a list, depending on if a name is
        provided. If no name is provided and start_at is specified, returns a
        list starting at that cursor position.
        """
        if name:
            name = quote_plus(name)
            return _get_by_name(entity, name)
        if start_at is not None and type(start_at) is int:
            return _get_with_cursor(entity, start_at)
        url = "%s/%s/%s" % (API_URL, VERSION, entity)
        result = json_from_url(url)
        return result

    @classmethod
    def next(cls, response):
        """Returns the next set of entities based on a previous API response."""
        start_at = _get_next(response)
        entity = _get_entity(response)
        return cls.get(entity=entity, start_at=start_at)

    @classmethod
    def download(cls, datasets=None, what=None):
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
            resp = input_shim(warning)
            if resp.lower() != "y":
                return

        for dataset in datasets:
            if dataset not in cls.DATASETS:
                msg = (
                    '"%s" is not a valid dataset name. Check the `DATASETS`'
                    " property for a complete list of names." % dataset
                )
                raise AttributeError(msg)
            if not os.path.exists(dataset):
                os.mkdir(dataset)

            if what is None:
                what = DOWNLOADS

            for dl in what:
                path = DATASET_TO_PATH[dataset]
                url = "%s/%s/%s" % (DOWNLOAD_URL, path, dl)

                try:
                    response = download_from_url(url)
                except HTTPError as e:
                    # Not every dataset has all downloads.
                    logging.warning("Skipping %s: HTTP Error %s" % (dl, e.code))
                    continue
                except Exception as e:
                    logging.warning("Skipping %s: %s" % (dl, e))
                    continue

                filename = "%s/%s" % (dataset, dl)
                filename = filename.replace(".gz", "")

                if response.code != 200:
                    logging.warning("Skipping %s: HTTP status %s" % (dl, response.code))
                    continue

                if os.path.isfile(filename):
                    logging.info("Using cached `%s`" % (filename))
                else:
                    logging.info("Downloading `%s`" % (filename))
                    _download_and_decompress_file(response, filename)

                yield filename

    @classmethod
    def download_df(cls, datasets=None, what=None, sparse=False, **kwargs):
        for file in cls.download(datasets, what):
            if sparse:
                yield _read_as_sparse_dataframe(file, **kwargs)
            else:
                yield _read_as_dataframe(file, **kwargs)

    @classmethod
    def get_gene_functional_annotations(
        cls, gene_symbol: str, datasets: list = None
    ) -> dict:
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
        except Exception as e:
            logging.warning(f"Could not get gene info for {gene_symbol}: {e}")
            gene_info = {
                "symbol": gene_symbol,
                "name": "Unknown",
                "description": "No description available",
            }

        # Process associations from API response
        associations = gene_info.get("associations", [])

        # Group associations by dataset
        dataset_associations = {}

        for assoc in associations:
            gene_set = assoc.get("geneSet", {})
            gene_set_name, dataset_name = GeneData._parse_gene_set(assoc)
            dataset_name = dataset_name or "Unknown Dataset"

            # Filter by specific datasets if requested
            if datasets and dataset_name not in datasets:
                continue

            if dataset_name not in dataset_associations:
                dataset_associations[dataset_name] = {
                    "increased": [],
                    "decreased": [],
                    "summary": f"Associations for {gene_symbol} in {dataset_name}",
                }

            # Extract score and direction
            threshold_value = assoc.get("thresholdValue", 0)
            standardized_value = assoc.get("standardizedValue", 0)

            # Use standardized value if available, otherwise threshold
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

        return {
            "gene_info": {
                "symbol": gene_info.get("symbol", gene_symbol),
                "name": gene_info.get("name", "Unknown"),
                "description": gene_info.get("description", "No description available"),
                "ncbi_id": gene_info.get("ncbiEntrezGeneId", "Unknown"),
            },
            "functional_associations": cls._format_functional_associations_from_grouped_data(
                dataset_associations
            ),
        }

    @classmethod
    def get_gene_with_associations(cls, gene_symbol: str) -> dict:
        """Get gene information with associations using the API.

        This uses the showAssociations=true parameter to get functional
        associations directly from the API without downloading files.

        Args:
            gene_symbol: Gene symbol

        Returns:
            Dictionary with gene info and associations
        """
        name = quote_plus(gene_symbol)
        url = "%s/%s/gene/%s?showAssociations=true" % (API_URL, VERSION, name)
        return json_from_url(url)

    @classmethod
    def _get_gene_dataset_annotations_from_api(
        cls, gene_symbol: str, dataset: str
    ) -> dict:
        """Get functional annotations for a gene in a specific dataset using API.

        Note: This method currently has limited functionality due to API constraints.
        For complete functional associations, use the download method.

        Args:
            gene_symbol: Gene symbol
            dataset: Dataset name

        Returns:
            Dictionary with annotation data for the dataset (may be empty)
        """
        try:
            # Get dataset info
            dataset_info = cls.get("dataset", dataset)
            summary = dataset_info.get(
                "description", f"Associations for {gene_symbol} in {dataset}"
            )

            # The API doesn't currently provide direct access to gene-gene set associations
            # This would require downloading the gene-attribute matrices
            logging.debug(
                f"API access to functional associations for {gene_symbol} in {dataset} is not available"
            )

            return None

        except Exception as e:
            logging.debug(
                f"Error processing dataset '{dataset}' for gene '{gene_symbol}': {e}"
            )
            return None

    @classmethod
    def download_gene_functional_annotations(
        cls, gene_symbol: str, datasets: list = None
    ) -> dict:
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
            try:
                dataset_annotations = cls._get_gene_dataset_annotations_from_files(
                    gene_symbol, dataset
                )
                if dataset_annotations:
                    results[dataset] = dataset_annotations
            except Exception as e:
                logging.debug(
                    f"Error processing dataset '{dataset}' for gene '{gene_symbol}': {e}"
                )
                continue

        return results

    @classmethod
    def _get_gene_dataset_annotations_from_files(
        cls, gene_symbol: str, dataset: str
    ) -> dict:
        """Get functional annotations for a gene in a specific dataset using downloaded files.

        Args:
            gene_symbol: Gene symbol
            dataset: Dataset name

        Returns:
            Dictionary with annotation data for the dataset
        """
        try:
            # Get dataset info
            dataset_info = cls.get("dataset", dataset)
            summary = dataset_info.get(
                "description", f"Associations for {gene_symbol} in {dataset}"
            )

            # Download the necessary files for this dataset
            gene_matrix_file = None
            attribute_list_file = None

            for filename in cls.download(
                [dataset],
                ["gene_attribute_matrix.txt.gz", "attribute_list_entries.txt.gz"],
            ):
                if "gene_attribute_matrix.txt" in filename:
                    gene_matrix_file = filename
                elif "attribute_list_entries.txt" in filename:
                    attribute_list_file = filename

            if not gene_matrix_file:
                logging.debug(
                    f"Could not download gene-attribute matrix for dataset '{dataset}'"
                )
                return None

            # Read the gene-attribute matrix
            import pandas as pd

            # Read the matrix file
            matrix_df = _read_as_dataframe(gene_matrix_file)

            # Check if the gene exists in the matrix
            # Genes are stored as JSON strings like ["GENE_NAME", "na", "ID"]
            gene_found = False
            gene_row = None

            for idx in matrix_df.index:
                try:
                    import json

                    gene_data = json.loads(idx)
                    if gene_data[0] == gene_symbol:
                        gene_found = True
                        gene_row = matrix_df.loc[idx]
                        break
                except (json.JSONDecodeError, IndexError):
                    continue

            if not gene_found:
                logging.debug(f"Gene '{gene_symbol}' not found in dataset '{dataset}'")
                return None

            # Read attribute list if available
            attribute_names = {}
            if attribute_list_file:
                try:
                    attr_df = pd.read_csv(
                        attribute_list_file, sep="\t", encoding="latin-1"
                    )
                    if len(attr_df.columns) >= 2:
                        # Attributes are also stored as JSON strings
                        for _, row in attr_df.iterrows():
                            try:
                                attr_data = json.loads(row.iloc[0])
                                attr_id = row.iloc[0]  # Use the full JSON string as ID
                                attr_name = (
                                    attr_data[0]
                                    if len(attr_data) > 0
                                    else str(row.iloc[0])
                                )
                                attribute_names[attr_id] = attr_name
                            except (json.JSONDecodeError, IndexError):
                                # Fallback to plain text
                                attribute_names[str(row.iloc[0])] = str(row.iloc[1])
                except Exception as e:
                    logging.debug(
                        f"Could not read attribute list for dataset '{dataset}': {e}"
                    )

            # Extract associations
            associations = []
            increased_associations = []
            decreased_associations = []

            for attr_id, score in gene_row.items():
                if pd.notna(score) and score != 0:  # Skip missing or zero values
                    try:
                        score_float = float(score)
                        attr_name = attribute_names.get(attr_id, attr_id)

                        association = {
                            "name": attr_name,
                            "score": score_float,
                            "attribute_id": attr_id,
                        }

                        associations.append(association)

                        if score_float > 0:
                            increased_associations.append(association)
                        elif score_float < 0:
                            decreased_associations.append(association)

                    except (ValueError, TypeError):
                        continue

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

        except Exception as e:
            logging.debug(
                f"Error processing dataset '{dataset}' for gene '{gene_symbol}': {e}"
            )
            return None

    @classmethod
    def _extract_gene_set_associations(
        cls, gene_set_detail: dict, gene_symbol: str
    ) -> list:
        """Extract associations from gene set details.

        Args:
            gene_set_detail: Detailed gene set information from API
            gene_symbol: Gene symbol to find associations for

        Returns:
            List of association dictionaries with name and score
        """
        associations = []

        # Look for attributes in the gene set
        attributes = gene_set_detail.get("attributes", [])

        for attr in attributes:
            attr_name = attr.get("name", "Unknown")
            attr_id = attr.get("id")

            if attr_id:
                try:
                    # Get attribute details to find the association score
                    attr_detail = cls.get("attribute", attr_id)

                    # Look for the gene in the attribute's gene associations
                    genes = attr_detail.get("genes", [])

                    for gene in genes:
                        if gene.get("symbol") == gene_symbol:
                            # Extract score from the association
                            score = cls._extract_association_score(gene, attr_detail)
                            if score is not None:
                                associations.append(
                                    {
                                        "name": attr_name,
                                        "score": score,
                                        "attribute_id": attr_id,
                                        "gene_set": gene_set_detail.get(
                                            "name", "Unknown"
                                        ),
                                    }
                                )
                            break

                except Exception as e:
                    logging.debug(
                        f"Error getting attribute details for {attr_name}: {e}"
                    )
                    continue

        return associations

    @classmethod
    def _extract_association_score(cls, gene_assoc: dict, attr_detail: dict) -> float:
        """Extract the association score from gene-attribute association.

        Args:
            gene_assoc: Gene association object
            attr_detail: Attribute detail object

        Returns:
            Association score as float, or None if not found
        """
        # Try different possible locations for the score
        score = gene_assoc.get("score")
        if score is not None:
            try:
                return float(score)
            except (ValueError, TypeError):
                pass

        # Look in the association object itself
        score = gene_assoc.get("association", {}).get("score")
        if score is not None:
            try:
                return float(score)
            except (ValueError, TypeError):
                pass

        # Look in attribute details for gene-specific scores
        genes = attr_detail.get("genes", [])
        for g in genes:
            if g.get("symbol") == gene_assoc.get("symbol"):
                score = g.get("score")
                if score is not None:
                    try:
                        return float(score)
                    except (ValueError, TypeError):
                        pass

        return None

    @classmethod
    def get_gene_associations_summary(
        cls, gene_symbol: str, datasets: list = None, use_download: bool = False
    ) -> dict:
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
        cls, gene_symbol: str, datasets: list = None, use_download: bool = False
    ) -> dict:
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
        # Get gene info
        try:
            gene_info = cls.get("gene", gene_symbol)
        except Exception as e:
            logging.warning(f"Could not get gene info for {gene_symbol}: {e}")
            gene_info = {
                "symbol": gene_symbol,
                "name": "Unknown",
                "description": "No description available",
            }

        # Get functional annotations
        if use_download:
            annotations = cls.download_gene_functional_annotations(
                gene_symbol, datasets
            )
        else:
            return cls.get_gene_functional_annotations(gene_symbol, datasets)

        return {
            "gene_info": {
                "symbol": gene_info.get("symbol", gene_symbol),
                "name": gene_info.get("name", "Unknown"),
                "description": gene_info.get("description", "No description available"),
                "ncbi_id": gene_info.get("ncbiEntrezGeneId", "Unknown"),
            },
            "functional_associations": cls._format_functional_associations_from_grouped_data(
                annotations
            ),
        }

    @classmethod
    def get_gene_data(cls, gene_symbol: str, use_cache: bool = False):
        if use_cache:
            return GeneData(cls._get_gene_with_associations_cached(gene_symbol))
        else:
            return GeneData(cls.get_gene_with_associations(gene_symbol))

    @staticmethod
    @cache_to_file
    def _get_gene_with_associations_cached(gene_symbol: str):
        return Harmonizome.get_gene_with_associations(gene_symbol)


# Utility functions
# -------------------------------------------------------------------------


def _get_with_cursor(entity, start_at):
    """Returns a list of entities based on cursor position."""
    url = "%s/%s/%s?cursor=%s" % (API_URL, VERSION, entity, str(start_at))
    result = json_from_url(url)
    return result


def _get_by_name(entity, name):
    """Returns a single entity based on name."""
    url = "%s/%s/%s/%s" % (API_URL, VERSION, entity, name)
    return json_from_url(url)


def _get_entity(response):
    """Returns the entity from an API response."""
    path = response["next"].split("?")[0]
    return path.split("/")[3]


def _get_next(response):
    """Returns the next property from an API response."""
    if response["next"]:
        return int(response["next"].split("=")[1])
    return None


# This function was adopted from here: http://stackoverflow.com/a/15353312.
def _download_and_decompress_file(response, filename):
    """Downloads and decompresses a single file from a response object."""
    compressed_file = BytesIO(response.read())
    decompressed_file = gzip.GzipFile(fileobj=compressed_file)

    with open(filename, "wb+") as outfile:
        outfile.write(decompressed_file.read())


def _getfshape(fn, row_sep="\n", col_sep="\t", open_args={}):
    """Fast and efficient way of finding row/col height of file"""
    with open(fn, "r", newline=row_sep, **open_args) as f:
        col_size = f.readline().count(col_sep) + 1
        row_size = sum(1 for line in f) + 1
        return (row_size, col_size)


def _parse(
    fn,
    column_size=3,
    index_size=3,
    shape=None,
    index_fmt=None,
    data_fmt=None,
    index_dtype=None,
    data_dtype=None,
    col_sep="\t",
    row_sep="\n",
    open_args={},
):
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
    import numpy as np

    if index_fmt is None:
        index_fmt = np.ndarray
    if data_fmt is None:
        data_fmt = np.ndarray
    if index_dtype is None:
        index_dtype = np.object
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

    with open(fn, "r", newline=row_sep, **open_args) as fh:
        header = np.array([next(fh).strip().split(col_sep) for _ in range(column_size)])

        column_names = header[:column_size, index_size - 1]
        index_names = header[column_size - 1, :index_size]

        columns[:, :] = header[:column_size, index_size:]

        for ind, line in enumerate(fh):
            lh = line.strip().split(col_sep)
            index[ind, :] = lh[:index_size]
            data[ind, :] = lh[index_size:]

        return (column_names, columns, index_names, index, data)


def _parse_df(
    fn,
    sparse=False,
    default_fill_value=None,
    column_apply=None,
    index_apply=None,
    df_args={},
    **kwargs,
):
    import numpy as np
    import pandas as pd
    from scipy.sparse import lil_matrix

    data_fmt = lil_matrix if sparse else np.ndarray
    df_type = pd.SparseDataFrame if sparse else pd.DataFrame
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

    return df_type(
        data=data.tocsr() if sparse else data,
        index=pd.Index(
            data=index,
            name=str(index_names),
            dtype=np.object,
        ),
        columns=pd.Index(
            data=columns,
            name=str(column_names),
            dtype=np.object,
        ),
        **df_args,
    )


def _df_column_uniquify(df):
    df_columns = df.columns
    new_columns = []
    for item in df_columns:
        counter = 0
        newitem = item
        while newitem in new_columns:
            counter += 1
            newitem = "{}_{}".format(item, counter)
        new_columns.append(newitem)
    df.columns = new_columns
    return df


def _json_ind_no_slash(ind_names, ind):
    return (
        json.dumps([ind_name.replace("/", "|") for ind_name in ind_names]),
        [json.dumps([ii.replace("/", "|") for ii in i]) for i in ind],
    )


def _read_as_dataframe(fn):
    """Standard loading of dataframe"""
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
        import pandas as pd

        return pd.read_table(fn, encoding="latin-1", index_col=None)
    else:
        raise Exception("Unable to parse this file into a dataframe.")


def _read_as_sparse_dataframe(fn, blocksize=10e6, fill_value=0):
    """Efficient loading sparse dataframe"""
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
        raise Exception("Unable to parse this file into a dataframe.")
