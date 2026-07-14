# Harmonizome Python Wrapper

A Python wrapper around the [Ma'ayan Lab Harmonizome API and download endpoints](https://maayanlab.cloud/Harmonizome/). The package provides a Python interface for querying Harmonizome entities, downloading Harmonizome datasets, and loading the downloaded artifacts into pandas DataFrames.

## Installation

```bash
pip install harmonizome
```

## Quick Start

```python
from harmonizome import Harmonizome

# Inspect one gene from the API
gene_info = Harmonizome.get('gene', 'BRCA1')
print(f"Gene: {gene_info['name']}")

# Download one dataset artifact
for filename in Harmonizome.download(['ENCODE']):
    print(f"Downloaded: {filename}")
```

## Features

- **API Wrapper**: Query genes, gene sets, attributes, and datasets from the Ma'ayan Lab Harmonizome API
- **Dataset Downloads**: Download complete Harmonizome datasets (~30GB total) from the Ma'ayan Lab service
- **DataFrame Support**: Load data directly as pandas DataFrames
- **Sparse Matrix Support**: Efficient handling of large sparse datasets
- **Python 3.9+ Support**: Targets the package minimum version and newer

## Quick API Reference

### Core Methods

#### `Harmonizome.get(entity, name=None, start_at=None)`

Retrieve entities from the Ma'ayan Lab Harmonizome API.

- `entity`: Type of entity ('gene', 'gene_set', 'attribute', etc.)
- `name`: Specific entity name (optional)
- `start_at`: Cursor position for pagination (optional)

#### `Harmonizome.download(datasets=None, what=None)`

Download Harmonizome dataset files from the Ma'ayan Lab service to local directories.

- `datasets`: List of dataset names (defaults to all datasets)
- `what`: List of file types to download (defaults to all types)

#### `Harmonizome.download_df(datasets=None, what=None, sparse=False)`

Download Harmonizome dataset files and load them as pandas DataFrames.

- `datasets`: List of dataset names
- `what`: List of file types to download
- `sparse`: Use sparse matrices for memory efficiency

### Entity Types

- `DATASET`: Dataset information
- `GENE`: Gene information
- `GENE_SET`: Gene set collections
- `ATTRIBUTE`: Gene attributes
- `GENE_FAMILY`: Gene family classifications
- `NAMING_AUTHORITY`: Naming authorities
- `PROTEIN`: Protein information
- `RESOURCE`: Data resources

## Examples

### Querying Genes

```python
brca1 = Harmonizome.get('gene', 'BRCA1')
print(f"BRCA1 description: {brca1['description']}")
```

### Downloading Datasets

```python
# Download one dataset
for filename in Harmonizome.download(['ENCODE']):
    print(f"Downloaded: {filename}")
```

### Working with DataFrames

```python
# Load one dataset into a pandas DataFrame
for df in Harmonizome.download_df(['ENCODE']):
    print(f"DataFrame: {df.shape}")
    break
```

### Working with Gene Associations as DataFrames

You can fetch all associations for a gene and convert them to a pandas DataFrame:

```python
from harmonizome import Harmonizome

gene = "STAT3"
gene_data = Harmonizome.get_gene_data(gene, use_cache=True)

# Get all associations as a DataFrame
df = gene_data.to_dataframe()

# Filter to one dataset
dataset_df = df[df["dataset"] == "ENCODE Transcription Factor Binding Site Profiles"]
print(dataset_df.head())
```

### API Reference: GeneData.to_dataframe()

```python
gene_data.to_dataframe(dataset: str = None) -> pandas.DataFrame
```

- Returns a DataFrame with columns: 'gene_set', 'dataset', 'thresholdValue', 'standardizedValue'.
- Optionally filter by dataset name.

## File Types

The following file types are available for download:

- `gene_attribute_matrix.txt.gz`: Gene-attribute association matrix
- `gene_list_terms.txt.gz`: List of genes with terms
- `attribute_list_entries.txt.gz`: List of attributes with entries

## Requirements

- Python 3.9+
- numpy >= 1.19.0
- pandas >= 1.3.0
- scipy >= 1.7.0

## Development

```bash
# Install development dependencies
pip install -e .

# Run tests
pytest
```

## License

This project is licensed under the [MIT License](LICENSE).

## Data Source and Citation

This wrapper depends on the public Harmonizome resource maintained by the Ma'ayan Lab. If you use this package in your research, cite the Harmonizome resource itself:

```txt
Rouillard AD, Gundersen GW, Fernandez NF, Wang Z, Monteiro CD, McDermott MG, Ma'ayan A. The harmonizome: a collection of processed datasets gathered to serve and mine knowledge about genes and proteins. Database (Oxford). 2016 Jul 3;2016:baw100. doi: 10.1093/database/baw100.
```

## Links

- [Harmonizome Website](https://maayanlab.cloud/Harmonizome/)
- [API Documentation](https://maayanlab.cloud/Harmonizome/api) 
