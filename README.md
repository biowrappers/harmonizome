# Harmonizome Python Client

A Python client for the [Harmonizome API](https://maayanlab.cloud/Harmonizome/), a resource for exploring gene sets and their attributes across multiple datasets.

## Installation

```bash
pip install harmonizome
```

## Quick Start

```python
from harmonizome import Harmonizome

# Get all available datasets
datasets = Harmonizome.DATASETS
print(f"Available datasets: {len(datasets)}")

# Get information about a specific gene
gene_info = Harmonizome.get('gene', 'BRCA1')
print(f"Gene: {gene_info['name']}")

# Download data for a specific dataset
for filename in Harmonizome.download(['ENCODE']):
    print(f"Downloaded: {filename}")

# Load data as pandas DataFrames
for df in Harmonizome.download_df(['ENCODE'], sparse=False):
    print(f"DataFrame shape: {df.shape}")
```

## Features

- **API Access**: Query genes, gene sets, attributes, and datasets
- **Data Download**: Download complete datasets (~30GB total)
- **DataFrame Support**: Load data directly as pandas DataFrames
- **Sparse Matrix Support**: Efficient handling of large sparse datasets
- **Python 3.9+ Support**: Targets the package minimum version and newer

## API Reference

### Core Methods

#### `Harmonizome.get(entity, name=None, start_at=None)`
Retrieve entities from the Harmonizome API.

- `entity`: Type of entity ('gene', 'gene_set', 'attribute', etc.)
- `name`: Specific entity name (optional)
- `start_at`: Cursor position for pagination (optional)

#### `Harmonizome.download(datasets=None, what=None)`
Download dataset files to local directories.

- `datasets`: List of dataset names (defaults to all datasets)
- `what`: List of file types to download (defaults to all types)

#### `Harmonizome.download_df(datasets=None, what=None, sparse=False)`
Download and load data as pandas DataFrames.

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
# Get all genes (paginated)
genes = Harmonizome.get('gene')
print(f"Found {len(genes['entities'])} genes")

# Get next page
next_genes = Harmonizome.next(genes)

# Get specific gene
brca1 = Harmonizome.get('gene', 'BRCA1')
print(f"BRCA1 description: {brca1['description']}")
```

### Downloading Datasets

```python
# Download all data (requires confirmation)
for filename in Harmonizome.download():
    print(f"Downloaded: {filename}")

# Download specific datasets
datasets = ['ENCODE', 'GTEx']
for filename in Harmonizome.download(datasets):
    print(f"Downloaded: {filename}")

# Download specific file types
file_types = ['gene_attribute_matrix.txt.gz']
for filename in Harmonizome.download(['ENCODE'], file_types):
    print(f"Downloaded: {filename}")
```

### Working with DataFrames

```python
# Load as regular DataFrames
for df in Harmonizome.download_df(['ENCODE']):
    print(f"DataFrame: {df.shape}")
    print(f"Columns: {df.columns[:5].tolist()}")
    break

# Load as sparse DataFrames (memory efficient)
for df in Harmonizome.download_df(['ENCODE'], sparse=True):
    print(f"Sparse DataFrame: {df.shape}")
    print(f"Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
    break
```

### Working with Gene Associations as DataFrames

You can fetch all associations for a gene and convert them directly to a pandas DataFrame for easy filtering and export:

```python
from harmonizome import Harmonizome

gene = "STAT3"
gene_data = Harmonizome.get_gene_data(gene, use_cache=True)

# Get all associations as a DataFrame
df = gene_data.to_dataframe()

# List all unique dataset names
print("Available datasets for this gene:")
for i, dataset in enumerate(df["dataset"].unique(), 1):
    print(f"{i}. {dataset}")

# Select a dataset by name
selected_dataset = df["dataset"].unique()[0]  # or set to any dataset name from the list
print(f"\nSelected dataset: {selected_dataset}")

# Filter associations for the selected dataset
dataset_df = df[df["dataset"] == selected_dataset]
print(dataset_df)

# Save to CSV
safe_name = selected_dataset.replace(" ", "_").replace("/", "_")
dataset_df.to_csv(f"{gene.lower()}_{safe_name}_associations.csv", index=False)
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

## Citation

If you use this package in your research, please cite:

```txt
Rouillard AD, Gundersen GW, Fernandez NF, Wang Z, Monteiro CD, McDermott MG, Ma'ayan A. The harmonizome: a collection of processed datasets gathered to serve and mine knowledge about genes and proteins. Database (Oxford). 2016 Jul 3;2016:baw100. doi: 10.1093/database/baw100.
```

## Links

- [Harmonizome Website](https://maayanlab.cloud/Harmonizome/)
- [API Documentation](https://maayanlab.cloud/Harmonizome/api)
- [GitHub Repository](https://github.com/maayanlab/harmonizome) 
