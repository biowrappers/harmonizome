import pytest

from harmonizome import Harmonizome


def test_file_based_gene_annotations():
    gene = "STAT3"
    dataset = "Achilles Cell Line Gene Essentiality Profiles"
    annotations = Harmonizome.get_gene_functional_annotations(gene, [dataset])
    assert isinstance(annotations, dict)
    assert "functional_associations" in annotations
    datasets = annotations["functional_associations"].get("datasets", [])
    if not datasets:
        pytest.skip(
            "No datasets returned by API for this gene/dataset (may be API/data issue)"
        )
