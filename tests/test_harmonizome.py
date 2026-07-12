"""Tests for the Harmonizome class (simplified for core functionality)."""

import os
import shutil
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from harmonizome import Entity, Harmonizome
from harmonizome.harmonizome import GeneData


class TestEntity:
    def test_entity_values(self):
        assert hasattr(Entity, "DATASET")
        assert hasattr(Entity, "GENE")
        assert hasattr(Entity, "GENE_SET")
        assert hasattr(Entity, "ATTRIBUTE")

    def test_entity_getattr(self):
        assert Entity.DATASET == "dataset"
        assert Entity.GENE == "gene"
        assert Entity.GENE_SET == "gene_set"

    def test_entity_invalid_attr(self):
        with pytest.raises(AttributeError):
            _ = Entity.INVALID_ENTITY


class TestHarmonizome:
    def test_version(self):
        assert hasattr(Harmonizome, "__version__")
        assert Harmonizome.__version__ == "1.0"

    def test_datasets_property(self):
        assert hasattr(Harmonizome, "DATASETS")
        # Accept dict_keys or similar iterable
        assert hasattr(Harmonizome.DATASETS, "__iter__")

    @patch("harmonizome.harmonizome.json_from_url")
    def test_get_entity_by_name(self, mock_json):
        mock_response = {"name": "BRCA1", "description": "Test gene"}
        mock_json.return_value = mock_response
        result = Harmonizome.get("gene", "BRCA1")
        assert result == mock_response
        mock_json.assert_called_once()

    @patch("harmonizome.harmonizome.input_shim")
    def test_download_one_dataset(self, mock_input):
        mock_input.return_value = "y"
        test_dir = "ENCODE"
        try:
            with patch.object(
                Harmonizome, "DATASETS", new={"ENCODE": "encode/path"}
            ), patch(
                "harmonizome.harmonizome.DATASET_TO_PATH", {"ENCODE": "encode/path"}
            ), patch(
                "harmonizome.harmonizome.DOWNLOADS", ["gene_attribute_matrix.txt.gz"]
            ), patch(
                "harmonizome.harmonizome.os.path.exists", return_value=True
            ), patch(
                "harmonizome.harmonizome.os.mkdir"
            ) as mock_mkdir:
                filenames = list(Harmonizome.download(["ENCODE"]))
                assert isinstance(filenames, list)
        finally:
            # Clean up the ENCODE directory if it was created
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)


class TestGeneData:
    def test_to_dataframe_all_and_filter(self):
        gene_info = {
            "associations": [
                {
                    "geneSet": {
                        "name": "SetA/Dataset1",
                        "dataset": {"name": "Dataset1"},
                    },
                    "thresholdValue": 1.2,
                    "standardizedValue": 0.8,
                },
                {
                    "geneSet": {
                        "name": "SetB/Dataset2",
                        "dataset": {"name": "Dataset2"},
                    },
                    "thresholdValue": 2.1,
                    "standardizedValue": 1.5,
                },
            ]
        }
        gene_data = GeneData(gene_info)
        df = gene_data.to_dataframe()
        assert set(df.columns) == {
            "gene_set",
            "dataset",
            "thresholdValue",
            "standardizedValue",
        }
        assert len(df) == 2
        assert set(df["dataset"]) == {"Dataset1", "Dataset2"}
        df1 = gene_data.to_dataframe(dataset="Dataset1")
        assert len(df1) == 1
        assert df1.iloc[0]["dataset"] == "Dataset1"
        assert df1.iloc[0]["gene_set"] == "SetA"
        df_none = gene_data.to_dataframe(dataset="Nonexistent")
        assert df_none.empty

    def test_to_dataframe_and_filter_with_dataset_embedded_in_gene_set_name(self):
        gene_info = {
            "associations": [
                {
                    "geneSet": {"name": "SetA/Dataset1"},
                    "thresholdValue": 1.2,
                    "standardizedValue": 0.8,
                },
                {
                    "geneSet": {"name": "SetB/Dataset2"},
                    "thresholdValue": 2.1,
                    "standardizedValue": 1.5,
                },
            ]
        }
        gene_data = GeneData(gene_info)
        df = gene_data.to_dataframe(dataset="Dataset1")
        assert len(df) == 1
        assert df.iloc[0]["dataset"] == "Dataset1"
        assert df.iloc[0]["gene_set"] == "SetA"


class TestFunctionalAssociationFormatting:
    @patch("harmonizome.harmonizome.Harmonizome.get_gene_with_associations")
    def test_get_gene_functional_annotations_formats_groups_once(
        self, mock_gene_with_associations
    ):
        mock_gene_with_associations.return_value = {
            "symbol": "STAT3",
            "name": "signal transducer and activator of transcription 3",
            "description": "Test gene",
            "ncbiEntrezGeneId": 6774,
            "associations": [
                {
                    "geneSet": {
                        "name": "SetA/Dataset1",
                        "dataset": {"name": "Dataset1"},
                    },
                    "thresholdValue": 1.0,
                    "standardizedValue": 2.5,
                },
                {
                    "geneSet": {
                        "name": "SetB/Dataset1",
                        "dataset": {"name": "Dataset1"},
                    },
                    "thresholdValue": -1.0,
                    "standardizedValue": -3.0,
                },
            ],
        }

        result = Harmonizome.get_gene_functional_annotations("STAT3")

        assert result["gene_info"]["symbol"] == "STAT3"
        functional_associations = result["functional_associations"]
        assert functional_associations["total_datasets"] == 1
        assert functional_associations["total_associations"] == 2
        assert functional_associations["total_increased"] == 1
        assert functional_associations["total_decreased"] == 1
        assert functional_associations["datasets"][0]["dataset"] == "Dataset1"
        assert [
            group["type"]
            for group in functional_associations["datasets"][0]["associations"]
        ] == [
            "increased",
            "decreased",
        ]

    @patch("harmonizome.harmonizome.Harmonizome.get_gene_with_associations")
    def test_get_gene_functional_annotations_uses_gene_set_name_when_dataset_field_missing(
        self,
        mock_gene_with_associations,
    ):
        mock_gene_with_associations.return_value = {
            "symbol": "STAT3",
            "name": "signal transducer and activator of transcription 3",
            "description": "Test gene",
            "ncbiEntrezGeneId": 6774,
            "associations": [
                {
                    "geneSet": {"name": "SetA/Dataset1"},
                    "thresholdValue": 1.0,
                    "standardizedValue": 2.5,
                },
                {
                    "geneSet": {"name": "SetB/Dataset1"},
                    "thresholdValue": -1.0,
                    "standardizedValue": -3.0,
                },
            ],
        }

        result = Harmonizome.get_gene_functional_annotations("STAT3", ["Dataset1"])

        functional_associations = result["functional_associations"]
        assert functional_associations["total_datasets"] == 1
        assert functional_associations["total_associations"] == 2
        assert functional_associations["datasets"][0]["dataset"] == "Dataset1"
        assert [
            item["name"]
            for group in functional_associations["datasets"][0]["associations"]
            for item in group["items"]
        ] == ["SetA", "SetB"]

    @patch("harmonizome.harmonizome.Harmonizome.get")
    @patch("harmonizome.harmonizome.Harmonizome.download_gene_functional_annotations")
    def test_get_gene_functional_associations_formatted_uses_shared_formatter_for_downloads(
        self,
        mock_download_annotations,
        mock_get,
    ):
        mock_get.return_value = {
            "symbol": "STAT3",
            "name": "signal transducer and activator of transcription 3",
            "description": "Test gene",
            "ncbiEntrezGeneId": 6774,
        }
        mock_download_annotations.return_value = {
            "Dataset1": {
                "summary": "Associations for STAT3 in Dataset1",
                "increased_associations": [{"name": "SetA", "score": 2.5}],
                "decreased_associations": [{"name": "SetB", "score": -3.0}],
            }
        }

        result = Harmonizome.get_gene_functional_associations_formatted(
            "STAT3",
            use_download=True,
        )

        functional_associations = result["functional_associations"]
        assert functional_associations["total_datasets"] == 1
        assert functional_associations["total_associations"] == 2
        assert functional_associations["total_increased"] == 1
        assert functional_associations["total_decreased"] == 1
        assert functional_associations["datasets"][0]["dataset"] == "Dataset1"
        assert [
            group["type"]
            for group in functional_associations["datasets"][0]["associations"]
        ] == [
            "increased",
            "decreased",
        ]
