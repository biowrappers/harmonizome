"""Tests for the Harmonizome class (simplified for core functionality)."""

import gzip
import json
import os
import shutil
from io import BytesIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import lil_matrix

from harmonizome import Entity, Harmonizome
import harmonizome.harmonizome as harmonizome_module
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
        assert Harmonizome.__version__ == "1.0.1"

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

    @patch("harmonizome.harmonizome.json_from_url")
    def test_lazy_config_loads_on_demand(self, mock_json, monkeypatch):
        mock_json.return_value = {
            "downloads": ["gene_attribute_matrix.txt.gz"],
            "datasets": {"ENCODE": "encode/path"},
        }
        monkeypatch.setattr(harmonizome_module, "_CONFIG", None)
        monkeypatch.setattr(harmonizome_module, "DOWNLOADS", [])
        monkeypatch.setattr(harmonizome_module, "DATASET_TO_PATH", {})

        downloads = harmonizome_module._get_downloads()
        dataset_to_path = harmonizome_module._get_dataset_to_path()

        assert downloads == ["gene_attribute_matrix.txt.gz"]
        assert dataset_to_path == {"ENCODE": "encode/path"}
        assert "ENCODE" in Harmonizome.DATASETS
        mock_json.assert_called_once()

    @patch("harmonizome.harmonizome.json_from_url")
    def test_load_config_falls_back_on_http_error(self, mock_json):
        mock_json.side_effect = HTTPError(
            url="https://example.org/config",
            code=500,
            msg="broken",
            hdrs=None,
            fp=None,
        )

        config = harmonizome_module._load_config()

        assert config["datasets"] == {"ENCODE": "encode", "GTEx": "gtex"}
        assert "gene_attribute_matrix.txt.gz" in config["downloads"]

    @patch("builtins.input")
    def test_download_one_dataset(self, mock_input):
        mock_input.return_value = "y"
        test_dir = "ENCODE"
        mock_response = MagicMock()
        mock_response.code = 200
        try:
            with patch.object(
                Harmonizome, "DATASETS", new={"ENCODE": "encode/path"}
            ), patch(
                "harmonizome.harmonizome._get_dataset_to_path",
                return_value={"ENCODE": "encode/path"},
            ), patch(
                "harmonizome.harmonizome._get_downloads",
                return_value=["gene_attribute_matrix.txt.gz"],
            ), patch(
                "harmonizome.harmonizome.Path.mkdir"
            ), patch(
                "pathlib.Path.is_file", return_value=True
            ), patch(
                "harmonizome.harmonizome.download_from_url", return_value=mock_response
            ):
                filenames = list(Harmonizome.download(["ENCODE"]))
                assert isinstance(filenames, list)
                assert filenames == [str(Path(test_dir) / "gene_attribute_matrix.txt")]
        finally:
            # Clean up the ENCODE directory if it was created
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)

    def test_download_skips_http_errors(self):
        http_error = HTTPError(
            url="https://example.org/file.gz",
            code=404,
            msg="missing",
            hdrs=None,
            fp=None,
        )

        with patch.object(Harmonizome, "DATASETS", new={"ENCODE"}), patch(
            "harmonizome.harmonizome._get_dataset_to_path",
            return_value={"ENCODE": "encode/path"},
        ), patch(
            "harmonizome.harmonizome._get_downloads",
            return_value=["gene_attribute_matrix.txt.gz"],
        ), patch(
            "harmonizome.harmonizome.Path.mkdir"
        ), patch(
            "harmonizome.harmonizome.download_from_url",
            side_effect=http_error,
        ):
            filenames = list(Harmonizome.download(["ENCODE"]))

        assert filenames == []

    def test_download_skips_url_errors(self):
        with patch.object(Harmonizome, "DATASETS", new={"ENCODE"}), patch(
            "harmonizome.harmonizome._get_dataset_to_path",
            return_value={"ENCODE": "encode/path"},
        ), patch(
            "harmonizome.harmonizome._get_downloads",
            return_value=["gene_attribute_matrix.txt.gz"],
        ), patch(
            "harmonizome.harmonizome.Path.mkdir"
        ), patch(
            "harmonizome.harmonizome.download_from_url",
            side_effect=URLError("offline"),
        ):
            filenames = list(Harmonizome.download(["ENCODE"]))

        assert filenames == []

    def test_download_raises_for_invalid_dataset(self):
        with patch.object(Harmonizome, "DATASETS", new={"ENCODE"}):
            with pytest.raises(AttributeError, match="not a valid dataset name"):
                list(Harmonizome.download(["INVALID"]))

    @patch("builtins.input", return_value="n")
    def test_download_returns_without_confirmation(self, mock_input):
        with patch.object(Harmonizome, "DATASETS", new={"ENCODE"}):
            assert list(Harmonizome.download()) == []

    def test_download_skips_non_200_response(self):
        mock_response = MagicMock()
        mock_response.code = 500
        with patch.object(Harmonizome, "DATASETS", new={"ENCODE"}), patch(
            "harmonizome.harmonizome._get_dataset_to_path",
            return_value={"ENCODE": "encode/path"},
        ), patch(
            "harmonizome.harmonizome._get_downloads",
            return_value=["gene_attribute_matrix.txt.gz"],
        ), patch(
            "harmonizome.harmonizome.Path.mkdir"
        ), patch(
            "pathlib.Path.is_file", return_value=False
        ), patch(
            "harmonizome.harmonizome.download_from_url",
            return_value=mock_response,
        ):
            filenames = list(Harmonizome.download(["ENCODE"]))

        assert filenames == []

    @patch("harmonizome.harmonizome.Harmonizome.get_gene_with_associations")
    def test_get_gene_functional_annotations_falls_back_on_http_error(
        self, mock_gene_with_associations
    ):
        mock_gene_with_associations.side_effect = HTTPError(
            url="https://example.org/gene",
            code=500,
            msg="broken",
            hdrs=None,
            fp=None,
        )

        result = Harmonizome.get_gene_functional_annotations("STAT3")

        assert result["gene_info"] == {
            "symbol": "STAT3",
            "name": "Unknown",
            "description": "No description available",
            "ncbi_id": "Unknown",
        }
        assert result["functional_associations"]["datasets"] == []

    def test_download_df_uses_requested_loader(self):
        with patch.object(
            Harmonizome,
            "download",
            return_value=["file1.txt"],
        ), patch(
            "harmonizome.harmonizome._read_as_dataframe",
            return_value="dense",
        ) as mock_dense, patch(
            "harmonizome.harmonizome._read_as_sparse_dataframe",
            return_value="sparse",
        ) as mock_sparse:
            dense_result = list(Harmonizome.download_df(["ENCODE"]))
            sparse_result = list(Harmonizome.download_df(["ENCODE"], sparse=True))

        assert dense_result == ["dense"]
        assert sparse_result == ["sparse"]
        mock_dense.assert_called_once_with("file1.txt")
        mock_sparse.assert_called_once_with("file1.txt")

    @patch("harmonizome.harmonizome.Harmonizome.download_gene_functional_annotations")
    def test_get_gene_associations_summary_for_downloads(self, mock_download_annotations):
        mock_download_annotations.return_value = {
            "Dataset1": {
                "total_associations": 2,
                "increased_count": 1,
                "decreased_count": 1,
            }
        }

        result = Harmonizome.get_gene_associations_summary("STAT3", use_download=True)

        assert result["total_datasets"] == 1
        assert result["total_associations"] == 2
        assert result["total_increased_associations"] == 1
        assert result["total_decreased_associations"] == 1

    @patch("harmonizome.harmonizome.Harmonizome._get_gene_with_associations_cached")
    def test_get_gene_data_uses_cache_when_requested(self, mock_cached):
        mock_cached.return_value = {"associations": []}

        result = Harmonizome.get_gene_data("STAT3", use_cache=True)

        assert isinstance(result, GeneData)
        mock_cached.assert_called_once_with("STAT3")

    @patch("harmonizome.harmonizome.Harmonizome.get_gene_with_associations")
    def test_get_gene_data_uses_live_fetch_by_default(self, mock_live):
        mock_live.return_value = {"associations": []}

        result = Harmonizome.get_gene_data("STAT3")

        assert isinstance(result, GeneData)
        mock_live.assert_called_once_with("STAT3")


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


class TestAnnotationHelpers:
    def test_get_gene_row_from_matrix_finds_requested_gene(self):
        matrix_df = pd.DataFrame(
            {"attr1": [1.0, 0.0]},
            index=['["STAT3", "na", "1"]', '["TP53", "na", "2"]'],
        )

        gene_row = Harmonizome._get_gene_row_from_matrix(matrix_df, "STAT3", "Test")

        assert gene_row is not None
        assert gene_row["attr1"] == 1.0

    def test_load_attribute_names_returns_empty_for_missing_file(self):
        assert Harmonizome._load_attribute_names(None) == {}
        assert Harmonizome._load_attribute_names("does-not-exist.tsv") == {}

    def test_build_annotation_summary_groups_scores(self):
        gene_row = pd.Series(
            {
                "attr_positive": 2.5,
                "attr_negative": -1.0,
                "attr_zero": 0,
                "attr_bad": "not-a-number",
            }
        )
        attribute_names = {
            "attr_positive": "Positive association",
            "attr_negative": "Negative association",
        }

        summary = Harmonizome._build_annotation_summary(
            "Associations for STAT3", gene_row, attribute_names
        )

        assert summary is not None
        assert summary["total_associations"] == 2
        assert summary["increased_count"] == 1
        assert summary["decreased_count"] == 1
        assert summary["increased_associations"][0]["name"] == "Positive association"
        assert summary["decreased_associations"][0]["name"] == "Negative association"

    def test_read_as_dataframe_raises_for_unsupported_file(self):
        with pytest.raises(ValueError, match="Unable to parse this file"):
            harmonizome_module._read_as_dataframe("unsupported_file.txt")


class TestUtilityFunctions:
    def test_json_from_url_falls_back_to_latin_1(self):
        mock_response = MagicMock()
        mock_response.read.return_value = b'{"name":"caf\xe9"}'

        with patch("harmonizome.harmonizome.urlopen", return_value=mock_response):
            result = harmonizome_module.json_from_url("https://example.org")

        assert result == {"name": "café"}

    def test_json_from_url_logs_and_raises_decode_error(self):
        mock_response = MagicMock()
        mock_response.read.return_value = b"not json"

        with patch("harmonizome.harmonizome.urlopen", return_value=mock_response):
            with pytest.raises(json.JSONDecodeError):
                harmonizome_module.json_from_url("https://example.org")

    def test_download_from_url_raises_url_error(self):
        with patch(
            "harmonizome.harmonizome.urlopen", side_effect=URLError("offline")
        ):
            with pytest.raises(URLError):
                harmonizome_module.download_from_url("https://example.org")

    @patch("harmonizome.harmonizome.json_from_url")
    def test_get_with_cursor_uses_expected_url(self, mock_json):
        mock_json.return_value = {"cursor": 2}

        result = harmonizome_module._get_with_cursor("gene", 2)

        assert result == {"cursor": 2}
        mock_json.assert_called_once_with(
            f"{harmonizome_module.API_URL}/{harmonizome_module.VERSION}/gene?cursor=2"
        )

    def test_get_entity_and_next_helpers(self):
        response = {"next": "https://host/api/v1/gene?cursor=5"}

        assert harmonizome_module._get_entity(response) == "gene"
        assert harmonizome_module._get_next(response) == 5
        assert harmonizome_module._get_next({"next": None}) is None

    def test_download_and_decompress_file(self, tmp_path):
        expected = b"gene\tvalue\nSTAT3\t1\n"
        compressed_buffer = BytesIO()
        with gzip.GzipFile(fileobj=compressed_buffer, mode="wb") as gz_file:
            gz_file.write(expected)

        response = MagicMock()
        response.read.return_value = compressed_buffer.getvalue()
        output_file = tmp_path / "output.txt"

        harmonizome_module._download_and_decompress_file(response, output_file)

        assert output_file.read_bytes() == expected

    def test_getfshape_counts_rows_and_columns(self, tmp_path):
        matrix_file = tmp_path / "matrix.txt"
        matrix_file.write_text("a\tb\tc\n1\t2\t3\n4\t5\t6\n", encoding="utf-8")

        assert harmonizome_module._getfshape(str(matrix_file)) == (3, 3)

    def test_df_column_uniquify_appends_suffixes(self):
        dataframe = pd.DataFrame([[1, 2]], columns=["dup", "dup"])

        result = harmonizome_module._df_column_uniquify(dataframe)

        assert list(result.columns) == ["dup", "dup_1"]

    def test_json_ind_no_slash_rewrites_names(self):
        column_name, values = harmonizome_module._json_ind_no_slash(
            ["a/b"], [["x/y", "z"]]
        )

        assert column_name == '["a|b"]'
        assert values == ['["x|y", "z"]']

    def test_parse_df_raises_for_unexpected_sparse_kwargs(self):
        with patch("harmonizome.harmonizome._parse") as mock_parse:
            sparse_matrix = lil_matrix([[1.0]])
            mock_parse.return_value = (
                np.array(["column"]),
                np.array([["column"]], dtype=object),
                np.array(["index"]),
                np.array([["row"]], dtype=object),
                sparse_matrix,
            )

            with pytest.raises(TypeError, match="Unexpected DataFrame arguments"):
                harmonizome_module._parse_df(
                    "unused.txt",
                    sparse=True,
                    df_args={"unexpected": "value"},
                )

    def test_read_as_sparse_dataframe_raises_for_unsupported_file(self):
        with pytest.raises(ValueError, match="Unable to parse this file"):
            harmonizome_module._read_as_sparse_dataframe("unsupported_file.txt")


class TestFunctionalAssociationFormattingAdditional:
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
