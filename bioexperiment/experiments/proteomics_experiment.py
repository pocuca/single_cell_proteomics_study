from functools import cached_property
from typing import Dict, List

import anndata
import numpy as np
import pandas as pd
from anndata import AnnData
from base_experiment import Experiment

import bioexperiment.data.loader as load
from bioexperiment.data.loader_definitions import DATA_FOLDER
from bioexperiment.data.preprocess import (
    ann_data_map_from_observation_names,
    apply_data_filters,
    apply_data_transformations,
    build_cell_cycle_markers,
    filter_out_category,
    get_data_intersection,
    score_genes_with_cell_markers,
)
from bioexperiment.visualizations.plots import (
    cell_cycle_stage_box_plot,
    cell_cycle_stage_prediction_plot,
    completeness_distribution_plot,
    pca_plot,
    pearson_corr_analysis_plot,
    variance_analysis_plot,
    variance_comparison_box_plot,
)


class ProteomicsExperiment(Experiment):

    """
    Contains necessary methods for loading data, transforming and
    creating plots for the analysis.
    """

    cell_cycle_stage_map = {
        "_G1_": "G1",
        "_TB_": "G1-S",
        "_G2_": "G2",
        "_NB_": "G2-M",
    }

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self._filtered_raw_protein_data = None
        self._filtered_protein_data = None
        self._cell_cycle_markers: Dict[str, List[str]] = {}

    # Cached properties represent the lazy evaluation
    # of the data structures used.
    @cached_property
    def cc_markers(self) -> pd.DataFrame:
        return pd.read_excel(DATA_FOLDER / "CellCycleMarker.xlsx", "Tami_Geiger")

    @cached_property
    def smartseq_data(self) -> AnnData:
        return load.read_smartseq_data(DATA_FOLDER)

    @cached_property
    def dropseq_data(self) -> AnnData:
        return load.read_dropseq_data(DATA_FOLDER)

    @cached_property
    def full_protein_data(self) -> Dict[str, AnnData]:
        return load.read_protein_data(DATA_FOLDER)

    @cached_property
    def core_proteome_data(self) -> List[str]:
        return load.read_proteome_file(DATA_FOLDER)

    @cached_property
    def full_data_dictionary(self) -> Dict[str, AnnData]:
        return {
            "Proteins": self.protein_data,
            "SMARTseq2": self.smartseq_data,
            "Drop-Seq": self.dropseq_data,
        }

    @property
    def protein_data(self) -> AnnData:
        return self.full_protein_data["norm"]

    @property
    def raw_protein_data(self) -> AnnData:
        return self.full_protein_data["raw"]

    @staticmethod
    def _copy_all_ann_data(ann_data_dict: Dict[str, AnnData]):
        return {key: value.copy() for key, value in ann_data_dict.items()}

    def preprocess_data(self):
        # 0. Mapping cell cycle stage
        ann_data_map_from_observation_names(
            self.raw_protein_data,
            "cell cycle stage",
            self.cell_cycle_stage_map,
        )
        ann_data_map_from_observation_names(
            self.protein_data,
            "cell cycle stage",
            self.cell_cycle_stage_map,
        )
        # 1. Filtering data
        self._filtered_raw_protein_data = self.raw_protein_data.copy()
        apply_data_filters(self._filtered_raw_protein_data)

        self._filtered_protein_data = self.protein_data.copy()
        apply_data_filters(self._filtered_protein_data)

        # 2. Preparing markers
        self._cell_cycle_markers = build_cell_cycle_markers(
            self._filtered_protein_data, self.cc_markers
        )

    def run_true_single_cell_analysis(self):
        # Fig 4C
        cell_cycle_stage_box_plot(
            self._filtered_raw_protein_data, str(self.results_full_path / "Fig 4C.png")
        )
        # Fig 4D
        pca_plot(
            filter_out_category(
                apply_data_transformations(
                    self._filtered_protein_data,
                    log_transform_data=True,
                ),
                "cell cycle stage",
                filter_out_categories=["other"],
            ),
            color_column="cell cycle stage",
            plot_path=str(self.results_full_path / "Fig 4D.png"),
        )
        # Fig 4E
        # This requires a special filter, with higher number of minimum cells
        data_for_prediction = self.protein_data.copy()
        apply_data_filters(data_for_prediction, min_part_cells=0.7, min_num_genes=600)
        scored_data_for_prediction = score_genes_with_cell_markers(
            filter_out_category(
                apply_data_transformations(
                    data_for_prediction,
                    log_transform_data=True,
                    fill_nans=True,
                ),
                "cell cycle stage",
                filter_out_categories=["other"],
            ),
            self._cell_cycle_markers,
        )
        cell_cycle_stage_prediction_plot(
            scored_data_for_prediction,
            ("G2-M", "G1-S"),
            plot_path=str(self.results_full_path / "Fig 4E.png"),
        )

    def _run_transcriptomics_completeness_analysis(self):
        """
        Produces completeness distribution plot, figure 5A.
        """
        full_data_copy = self._copy_all_ann_data(self.full_data_dictionary)
        for adata in full_data_copy.values():
            apply_data_filters(adata, min_part_cells=None)
        completeness_distribution_plot(
            full_data_copy, str(self.results_full_path / "Fig 5A.png")
        )

    def _run_transcriptomics_pca_analysis(self):
        """
        Produces PCA visualisation of transcriptomics data, figure 5B.
        """
        # Fig 5B
        full_data_copy = self._copy_all_ann_data(self.full_data_dictionary)

        for adata in full_data_copy.values():
            apply_data_filters(adata)

        intersection_data = get_data_intersection(full_data_copy)

        for adata in intersection_data.values():
            adata.X = np.nan_to_num(adata.X)
            apply_data_transformations(
                adata,
                log_transform_data=True,
                fill_nans=True,
                normalize=True,
                copy_data=False,
            )

        for datasetname, adata in intersection_data.items():
            adata.obs["dataset"] = datasetname

        joined = anndata.concat([adata for adata in intersection_data.values()])

        pca_plot(
            joined,
            color_column="dataset",
            plot_path=str(self.results_full_path / "Fig 5B.png"),
            impute_downshifted=False,
            invert_axis=True,
        )

    def _run_transcriptomics_variance_comparison_and_pearson(self):
        """
        Produces variance comparison plot and pearson correlation. Plots are done
        together in the same method because they share part of the data pipeline.
        """
        # Fig 5C & # Fig 5E

        # These two figures are put in the same pipeline because
        # both of these plots share most of the
        # transformations necessary for their computation
        full_data_copy = self._copy_all_ann_data(self.full_data_dictionary)

        for adata in full_data_copy.values():
            apply_data_filters(adata)

        for name, adata in full_data_copy.items():
            if name != "Proteins":
                apply_data_transformations(
                    adata,
                    normalize=True,
                    norm_target_sum=np.mean(np.nansum(adata.X, axis=1)),
                    copy_data=False,
                )
        variance_comparison_box_plot(
            self.core_proteome_data,
            full_data_copy,
            self.results_full_path / "Fig 5E.png",
        )

        for _, adata in full_data_copy.items():
            apply_data_transformations(
                adata.X, log_transform_data=True, copy_data=False
            )

        pearson_corr_analysis_plot(
            full_data_copy, plot_path=str(self.results_full_path / "Fig 5C.png")
        )

    def _run_variance_analysis_plot(self):
        """
        Produces variance analysis plot, figure 5D.
        """
        # Fig 5D
        variance_analysis_plot(
            self.core_proteome_data,
            self._filtered_protein_data,
            plot_path=str(self.results_full_path / "Fig 5D.png"),
        )

    def run_comparison_transcriptomics_data(self):
        self._run_transcriptomics_completeness_analysis()
        self._run_transcriptomics_pca_analysis()
        self._run_transcriptomics_variance_comparison_and_pearson()
        self._run_variance_analysis_plot()

    # To be done, in the similar way as run_comparison_transcriptomics_data
    def run_supplementary_analysis(self):
        pass

    def run(self):
        self.preprocess_data()
        self.run_true_single_cell_analysis()
        self.run_comparison_transcriptomics_data()
        self.run_supplementary_analysis()


if __name__ == "__main__":
    ProteomicsExperiment("Proteomics_study").run()
