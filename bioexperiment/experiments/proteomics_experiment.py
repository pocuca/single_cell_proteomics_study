from functools import cached_property
from typing import Dict, List

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
    score_genes_with_cell_markers,
)
from bioexperiment.visualizations.plots import (
    cell_cycle_stage_box_plot,
    cell_cycle_stage_prediction_plot,
    pca_plot,
)


class ProteomicsExperiment(Experiment):

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
    # of the data structures used
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

    @property
    def protein_data(self) -> AnnData:
        return self.full_protein_data["norm"]

    @property
    def raw_protein_data(self) -> AnnData:
        return self.full_protein_data["raw"]

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
                    fill_nans=False,
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
                    copy_data=False,
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

    def run_comparison_transcriptomics_data(self):
        pass

    def run_supplementary_analysis(self):
        pass

    def run(self):
        self.preprocess_data()
        self.run_true_single_cell_analysis()
        self.run_comparison_transcriptomics_data()
        self.run_supplementary_analysis()


ProteomicsExperiment("test").run()
