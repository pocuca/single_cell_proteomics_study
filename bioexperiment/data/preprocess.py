#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def impute_downshifted_normal_global(
    adata: AnnData, scale: float = 0.3, shift: float = 1.8, random_seed=42
) -> AnnData:
    """
    A function that prepares data for a Principal Component Analysis plot.

    Args:
        adata (AnnData): AnnData object
        scale (float, optional): Scale value. Defaults to 0.3.
        shift (float, optional): Shift value. Defaults to 1.8.
        random_seed (int, optional): Random seed, for reproducibility. Defaults to 42.

    Returns:
        AnnData: AnnData object
    """
    np.random.seed(random_seed)
    mean = np.nanmean(adata.X)
    std = np.nanstd(adata.X)
    draws = np.random.normal(
        loc=mean - shift * std, scale=scale * std, size=np.sum(np.isnan(adata.X))
    )
    adata.X[np.isnan(adata.X)] = draws
    return adata


def apply_data_filters(
    data: AnnData,
    min_num_genes: Optional[int] = 600,
    min_part_cells: Optional[float] = 0.15,
    **kwargs
) -> AnnData:
    """
    Apply data filters. Filters are done inplace to prevent as much as possible
     copying the whole data.

    Args:
        data (AnnData): AnnData object
        min_num_genes (Optional[int], optional): The lower bound
        for number of genes. Defaults to 600.
        min_part_cells (Optional[float], optional): The lower bound
        for part cells. Defaults to 0.15.

    Returns:
        AnnData: filtered data
    """
    if min_num_genes:
        sc.pp.filter_cells(data, min_genes=min_num_genes, **kwargs)
    if min_part_cells:
        sc.pp.filter_genes(data, min_cells=min_part_cells * data.shape[0], **kwargs)

    return data


def apply_data_transformations(
    data: AnnData,
    log_transform_data: bool = False,
    fill_nans: bool = False,
    normalize: bool = False,
    norm_target_sum: Union[float, np.ndarray] = 1e6,
    copy_data: bool = True,
) -> AnnData:
    """
    This function prepares the data in the way to be simple to analyse and visualize

    Args:
        data (object): Data that gets filtered and prepared
        plot_type (str): The plot type needs to be specified

    Returns:
        AnnData: The filtered AnnData object

    """
    # Data is copied because all transformations will be done in
    # a new data structure, preventing the original full data to be modified.
    # In some cases, copy can be skipped if the idea is to modify the original data.
    if copy_data:
        data = data.copy()
    if normalize:
        sc.pp.normalize_total(data, target_sum=norm_target_sum)
    if log_transform_data:
        sc.pp.log1p(data)
    if fill_nans:
        data.X = np.nan_to_num(data.X)

    return data


def ann_data_map_from_observation_names(
    data: AnnData, output_column_name: str, keywords_dictionary: Dict[str, str]
) -> AnnData:
    """
    Mapping AnnData from names of observations

    Args:
        data (AnnData): AnnData object
        output_column_name (str): New column name
        keywords_dictionary (Dict[str, str]): Dictionary keywords used for the mapping

    Returns:
        AnnData: AnnData object
    """
    # Defaults in case of no mapping
    name_to_map = {name: "other" for name in data.obs_names}
    for name in data.obs_names:
        for keyword, mapped_name in keywords_dictionary.items():
            if keyword in name:
                name_to_map[name] = mapped_name
                break
    data.obs[output_column_name] = data.obs_names.map(name_to_map)

    return data


def build_cell_cycle_markers(
    data: AnnData, cc_markers: pd.DataFrame
) -> Dict[str, list]:
    """
    A function that builds markers for the cell cycle stage prediction plot

    Args:
        data (AnnData): AnnData whose variable names are being used for the mappign
        cc_markers (pd.DataFrame): Markers being used for the mapping

    Returns:
        Dict[str, list]: A dictionary
    """
    g1_list = cc_markers["G1"][1:-1]
    g1_list = [g for g in g1_list if g in data.var_names]

    s_list = cc_markers["S"][1:-1]
    s_list = [g for g in s_list if g in data.var_names]

    g2m_list = cc_markers["G2M"][1:-1]
    g2m_list = [g for g in g2m_list if g in data.var_names]

    return {"G1": g1_list, "S": s_list, "G2-M": g2m_list}


def filter_out_category(
    data: AnnData, column_name: str, filter_out_categories: List[str]
):
    """Function that filters out the cathegory 'other' before plotting

    Args:
        data (AnnData): AnnData object
        column_name (str): The name of the column to be filtered
        filter_out_categories (List[str]): The categories to be filtered out

    Returns:
        _type_: the filtered data
    """
    filtered_data = data
    for category in filter_out_categories:
        filtered_data = filtered_data[data.obs[column_name] != category]
    return filtered_data


def score_genes_with_cell_markers(data: AnnData, markers: Dict[str, list]) -> AnnData:
    """
    Obtaining the scores for specific genes that are later used for a plot

    Args:
        data (AnnData): AnnData object
        markers (Dict[str, list]): Markers used to compute the score

    Returns:
        AnnData: AnnData object
    """
    for marker in markers:
        sc.tl.score_genes(data, markers[marker], score_name=marker)
    return data


def get_data_intersection(full_data: Dict[str, AnnData]) -> AnnData:
    """
    Intersecting the shared genes for a PCA plot

    Args:
        full_data (Dict[str, AnnData]): The dataset being passed to be intersected

    Returns:
        AnnData: AnnData object
    """
    shared_genes = set(full_data["Proteins"].var_names)
    for adata in full_data.values():
        shared_genes = shared_genes.intersection(adata.var_names)

    intersection_data = {
        datasetname: adata[:, list(shared_genes)]
        for datasetname, adata in full_data.items()
    }

    return intersection_data
