from typing import Dict, List, Tuple

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import seaborn as sns
import sklearn.metrics as metrics
from anndata import AnnData

from bioexperiment.data.preprocess import impute_downshifted_normal_global
from bioexperiment.visualizations.config import MatplotlibConfig, color


def cell_cycle_stage_box_plot(data: AnnData, plot_path: str):
    """
    This function plots a boxplot.

    Args:
        data (AnnData): Data whose datapoints are plotted
        ploth_path (str): The directory path where the plot will be saved
    """

    size = {
        cond: data[data.obs["cell cycle stage"] == cond].X.sum(axis=1)
        for cond in ["G1", "G1-S", "G2", "G2-M"]
    }
    x = np.concatenate([[cond] * len(cells) for cond, cells in size.items()])
    y = np.concatenate([cells for cells in size.values()])
    sns.boxplot(x=x, y=y)
    plt.xlabel("Cell cycle stage")
    plt.ylabel("Raw MS signal")

    # Write descriptions to file so we do not overcrowd the image
    text_box = []
    for cond, d in size.items():
        text_box.append(f"{cond}: n = {len(d)}, median = {np.median(d)}")

    dump_text("\n".join(text_box), f"{plot_path}.description.txt")
    plt.savefig(plot_path)
    plt.close()


def dump_text(text: str, path: str):
    """
    Writes description text from plots in a file

    Args:
        text (str): The description
        path (str): The path to the file
    """
    with open(path, "w") as f:
        f.write(text)


def pca_plot(
    data: AnnData,
    color_column: str,
    plot_path: str,
    impute_downshifted: bool = True,
    invert_axis: bool = False,
    **kwargs,
):
    """
    A function that builds a Principal Component Analysis plot

    Args:
        data (AnnData): AnnData object
        color_column (str): The color scheme the plot will follow
        ploth_path (str): The directory path where the plot will be saved
        impute_downshifted (bool, optional): Does the function need to be used or not. Defaults to True.
        invert_axis (bool, optional): Should the axes be inverted. Defaults to False.
    """

    if impute_downshifted:
        impute_downshifted_normal_global(data, **kwargs)

    sc.pp.pca(data)
    ax = sc.pl.pca(data, color=color_column, annotate_var_explained=True, show=False)

    if invert_axis:
        ax.invert_xaxis()

    plt.savefig(plot_path)
    plt.close()


def cell_cycle_stage_prediction_plot(
    data: AnnData,
    classes_to_compare: Tuple[str, str],
    plot_path: str,
    score_keys: List[str] = ["G1", "S", "G2-M"],
):
    """Building of a cell cycle stage visualization plot

    Args:
        data (AnnData): AnnData object
        classes_to_compare (Tuple[str, str]): Classess of Genes which will be compared in the plot
        ploth_path (str): The directory path where the plot will be saved
        score_keys (List[str], optional): Defaults to ["G1", "S", "G2-M"].
    """

    protein_data_nb_tb = data[
        np.array(data.obs["cell cycle stage"] == classes_to_compare[0])
        | np.array(data.obs["cell cycle stage"] == classes_to_compare[1])
    ]

    text_description = []
    for phase in score_keys:
        true_label = np.array(
            [phase in a for a in protein_data_nb_tb.obs["cell cycle stage"].values]
        )
        scores = protein_data_nb_tb.obs[phase].values
        fp, tp, _ = metrics.roc_curve(true_label, scores)
        auc_text = (
            f"AUC using the {phase} score: "
            f"{metrics.roc_auc_score(true_label, scores):.2f}"
        )
        text_description.append(auc_text)
        plt.plot(fp, tp, label=phase)
    plt.plot([0, 1], [0, 1], c="black")
    plt.gca().set_aspect("equal")
    plt.xlabel("False-positive")
    plt.ylabel("True-positive")
    plt.legend()
    plt.title(f"{classes_to_compare[0]} vs {classes_to_compare[1]}")
    plt.savefig(plot_path)
    plt.close()

    dump_text("\n".join(text_description), f"{plot_path}.description.txt")


def completeness_distribution_plot(data: Dict[str, AnnData], plot_path: str):
    """
    This function plots a Violin plot.

    Args:
        data (AnnData): Data whose datapoints are plotted
        ploth_path (str): The directory path where the plot will be saved
    """

    perc_expr = []
    dataset_name = []
    for datasetname, data in data.items():
        perc_expr += list((data.X > 0).mean(axis=1))  # type: ignore
        dataset_name += [datasetname] * data.X.shape[0]  # type: ignore
    sns.violinplot(y=perc_expr, x=dataset_name)
    plt.ylabel("Gene/Protein expression completeness per cell")
    plt.savefig(plot_path)
    plt.close()


# flake8: noqa: C901
def pearson_corr_analysis_plot(
    full_data: Dict[str, AnnData],
    plot_path: str,
):
    """
    Function that visualized the Pearson Correlation analysis

    Args:
        full_data (Dict[str, AnnData]): All three datasets saved in a dictionary
        ploth_path (str): The directory path where the plot will be saved
    """
    for adata in full_data.values():
        np.random.shuffle(adata.X)

    joined = anndata.concat(full_data.values())  # performs inner join
    dataset_names = np.concatenate(
        [[name] * adata.shape[0] for name, adata in full_data.items()]
    )

    # print('number of shared genes:', joined.shape[1])
    c = np.array(pd.DataFrame(joined.X.T).corr())

    fig, ax = plt.subplots(
        3, 3, figsize=(5, 5), gridspec_kw=dict(wspace=0.03, hspace=0.03)
    )
    for i, dataset_1 in enumerate(full_data.keys()):
        for j, dataset_2 in enumerate(full_data.keys()):
            subdata = c[dataset_names == dataset_1][:, dataset_names == dataset_2]
            color_obj = ax[i, j].imshow(
                subdata, vmin=-1, vmax=1, aspect="auto", cmap=MatplotlibConfig.newcmp
            )
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])
            if i == 0:
                ax[i, j].set_title(dataset_2, fontsize=11)
            if j == 0:
                ax[i, j].set_ylabel(dataset_1, fontsize=11)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
    plt.colorbar(color_obj, cax=cbar_ax)
    plt.savefig(plot_path)
    plt.close()


def variance_analysis_plot(core_proteins: List, data: AnnData, plot_path: str):
    """
    A plot that analyses the variance of proteins and shows the results in a log 2 scale

    Args:
        core_proteins (List): Core proteins from the core proteome file
        data (AnnData): AnnData object
        ploth_path (str): The directory path where the plot will be saved
    """
    data.var["core/noncore"] = [
        "Core proteome genes" if prot in core_proteins else "Other proteins"
        for prot in data.var_names
    ]
    mean = np.nanmean(data.X, axis=0)
    coeff_of_var = scipy.stats.variation(data.X, nan_policy="omit").data

    g = sns.JointGrid(
        x=np.log2(mean),
        y=np.log2(coeff_of_var),
        hue=data.var["core/noncore"],
        palette=color,
    )
    g.plot_joint(sns.scatterplot, s=10, alpha=1, linewidth=0.1)
    g.plot_marginals(
        sns.histplot,
        kde=True,
        stat="density",
        common_norm=False,
        binwidth=0.2,
        alpha=0.8,
    )
    g.ax_joint.set_xlabel("log 2 mean")
    g.ax_joint.set_ylabel("log 2 coefficient of variation")
    g.ax_joint.set_xlim(-4, 18)
    g.ax_joint.set_ylim(-4, 4)
    g.ax_joint.plot([-6, 10], [3, -5], "--", c="black", label="Poisson distribution")
    g.ax_joint.legend(title=[])
    plt.savefig(plot_path)
    plt.close()


def variance_comparison_box_plot(
    core_proteins: List, full_data: Dict[str, AnnData], plot_path: str
):
    """
    A box plot that showes the comparison analysis between the three datasets

    Args:
        core_proteins (List): The core proteins from the proteome
        full_data (AnnData): AnnData object
        plot_path (str): The directory path where the plot will be saved
    """
    cov = []
    hue = []
    dataset_name = []

    for datasetname, adata in full_data.items():
        adata.var["core/noncore"] = [
            "Core proteome genes" if prot in core_proteins else "Other proteins"
            for prot in adata.var_names
        ]
        coeff_of_var = scipy.stats.variation(adata.X, nan_policy="omit").data
        cov += list(coeff_of_var)
        hue += list(adata.var["core/noncore"])
        dataset_name += [datasetname] * adata.shape[1]

    sns.boxplot(
        x=dataset_name,
        y=cov,
        hue=hue,
        palette=color,
        showfliers=False,
    )
    plt.ylabel("coefficient of variation")
    plt.ylim(top=6)
    plt.savefig(plot_path)
    plt.close()


def local_regression_norm_plot(data_raw: AnnData, data_norm: AnnData):
    """
    Visualization of the local regression normalization

    Args:
        data_raw (AnnData): AnnData raw data object
        data_norm (AnnData): AnnData raw data object
    """
    x = data_raw.X.sum(axis=1)
    y = (data_raw.X > 0).sum(axis=1)
    plt.scatter(x, y, s=1)
    plt.xlabel("Raw MS-Signal per cell")
    plt.ylabel("Identified proteins per cell")
    plt.show()

    x = data_norm.X.sum(axis=1)
    y = (data_norm.X > 0).sum(axis=1)
    plt.scatter(x, y, s=1)
    plt.xlabel("Normalized MS-Signal per cell")
    plt.ylabel("Identified proteins per cell")
    plt.show()


def ranked_protein_by_completenes_plot(data: AnnData, full_data: AnnData):
    """
    Plot ranking the proteins by completeness

    Args:
        data (AnnData): AnnData object
        full_data (AnnData): All the three datasets
    """
    for name, data in full_data.items():
        completeness = (data.X > 0).sum(axis=0) / data.X.shape[0]
        sorted_completeness = np.sort(completeness)[::-1]
        plt.plot(
            np.arange(len(sorted_completeness)) + 1, sorted_completeness, label=name
        )

    plt.legend()
    plt.ylabel("Gene/Protein completeness")
    plt.xlabel("Ranked protein (by completeness)")
    plt.xscale("log")
    plt.show()


def log_mean_protein_abundance_plot(data: AnnData, full_data: AnnData):
    """
    Log_mean_protein_abundance_plot

    Args:
        data (AnnData): AnnData object
        full_data (AnnData): All three datasets used for the plot
    """
    for datasetname, data in full_data.items():
        x = np.nanmean(data.X, axis=0)
        y = np.nanmean(data.X > 0, axis=0)
        plt.scatter(np.log(x), y, s=1, label=datasetname)

    x2 = np.arange(0, 1000000, 0.1)
    plt.plot(np.log(x2), 1 - np.exp(-x2), c="black", label="expected poisson dropout")
    plt.xlabel("log mean protein abundance")
    plt.ylabel("Completeness accross cells")
    plt.legend(markerscale=5.0, loc="upper left")
    plt.show()
