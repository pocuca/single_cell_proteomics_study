from pathlib import Path
from typing import Dict, List

import anndata
import pandas as pd
from anndata import AnnData
from numpy import nan

import bioexperiment.data.loader_definitions as ldef


def read_protein_data(data_folder_path: Path) -> Dict[str, AnnData]:
    """
    This function reads the two tsv files containing the raw protein dataset
    as well as the protein dataset with local regression normalization.

    Args:
        data_folder_path (Path): The path where the data is saved

    Returns:
        Dict[str, AnnData]: AnnData object
    """
    raw_data = read_tsv_file(data_folder_path, ldef.ProteinData.raw_file)
    norm_data = read_tsv_file(
        data_folder_path, ldef.ProteinData.normalized_file, fill_filtered_or_nan=nan
    )
    return {"raw": raw_data, "norm": norm_data}


def read_smartseq_data(data_folder_path: Path) -> AnnData:
    """
    This function loads and concatenates data from "HeLa-CCL2 cell heterogeneity
    studied by single-cell DNA and RNA sequencing"
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129447
    (SMART-Seq2 data from 9, 14 and 20 passages)
    Args:
        data_folder_path (Path): The path where the data is saved

    Returns:
        AnnData: AnnData object
    """
    data = [
        read_annotated_txt_file(data_folder_path, filename)
        for filename in ldef.Smartseq.files
    ]
    # TODO: point out this is bad
    data[2].var_names = data[0].var_names
    data[1] = data[1][:, data[0].var_names]  # without ERCC genes
    return anndata.concat(data)


def read_dropseq_data(data_folder_path: Path) -> AnnData:
    """
    This function loads and concatenates data from "The transcriptome
    dynamics of single cells during the cell cycle"
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142277

    Args:
        data_folder_path (Path): the path where the data is saved

    Returns:
        AnnData: AnnData object
    """
    ann_data_read = [
        read_intron_exon_txt_file(
            data_folder_path, file_pair["exon"], file_pair["intron"]
        )
        for file_pair in ldef.Dropseq.files.values()
    ]
    return anndata.concat(ann_data_read, join="outer", fill_value=0.0)


def read_tsv_file(
    data_folder_path: Path,
    filename: str,
    col_keyword_filter: str = "TSCP_DIA_SingleCell",
    annotation_column: str = "Genes",
    fill_filtered_or_nan: float = 0.0,
) -> anndata.AnnData:
    """
    This function reads a .tsv file, preprocess and returns an annotated data object.

    Args:
        data_folder_path (str): Path under which file can be found
        filename (str): The name of the file to be read and saved into memory
        col_keyword_filter: Substring used to filter columns in the file
        annotation_column: Column from annotation csv to be used
    Returns:
        _type_: AnnData matrix
    """
    cols = list(pd.read_csv(data_folder_path / filename, nrows=1, sep="\t"))
    prot = pd.read_csv(
        data_folder_path / filename,
        usecols=[i for i in cols if col_keyword_filter in i],
        sep="\t",
    )
    var_annotation = pd.read_csv(
        data_folder_path / filename,
        index_col=annotation_column,
        usecols=[annotation_column],
        sep="\t",
    )
    prot = prot.fillna(fill_filtered_or_nan)
    prot[prot == "Filtered"] = fill_filtered_or_nan

    adata_prot_raw = anndata.AnnData(prot.transpose())
    adata_prot_raw.var = var_annotation
    adata_prot_raw.var_names = adata_prot_raw.var_names.map(str)
    adata_prot_raw.var_names_make_unique()

    return adata_prot_raw


def read_annotated_txt_file(data_folder_path: Path, filename: str) -> AnnData:
    """
    This function reads specified .txt files and returns an annotated data object

    Args:
        data_folder_path (Path): Path to the data that is being read
        filename (str): The name of the file to be read

    Returns:
        AnnData: AnnData object
    """
    return anndata.read_text(data_folder_path / filename).transpose()


def read_intron_exon_txt_file(
    data_folder_path: Path, exon_filename: str, intron_filename: str
) -> anndata.AnnData:
    """
    This function reads specified .txt intron and exon files and returns
    an annotated data object

    Args:
        data_folder_path (str): Path under which file can be found
        exon_filename (str): The name of the exon file to be read and
            saved into memory
        intron_filename (str): The name of the intron file to be read and
            saved into memory

    Returns:
        _type_: AnnData matrix
    """
    adata_exon = read_annotated_txt_file(data_folder_path, exon_filename)
    adata_intron = read_annotated_txt_file(data_folder_path, intron_filename)

    intersection = list(set(adata_exon.var_names).intersection(adata_intron.var_names))
    intersection_part = adata_exon[:, intersection]
    intersection_part.X += adata_intron[:, intersection].X

    adata_exon_intron = anndata.concat(
        [
            intersection_part,
            adata_exon[
                :,
                [
                    gene
                    for gene in adata_exon.var_names
                    if gene not in adata_intron.var_names
                ],
            ],
            adata_intron[
                :,
                [
                    gene
                    for gene in adata_intron.var_names
                    if gene not in adata_exon.var_names
                ],
            ],
        ],
        axis=1,
    )

    return adata_exon_intron


def read_proteome_file(data_folder_path: Path, to_be_read: str = "Genes") -> List[str]:
    """
    Reading a txt proteome file.

    Args:
        filename (str): Name of file to be read
        to_be_read (str): Name of column to be read

    Returns:
        List[str]: a list of core proteomes read from file
    """
    dir_core_proteome = ldef.CoreProteome.file

    core_proteins = list(
        pd.read_csv(data_folder_path / dir_core_proteome, sep="\t")[to_be_read]
    )

    return core_proteins
