from pathlib import Path
from typing import Dict, List

import anndata
import pandas as pd
from anndata import AnnData
from numpy import nan

import bioexperiment.data.loader_definitions as ldef


def read_smartseq_data(data_folder_path: Path) -> AnnData:
    data = [
        read_annotated_txt_file(data_folder_path, filename)
        for filename in ldef.Smartseq.files
    ]
    # TODO: point out this is bad
    data[2].var_names = data[0].var_names
    data[1] = data[1][:, data[0].var_names]  # without ERCC genes
    return anndata.concat(data)


def read_dropseq_data(data_folder_path: Path) -> AnnData:

    return anndata.concat(
        [
            read_intron_exon_txt_file(
                data_folder_path, file_pair["exon"], file_pair["intron"]
            )
            for file_pair in ldef.Dropseq.files.values()
        ]
    )


def read_protein_data(data_folder_path: Path) -> Dict[str, AnnData]:
    raw_data = read_tsv_file(data_folder_path, ldef.ProteinData.raw_file)
    norm_data = read_tsv_file(
        data_folder_path, ldef.ProteinData.normalized_file, fill_filtered_or_nan=nan
    )
    return {"raw": raw_data, "norm": norm_data}


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
        data_folder_path (str): path under which file can be found
        filename (str): the name of the file to be read and saved into memory
        col_keyword_filter: substring used to filter columns in the file
        annotation_column: column from annotation csv to be used
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
    return anndata.read_text(data_folder_path / filename).transpose()


def read_intron_exon_txt_file(
    data_folder_path: Path, exon_filename: str, intron_filename: str
) -> anndata.AnnData:
    """
    This function reads specified .txt files and returns an annotated data object

    Args:
        data_folder_path (str): path under which file can be found
        exon_filename (str): the name of the exon file to be read and
            saved into memory
        intron_filename (str): the name of the intron file to be read and
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


def read_proteome_file(filename: str, to_be_read: str = "Genes") -> List[str]:

    """Reading a txt proteome file.

    Args:
        filename (str): name of file to be read
        to_be_read (str): name of column to be read

    Returns:
        _type_: a list of core proteomes read from file
    """

    dir_core_proteome = filename
    core_proteins = list(pd.read_csv(dir_core_proteome, sep="\t")[to_be_read])

    return core_proteins
