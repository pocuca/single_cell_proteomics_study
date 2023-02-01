from pathlib import Path

root = Path(__file__).parent.parent


class ProteinData:
    raw_file = "20210919_DIANN_SingleCellOutput.pg_matrix_notnormalized.tsv"
    normalized_file = "20210919_DIANN_SingleCellOutput.pg_matrix_cellcyclepred.tsv"


class Smartseq:
    files = [
        "GSE129447_RAW/GSM3713084_HeLa_1.txt",
        "GSE129447_RAW/GSM3713085_HeLa_2.txt",
        "GSE129447_RAW/GSM3713086_HeLa_3.txt",
    ]


class Dropseq:
    files = {
        1: {
            "exon": "GSE142277_RAW/"
            "GSM4224315_out_gene_exon_tagged.dge_exonssf002_WT.txt",
            "intron": "GSE142277_RAW/"
            "GSM4224315_out_gene_exon_tagged.dge_intronssf002_WT.txt",
        },
        2: {
            "exon": "GSE142277_RAW/"
            "GSM4224316_out_gene_exon_tagged.dge_exonssf002_KO.txt",
            "intron": "GSE142277_RAW/"
            "GSM4224316_out_gene_exon_tagged.dge_intronssf002_KO.txt",
        },
        3: {
            "exon": "GSM4226257_RAW/"
            "GSM4226257_out_gene_exon_tagged.dge_exonsds_046.txt",
            "intron": "GSM4226257_RAW/"
            "GSM4226257_out_gene_exon_tagged.dge_intronsds_046.txt",
        },
    }


DATA_FOLDER = root / "data" / "experiment_data"
