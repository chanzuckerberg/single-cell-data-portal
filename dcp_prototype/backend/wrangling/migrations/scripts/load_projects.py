from argparse import ArgumentParser
import sys

sys.path.insert(0, "")  # noqa
from dcp_prototype.backend.wrangling.migrations.scripts.squish_files import squish_files, reset_counts

projects = [
        "HumanMousePancreas",
        "1M Neurons",
        "MouseGastrulationAtlas",
        "HumanTissueTcellActivation",
        "WongAdultRetina",
        "HPSI human cerebral organoids",
        "BM_PC",
        "Multiplexed scRNA-seq with barcoded antibodies",
        "1M Immune Cells",
        "Drop-seq, DroNc-seq, Fluidigm C1 Comparison",
        "SingleCellLiverLandscape",
        "HDCA-Sweden-10x",
        "TissueStability",
        "Mouse Melanoma",
        "CD4+ cytotoxic T lymphocytes",
        "Human Hematopoietic Profiling",
        "Fetal/Maternal Interface",
        "Kidney biopsy scRNA-seq",
        "Reprogrammed_Dendritic_Cells",
        "Single cell transcriptome analysis of human pancreas",
        "Single cell RNAseq characterization of cell types produced over time in an in vitro model of human inhibitory interneuron differentiation",
        "Healthy and type 2 diabetes pancreas",
        "KidneySingleCellAtlas",
        "HumanColonicMesenchymeIBD",
        "scRNAseqSystemicComparison",
        "Tabula Muris",
        "Diabetic Nephropathy snRNA-seq",
        "Mouse Endoderm Project",
        "snRNA-seq_for_human_retina"]

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=True,
        help="A data directory containing bundles with files that contain both the "
        "metadata and the data of a single project.",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        nargs="+",
        required=True,
        help="An output directory to which the files from the input directory will be "
        "copied to and de-duped. File names will be retained except for the suffix"
        " that is a numbering of the file.",
    )
    parser.add_argument(
        "-c",
        "--count_file",
        nargs="+",
        required=True,
        help="A count file that keeps track of the number of a particular file that the"
        " script has processed. This file allows for multiple runs of the script "
        "while retaining previous counts.",
    )

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    if arguments.output_directory:
        output_directory = arguments.output_directory[0]
    if arguments.count_file:
        count_file = arguments.count_file[0]
    for project in projects:
        project = project.replace(" ", "-")
        full_input_directory = f"{input_directory}/{project}/bundles"
        full_output_directory = f"{output_directory}/{project}"
        reset_counts(count_file)
        squish_files(full_input_directory, full_output_directory, count_file)

