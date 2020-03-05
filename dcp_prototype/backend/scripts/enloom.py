"""
Create a loom file from a GEO matrix.

This is more of an ad hoc POC, but outlines the steps.
"""

import csv
import loompy
import sys

import pandas as pd

# Read in gene annotation information. mart_export is from EBI's
# biomart. This is a mouse dataset, so the mouse genes were selected
gene_info = {}
with open("mart_export.txt") as gene_f:
    reader = csv.DictReader(gene_f, delimiter="\t")
    for row in reader:
        gene_info[row["Gene name"]] = {
            "accession": row["Gene stable ID"],
            "chromosome": row["Chromosome/scaffold name"],
            "start": row["Gene start (bp)"],
            "end": row["Gene end (bp)"],
            "type": row["Gene type"],
        }

expr_df = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=0)

cell_id = expr_df.columns.to_numpy()

# Figure out the gene annotations
accession = []
gene_names = []
chromosome = []
featurestart = []
featureend = []
featuretype = []
genus_species = []

counter = 0
for gene_name in expr_df.index:
    gene = gene_info.get(gene_name)

    if not gene:
        print("Missing gene annotation info for", gene_name)
        accession.append("Not found " + str(counter))
        gene_names.append(gene_name)
        chromosome.append("N/A")
        featurestart.append("N/A")
        featureend.append("N/A")
        featuretype.append("N/A")
        genus_species.append("Mus musculus")
        counter += 1
    else:
        accession.append(gene["accession"])
        gene_names.append(gene_name)
        chromosome.append(gene["chromosome"])
        featurestart.append(gene["start"])
        featureend.append(gene["end"])
        featuretype.append(gene["type"])
        genus_species.append("Mus musculus")

# And now the cell annotations. GEO has a few formats for its metadata, this
# is the one called "matrix".
with open("GSE60361_series_matrix.txt") as gene_f:
    for line in gene_f:
        if line.startswith("!Sample_description"):
            cell_ids = [eval(c) for c in line.lstrip("!Sample_description\t").split("\t")]
        elif line.startswith("!Sample_characteristics_ch1") and "Sex:" in line:
            sexes = [eval(c).split(": ")[-1] for c in line.lstrip("!Sample_description_ch1\t").split("\t")]
        elif line.startswith("!Sample_characteristics_ch1") and "age:" in line:
            ages = [eval(c).split(": ")[-1] for c in line.lstrip("!Sample_description_ch1\t").split("\t")]
        elif line.startswith("!Sample_relation") and "SAMN" in line:
            donor_accessions = [eval(c).split("/")[-1] for c in line.lstrip("!Sample_relation\t").split("\t")]
        elif line.startswith("!Sample_source_name_ch1") and "sscortex" in line:
            organ_names = [eval(c) for c in line.lstrip("!Sample_source_name_ch1\t").split("\t")]
            organs = []
            # These are more of an educated guess
            for o in organ_names:
                if o == "sscortex":
                    organs.append("UBERON:0008930")
                elif o == "ca1hippocampus":
                    organs.append("UBERON:0003881")

    cell_metadata = {
        c[0]: {"sex": c[1], "age": c[2], "donor_accession": c[3], "organ": c[4]}
        for c in zip(cell_ids, sexes, ages, donor_accessions, organs)
    }

biosample_category = ["primary" for _ in cell_id]
col_attrs = {
    "CellID": cell_id,
    "biosample_preparation.biosample_category": ["primary" for _ in cell_id],
    "biosample_preparation.donor_accession": [cell_metadata[c]["donor_accession"] for c in cell_id],
    "biosample_preparation.donor_development_stage_at_collection": ["MmusDv:0000112" for _ in cell_id],
    "biosample_preparation.donor_diseases": ["PATO:0000461" for _ in cell_id],
    "biosample_preparation.donor_strain": ["EFO:0005180" for _ in cell_id],
    "biosample_preparation.donor_sex": [cell_metadata[c]["sex"] for c in cell_id],
    "biosample_preparation.donor_species": ["Mus musculus" for c in cell_id],
    "biosample_preparation.organ": [cell_metadata[c]["organ"] for c in cell_id],
    "library_preparation_protocol.library_construction_method": ["EFO:0009990" for _ in cell_id],
    "library_preparation_protocol.nucleic_acid_source": ["OBI:0000880" for _ in cell_id],
    "project.project_title": ["Single-cell RNA-seq of mouse cerebral cortex" for c in cell_id],
}


loompy.create(
    "matrix.loom",
    expr_df.to_numpy(),
    col_attrs=col_attrs,
    row_attrs={
        "Accession": accession,
        "Gene": gene_names,
        "chromosome": chromosome,
        "featurestart": featurestart,
        "featureend": featureend,
        "genus_species": genus_species,
    },
)
