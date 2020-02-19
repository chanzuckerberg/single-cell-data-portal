"""
(2/19/20) Temporary script used to load a preliminary artifact CSV sample.
          To be succeeded by ETL job(s) from a JSON artifact.
"""

import os
import pandas
import sys
from math import ceil

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.config.db_config import BrowserDbConfig

stage = os.environ["DEPLOYMENT_STAGE"]
db_name = f"browser_{stage}"

Base = declarative_base()
engine = create_engine(BrowserDbConfig().database_uri)
conn = engine.connect()
engine.execute(f"USE {db_name}")

# project
projects_df = pandas.read_csv("projects.csv")
for index, row in projects_df.iterrows():
    project_id = row["project_id"]
    title = row["project_title"]
    label = "Single cell transcriptome analysis of human pancreas"
    description = (
        "As organisms age, cells accumulate genetic and "
        "epigenetic changes that eventually lead to "
        "impaired organ function or catastrophic failure "
        "such as cancer. Here we describe a single-cell "
        "transcriptome analysis of 2544 human pancreas "
        "cells from donors, spanning six decades of life. "
        "We find that islet cells from older donors have "
        "increased levels of disorder as measured both "
        "by noise in the transcriptome and by the number "
        "of cells which display inappropriate hormone "
        "expression, revealing a transcriptional instability "
        "associated with aging. By analyzing the spectrum "
        "of somatic mutations in single cells from "
        "previously-healthy donors, we find a specific "
        "age-dependent mutational signature characterized "
        "by C to A and C to G transversions, indicators "
        "of oxidative stress, which is absent in single "
        "cells from human brain tissue or in a tumor cell "
        "line. Cells carrying a high load of such mutations "
        "also express higher levels of stress and senescence "
        "markers, including FOS, JUN, and the cytoplasmic "
        "superoxide dismutase SOD1, markers previously linked "
        "to pancreatic diseases with substantial age-dependent "
        "risk, such as type 2 diabetes mellitus and adenocarcinoma. "
        "Thus, our single-cell approach unveils gene expression "
        "changes and somatic mutations acquired in aging human "
        "tissue, and identifies molecular pathways induced by "
        "these genetic changes that could influence human "
        "disease. Also, our results demonstrate the feasibility "
        "of using single-cell RNA-seq data from primary "
        "cells to derive meaningful insights into the "
        "genetic processes that operate on aging human "
        "tissue and to determine which molecular mechanisms "
        "are coordinated with these processes. Examination "
        "of single cells from primary human pancreas tissue"
    )
    category = row["category"]
    developmental_stage = row["developmental_stage"]
    disease_ontology = row["disease_ontology"]
    sample_type = "specimens"
    organ_part = "islet of Langerhans"
    analysis_protocol = "smartseq2_v2.3.0,smartseq2_v2.4.0"
    cell_count = 2544
    donor_count = 8
    publication_title = row["publication_title"]
    publication_doi = row["publication_doi"]
    contact_name = row["contributor_name"]
    contact_institution = row["contributor_lab"]
    contact_email = row["contributor_email"]

    organ_ontology = row["organ_ontology"]
    species = row["species"]
    external_accessions = row["external_accessions"]
    library_construction_method_ontology = row["library_construction_method_ontology"]
    nucleic_acid_source = row["nucleic_acid_source"]
    end_bias = row["end_bias"]

    engine.execute(
        f"INSERT INTO project VALUES ("
        f"'{project_id}', '{title}', '{label}', '{description}', '{category}', "
        f"'{developmental_stage}', '{disease_ontology}', '{sample_type}', "
        f"'{organ_part}', '{analysis_protocol}', {cell_count}, {donor_count}, "
        f"'{publication_title}', '{publication_doi}', '{contact_name}', "
        f"'{contact_institution}', '{contact_email}')"
    )

    # library prep protocol
    engine.execute(
        f"INSERT INTO library_prep_protocol VALUES ("
        f"NULL, '{library_construction_method_ontology}', "
        f"'{end_bias}', '{nucleic_acid_source}')"
    )

    # tissue
    engine.execute(f"INSERT INTO tissue VALUES (NULL, '{organ_ontology}')")

    # species
    engine.execute(f"INSERT INTO species VALUES (NULL, '{species}')")

    # data repository

    # contributor

    # external accession

    # lpp x project
    engine.execute(
        f"INSERT INTO library_prep_protocol_join_project VALUES "
        f"(NULL, 1, '{project_id}')"
    )

    # tissue x project
    engine.execute(f"INSERT INTO tissue_join_project VALUES (NULL, 1, '{project_id}')")

    # species x project
    engine.execute(f"INSERT INTO species_join_project VALUES (NULL, 1, '{project_id}')")

    # contributor x project


# file
files_df = pandas.read_csv("files.csv")
files_df["file_size"].fillna(0, inplace=True)
files_df.fillna("", inplace=True)
files_df.insert(
    1,
    "project_id",
    ["HCA-Project-cddab57b-6868-4be4-806f-395ed9dd635a"] * files_df.shape[0],
)

items = []
for index, row in files_df.iterrows():
    items.append(
        ", ".join([f"'{str(v)}'" if isinstance(v, str) else str(v) for v in row.array])
    )

chunk_size = 100
num_chunks = ceil(len(items) / chunk_size)
for i in range(num_chunks):
    chunk = items[
        i * chunk_size : len(items)
        if len(items) < i * chunk_size + chunk_size
        else i * chunk_size + chunk_size
    ]
    values = ", ".join([f"({item})" for item in chunk])
    engine.execute(f"INSERT INTO file VALUES {values}")
    print(i * chunk_size)
