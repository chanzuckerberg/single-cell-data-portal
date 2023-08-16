import json
import logging
import pathlib

import pandas as pd
import requests

from backend.cellguide.pipeline.canonical_marker_genes.utils import (
    clean_doi,
    get_title_and_citation_from_doi,
)
from backend.cellguide.pipeline.constants import ASCTB_MASTER_SHEET_URL, ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME
from backend.wmg.data.ontology_labels import ontology_term_label

logger = logging.getLogger(__name__)


class CanonicalMarkerGenesCompiler:
    def __init__(self, wmg_tissues, wmg_human_genes):
        logger.info("Fetching ASCTB data...")
        self.asctb_data = requests.get(ASCTB_MASTER_SHEET_URL).json()
        self.wmg_tissues = [i for i in wmg_tissues if i.startswith("UBERON:")]

        self.wmg_human_genes = wmg_human_genes

        file_path = self._get_ensembl_to_gene_descriptions_file_path()
        gene_id_to_name = pd.read_csv(file_path, sep="\t")
        self.gene_id_to_name = gene_id_to_name.set_index("Symbols")["Description"].to_dict()

    def _get_ensembl_to_gene_descriptions_file_path(self):
        return (
            pathlib.Path(__file__)
            .parent.absolute()
            .parent.joinpath("fixtures", ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME)
        )

    def _get_tissue_id(self, row):
        tissue_id = row["anatomical_structures"][0]["id"]
        for entry in row["anatomical_structures"]:
            if entry["id"] in self.wmg_tissues:
                tissue_id = entry["id"]
        return tissue_id

    def _get_gene_info(self, row):
        gene_symbols = []
        gene_names = []
        for gene in row["biomarkers_gene"]:
            symbol = gene["name"].upper()
            if symbol and symbol in self.wmg_human_genes:
                name = gene["rdfs_label"] or self.gene_id_to_name.get(symbol, symbol)
                gene_symbols.append(symbol)
                gene_names.append(name)
        return gene_symbols, gene_names

    def _get_references(self, row, doi_to_citation):
        refs = []
        titles = []
        for ref in row["references"]:
            doi = clean_doi(ref["doi"])
            if doi:
                if doi not in doi_to_citation:
                    title = get_title_and_citation_from_doi(doi)
                    doi_to_citation[doi] = title
                else:
                    title = doi_to_citation[doi]
                refs.append(doi)
                titles.append(title)
        return ";;".join(refs), ";;".join(titles)

    def get_processed_asctb_table_entries(self):

        # DOI to citation mapping
        doi_to_citation = {}

        hashed_entries_seen = []
        parsed_table_entries = []

        for i, tissue in enumerate(self.asctb_data):
            version = self.asctb_data[tissue]["metadata"]["version"]
            logger.info(
                f"Getting processed ASCTB table entry for {tissue} tissue ({version}) ({i+1}/{len(self.asctb_data)})"
            )
            data = self.asctb_data[tissue]["data"]

            for row in data:
                cell_types = [i["id"] for i in row["cell_types"] if i["id"].startswith("CL:")]
                if not cell_types or not row["biomarkers_gene"]:
                    continue

                tissue_id = self._get_tissue_id(row)
                gene_symbols, gene_names = self._get_gene_info(row)
                refs, titles = self._get_references(row, doi_to_citation)

                for cell_type in cell_types:
                    for i in range(len(gene_symbols)):
                        symbol = gene_symbols[i]
                        name = gene_names[i]

                        entry = {
                            "tissue": tissue_id,
                            "symbol": symbol,
                            "name": name,
                            "publication": refs,
                            "publication_titles": titles,
                            "cell_type_ontology_term_id": cell_type,
                        }
                        hashed_dict = hash(json.dumps(entry))
                        if hashed_dict not in hashed_entries_seen:
                            parsed_table_entries.append(entry)
                            hashed_entries_seen.append(entry)

        logger.info("Fetching tissue names from UBERON ontology...")
        tissues_in_parsed_table_entries = [i["tissue"] for i in parsed_table_entries]
        tissues_by_id = {t: ontology_term_label(t) for t in tissues_in_parsed_table_entries}

        gene_infos = {}
        for entry in parsed_table_entries:
            entry = entry.copy()
            ct = entry["cell_type_ontology_term_id"]
            del entry["cell_type_ontology_term_id"]

            a = gene_infos.get(ct, [])
            entry["tissue"] = tissues_by_id.get(entry["tissue"], entry["tissue"])
            a.append(entry)
            gene_infos[ct] = a

        logger.info("Aggregating gene biomarkers across tissues and publications across biomarkers and cell types...")

        def aggregate_publications(publication):
            # Remove empty publications and get unique ones with order preserved
            publications = [i for i in publication.values if i != ""]
            unique_publications = list(dict.fromkeys(publications))
            return ";;".join(unique_publications)

        def aggregate_gene_names(names, symbols):
            # If possible, pick a gene name that is not a gene symbol
            non_symbol_names = [name for name in names.values if name not in symbols]
            return next(iter(non_symbol_names), names.values[0])

        for cell_type in gene_infos:
            # Convert the list of gene info for the current cell type to a DataFrame
            gene_info_df = pd.DataFrame(gene_infos[cell_type])
            gene_symbols = gene_info_df["symbol"].values

            # Group the DataFrame by tissue and symbol, and aggregate the names and publications
            aggregated_gene_info = (
                gene_info_df.groupby(["tissue", "symbol"])
                .agg(
                    {
                        "name": lambda names, symbols=gene_symbols: aggregate_gene_names(names, symbols),
                        "publication": aggregate_publications,
                        "publication_titles": aggregate_publications,
                    }
                )
                .reset_index()
                .to_dict(orient="records")
            )

            # Create a DataFrame from the aggregated gene info
            aggregated_gene_info_df = pd.DataFrame(aggregated_gene_info)

            # Aggregate the gene info across all tissues
            all_tissues_gene_info = (
                aggregated_gene_info_df.groupby("symbol")
                .agg(
                    {
                        "name": lambda names, symbols=gene_symbols: aggregate_gene_names(names, symbols),
                        "publication": aggregate_publications,
                        "publication_titles": aggregate_publications,
                    }
                )
                .reset_index()
            )

            # Add a column indicating these records are for all tissues
            all_tissues_gene_info["tissue"] = "All Tissues"

            # Ensure the DataFrame has the same column order as the original
            all_tissues_gene_info = all_tissues_gene_info[aggregated_gene_info_df.columns]

            # Add the all-tissues gene info to the list of aggregated gene info
            aggregated_gene_info.extend(all_tissues_gene_info.to_dict(orient="records"))

            # Update the gene info for the current cell type
            gene_infos[cell_type] = aggregated_gene_info

        return gene_infos
