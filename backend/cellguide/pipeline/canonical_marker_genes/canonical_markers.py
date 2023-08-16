import json
import logging

import pandas as pd
import requests

from backend.cellguide.pipeline.canonical_marker_genes.utils import (
    clean_doi,
    get_gene_name_from_gene_info,
    get_tissue_name_from_uberon,
    get_title_and_citation_from_doi,
)
from backend.cellguide.pipeline.constants import ASCTB_MASTER_SHEET_URL
from backend.wmg.data.snapshot import WmgSnapshot

logger = logging.getLogger(__name__)


class CanonicalMarkerGenesCompiler:
    def __init__(self, snapshot: WmgSnapshot):
        self.asctb_data = requests.get(ASCTB_MASTER_SHEET_URL).json()

        cell_counts_df = snapshot.cell_counts_df[:]
        self.wmg_tissues = [i for i in cell_counts_df["tissue_ontology_term_id"].unique() if i.startswith("UBERON:")]

        self.wmg_human_genes = [
            next(iter(i.values())) for i in snapshot.primary_filter_dimensions["gene_terms"]["NCBITaxon:9606"]
        ]

        self.doi_to_citation = self._get_formatted_citations_per_doi()

    def _get_formatted_citations_per_doi(self):
        # get formatted citations per doi
        doi_to_citation = {}
        for i, tissue in enumerate(self.asctb_data):
            logger.info(f"Getting formatted citations for DOIs in {tissue} tissue ({i}/{len(self.asctb_data)})")

            data = self.asctb_data[tissue]["data"]
            for row in data:
                cell_types = [i["id"] for i in row["cell_types"] if i["id"].startswith("CL:")]
                if len(cell_types) and len(row["biomarkers_gene"]):
                    for ref in row["references"]:
                        doi = clean_doi(ref["doi"])

                        if doi != "" and doi not in doi_to_citation:
                            title = get_title_and_citation_from_doi(doi)
                            doi_to_citation[doi] = title
        return doi_to_citation

    def _get_tissue_id(self, row):
        tissue_id = row["anatomical_structures"][0]["id"]
        for entry in row["anatomical_structures"]:
            if entry["id"] in self.wmg_tissues:
                tissue_id = entry["id"]
        return tissue_id

    def _get_gene_info(self, row, gene_name_memory):
        gene_symbols = []
        gene_names = []
        for gene in row["biomarkers_gene"]:
            symbol = gene["name"].upper()
            if symbol and symbol in self.wmg_human_genes:
                name = gene["rdfs_label"] or self._get_gene_name(symbol, gene_name_memory)
                gene_symbols.append(symbol)
                gene_names.append(name)
        return gene_symbols, gene_names

    def _get_gene_name(self, symbol, gene_name_memory):
        if symbol in gene_name_memory:
            return gene_name_memory[symbol]
        name = get_gene_name_from_gene_info(symbol)
        gene_name_memory[symbol] = name
        return name

    def _get_references(self, row):
        refs = []
        titles = []
        for ref in row["references"]:
            doi = clean_doi(ref["doi"])
            if doi:
                title = self.doi_to_citation[doi]
                refs.append(doi)
                titles.append(title)
        return ";;".join(refs), ";;".join(titles)

    def get_processed_asctb_table_entries(self):
        gene_name_memory = {}
        hashed_entries_seen = []
        parsed_table_entries = []

        for i, tissue in enumerate(self.asctb_data):
            logger.info(f"Getting processed ASCTB table entry for {tissue} tissue ({i}/{len(self.asctb_data)})")
            data = self.asctb_data[tissue]["data"]

            for row in data:
                cell_types = [i["id"] for i in row["cell_types"] if i["id"].startswith("CL:")]
                if not cell_types or not row["biomarkers_gene"]:
                    continue

                tissue_id = self._get_tissue_id(row)
                gene_symbols, gene_names = self._get_gene_info(row, gene_name_memory)
                refs, titles = self._get_references(row)

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

        tissues_in_parsed_table_entries = [i["tissue"] for i in parsed_table_entries]
        tissues_by_id = {t: get_tissue_name_from_uberon(t) for t in tissues_in_parsed_table_entries}

        gene_infos = {}
        for entry in parsed_table_entries:
            entry = entry.copy()
            ct = entry["cell_type_ontology_term_id"]
            del entry["cell_type_ontology_term_id"]

            a = gene_infos.get(ct, [])
            entry["tissue"] = tissues_by_id.get(entry["tissue"], entry["tissue"])
            a.append(entry)
            gene_infos[ct] = a

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
