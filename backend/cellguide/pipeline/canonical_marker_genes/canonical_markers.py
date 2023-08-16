import json
import logging
import pathlib
from typing import Dict, Tuple

import pandas as pd
import requests

from backend.cellguide.pipeline.canonical_marker_genes.types import (
    AnatomicalStructure,
    GeneBiomarker,
    ParsedAsctbTableEntry,
    Reference,
)
from backend.cellguide.pipeline.canonical_marker_genes.utils import (
    clean_doi,
    get_title_and_citation_from_doi,
)
from backend.cellguide.pipeline.constants import ASCTB_MASTER_SHEET_URL, ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME
from backend.wmg.data.ontology_labels import ontology_term_label

logger = logging.getLogger(__name__)


class CanonicalMarkerGenesCompiler:
    def __init__(self, wmg_tissues: list[str], wmg_human_genes: list[str]):
        """
        Initializes the CanonicalMarkerGenesCompiler class.

        Arguments
        ---------
        wmg_tissues (list[str]): List of tissue ontology term IDs present in WMG.
        wmg_human_genes (list[str]): List of gene ensembl IDs present in WMG.
        """

        logger.info("Fetching ASCTB data...")
        self.asctb_data = requests.get(ASCTB_MASTER_SHEET_URL).json()
        self.wmg_tissues = [i for i in wmg_tissues if i.startswith("UBERON:")]

        self.wmg_human_genes = wmg_human_genes

        file_path = self._get_symbol_to_gene_descriptions_file_path()
        gene_id_to_name = pd.read_csv(file_path, sep="\t")
        self.gene_id_to_name = gene_id_to_name.set_index("Symbols")["Description"].to_dict()

    def _get_symbol_to_gene_descriptions_file_path(self) -> pathlib.Path:
        """
        Returns the file path for the table mapping gene symbols to gene descriptions.

        Returns
        -------
        pathlib.Path: The file path for the ensembl to gene descriptions file.
        """

        return (
            pathlib.Path(__file__)
            .parent.absolute()
            .parent.joinpath("fixtures", ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME)
        )

    def _get_tissue_id(self, anatomical_structures: list[AnatomicalStructure]) -> str:
        """
        Extracts the tissue ID from the given anatomical structures. By convention,
        the anatomical structures from ASCTB are sorted in order of increasing granularity.
        This function prioritizes tissue IDs that are present in WMG. It takes the most
        granular term that is present in WMG. Otherwise, it defaults to the least
        granular term.


        Arguments
        ---------
        anatomical_structures - list[AnatomicalStructure]
            a list of anatomical structures for a particular ASCTB table entry.

        Returns
        -------
        tissue_id - str
            The tissue ID extracted from the input anatomical structures.
        """

        tissue_id = anatomical_structures[0].id
        for entry in anatomical_structures:
            if entry.id in self.wmg_tissues:
                tissue_id = entry.id
        return tissue_id

    def _get_gene_info(self, gene_biomarkers: list[GeneBiomarker]) -> Tuple[list[str], list[str]]:
        """
        Extracts the gene information from the given gene biomarkers.
        This function only adds gene IDs that are present in WMG.

        Arguments
        ---------
        gene_biomarkers - list[GeneBiomarker]
            a list of gene biomarkers for a particular ASCTB table entry.

        Returns
        -------
        gene_symbols - list[str]
            The list of gene symbols extracted from the input gene biomarkers.
        gene_names - list[str]
            The list of gene names extracted from the input gene biomarkers.
        """

        gene_symbols = []
        gene_names = []
        for gene in gene_biomarkers:
            symbol = gene.name.upper()
            if symbol and symbol in self.wmg_human_genes:
                name = gene.rdfs_label or self.gene_id_to_name.get(symbol, symbol)
                gene_symbols.append(symbol)
                gene_names.append(name)
        return gene_symbols, gene_names

    def _get_references(
        self, references: list[Reference], doi_to_citation: Dict[str, str]
    ) -> Tuple[list[str], list[str]]:
        """
        Extracts the DOIs and citations from the given list of references.
        This function cleans the DOI and gets the title and formatted citation from the DOI.

        Arguments
        ---------
        references - list[Reference]
            a list of references for a particular ASCTB table entry.
        doi_to_citation - Dict[str,str]
            a dictionary mapping DOIs to citations.

        Returns
        -------
        refs - list[str]
            The list of cleaned DOIs extracted from the input references.
        titles - list[str]
            The list of titles extracted from the input references.
        """

        refs = []
        titles = []
        for ref in references:
            doi = clean_doi(ref.doi)
            if doi:
                if doi not in doi_to_citation:
                    title = get_title_and_citation_from_doi(doi)
                    doi_to_citation[doi] = title
                else:
                    title = doi_to_citation[doi]
                refs.append(doi)
                titles.append(title)
        return ";;".join(refs), ";;".join(titles)

    def get_processed_asctb_table_entries(self) -> Dict[str, ParsedAsctbTableEntry]:
        """
        Processes the ASCTB table entries and returns a dictionary mapping cell type ontology term IDs to
        ParsedAsctbTableEntry objects. The processing involves extracting relevant information from the ASCTB data
        such as tissue ID, gene symbols and names, DOIs and citations, and cell type ontology term IDs. It also
        involves cleaning and formatting the extracted data.

        Returns
        -------
        Dict[str, ParsedAsctbTableEntry]
            A dictionary mapping cell type ontology term IDs to ParsedAsctbTableEntry objects.
        """

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

                tissue_id = self._get_tissue_id(
                    [AnatomicalStructure(**entry) for entry in row["anatomical_structures"]]
                )
                gene_symbols, gene_names = self._get_gene_info(
                    [GeneBiomarker(**entry) for entry in row["biomarkers_gene"]]
                )
                refs, titles = self._get_references(
                    [Reference(**entry) for entry in row["references"]], doi_to_citation
                )

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
            gene_infos[cell_type] = [ParsedAsctbTableEntry(**entry) for entry in aggregated_gene_info]

        return gene_infos
