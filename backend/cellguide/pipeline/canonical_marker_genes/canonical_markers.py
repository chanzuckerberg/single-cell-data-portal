import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Tuple

import pandas as pd
from dask import compute, delayed
from dask.diagnostics import ProgressBar

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
from backend.cellguide.pipeline.constants import ASCTB_MASTER_SHEET_URL, CELLGUIDE_PIPELINE_NUM_CPUS
from backend.cellguide.pipeline.utils import get_gene_id_to_name_and_symbol
from backend.wmg.data.ontology_labels import ontology_term_label
from backend.wmg.data.utils import setup_retry_session

logger = logging.getLogger(__name__)


def get_asctb_master_sheet():
    session = setup_retry_session()
    asctb_data_response = session.get(ASCTB_MASTER_SHEET_URL)
    asctb_data_response.raise_for_status()
    return asctb_data_response.json()


class CanonicalMarkerGenesCompiler:
    def __init__(self, *, wmg_tissues: list[str], wmg_human_genes: list[str]):
        """
        Initializes the CanonicalMarkerGenesCompiler class.

        Arguments
        ---------
        wmg_tissues (list[str]): List of tissue ontology term IDs present in WMG.
        wmg_human_genes (list[str]): List of gene ensembl IDs present in WMG.
        """

        logger.info("Fetching ASCTB data...")
        self.asctb_data = get_asctb_master_sheet()

        # WMG tissues have some terms that start with "CL" because they're cell cultures.
        # We filter these out as they are specific to our platform and don't exist in ASCTB.
        self.wmg_tissues = [i for i in wmg_tissues if i.startswith("UBERON:")]

        self.wmg_human_genes = wmg_human_genes

        gene_metadata = get_gene_id_to_name_and_symbol()
        self.gene_id_to_name = gene_metadata.gene_id_to_name

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

    def _get_gene_info(self, gene_biomarkers: list[GeneBiomarker]) -> Tuple[str, str]:
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
        return ";;".join(gene_symbols), ";;".join(gene_names)

    def _get_references(self, references: list[Reference], doi_to_citation: dict[str, str]) -> Tuple[str, str]:
        """
        Extracts the DOIs and citations from the given list of references.
        This function cleans the DOI and gets the title and formatted citation from the DOI.

        Arguments
        ---------
        references - list[Reference]
            a list of references for a particular ASCTB table entry.
        doi_to_citation - dict[str,str]
            a dictionary mapping DOIs to citations.

        Returns
        -------
        refs - str
            The ';;'-concatenated DOIs extracted from the input references.
        titles - list[str]
            The ';;'-concatenated titles extracted from the input references.
        """

        def fetch_doi_info(ref):
            doi = clean_doi(ref.doi)
            if doi:
                if doi not in doi_to_citation:
                    title = get_title_and_citation_from_doi(doi)
                    doi_to_citation[doi] = title
                else:
                    title = doi_to_citation[doi]
                return doi, title

        with ThreadPoolExecutor() as executor:
            results = list(executor.map(fetch_doi_info, references))

        filtered_results = [result for result in results if result is not None]
        if not filtered_results:
            return "", ""

        refs, titles = zip(*filtered_results)
        return ";;".join(refs), ";;".join(titles)

    def get_processed_asctb_table_entries(self) -> dict[str, ParsedAsctbTableEntry]:
        """
        Processes the ASCTB table entries and returns a dictionary mapping cell type ontology term IDs to
        ParsedAsctbTableEntry objects. The processing involves extracting relevant information from the ASCTB data
        such as tissue ID, gene symbols and names, DOIs and citations, and cell type ontology term IDs. It also
        involves cleaning and formatting the extracted data.

        Returns
        -------
        dict[str, ParsedAsctbTableEntry]
            A dictionary mapping cell type ontology term IDs to ParsedAsctbTableEntry objects.
        """

        # DOI to citation mapping
        results = [delayed(self._process_asct_table__parallel)(tissue) for tissue in self.asctb_data]
        logger.info(f"Getting processed ASCTB table entries for {len(self.asctb_data)} tissues...")
        with ProgressBar():
            parsed_table_entries = sum(compute(*results, num_workers=CELLGUIDE_PIPELINE_NUM_CPUS), [])
        # Drop duplicate entries if they exist
        parsed_table_entries = pd.DataFrame(parsed_table_entries).drop_duplicates().to_dict("records")
        logger.info("Fetching tissue names from UBERON ontology...")
        tissues_in_parsed_table_entries = [i["tissue"] for i in parsed_table_entries]
        tissues_by_id = {t: ontology_term_label(t) for t in tissues_in_parsed_table_entries}

        gene_infos = {}
        for entry in parsed_table_entries:
            entry_copy = entry.copy()
            ct = entry_copy["cell_type_ontology_term_id"]
            del entry_copy["cell_type_ontology_term_id"]

            celltype_markers_list = gene_infos.get(ct, [])
            entry_copy["tissue"] = tissues_by_id.get(entry_copy["tissue"], entry_copy["tissue"])
            celltype_markers_list.append(entry_copy)
            gene_infos[ct] = celltype_markers_list

        for ct in gene_infos:
            gene_infos[ct] = [ParsedAsctbTableEntry(**entry) for entry in gene_infos[ct]]

        return gene_infos

    def _process_asct_table__parallel(self, tissue: str) -> list[dict[str, str]]:
        """
        Processes the ASCTB table entries for a given tissue in parallel and returns a list of dictionaries.
        Each dictionary is an entry containing a tissue, gene symbol, gene name, publication,
        publication title, and cell type ontology term ID. The processing involves extracting
        relevant information from the ASCTB data such as tissue ID, gene symbols and names, DOIs and citations,
        and cell type ontology term IDs. It also involves cleaning and formatting the extracted data.

        Arguments
        ---------
        tissue - str
            The tissue for which the ASCTB table entries are to be processed.

        Returns
        -------
        list[dict[str, str]]
            A list of dictionaries, each containing information about a unique combination of tissue, gene symbol,
            gene name, publication, publication title, and cell type ontology term ID.
        """

        doi_to_citation = {}

        data = self.asctb_data[tissue]["data"]

        parsed_table_entries = []
        for row in data:
            cell_types = [celltype["id"] for celltype in row["cell_types"] if celltype["id"].startswith("CL:")]
            if not cell_types or not row["biomarkers_gene"]:
                continue

            tissue_id = self._get_tissue_id([AnatomicalStructure(**entry) for entry in row["anatomical_structures"]])
            gene_symbols, gene_names = self._get_gene_info([GeneBiomarker(**entry) for entry in row["biomarkers_gene"]])
            refs, titles = self._get_references([Reference(**entry) for entry in row["references"]], doi_to_citation)

            for cell_type in cell_types:
                entry = {
                    "tissue": tissue_id,
                    "symbol": gene_symbols,
                    "name": gene_names,
                    "publication": refs,
                    "publication_titles": titles,
                    "cell_type_ontology_term_id": cell_type,
                }
                parsed_table_entries.append(entry)
        return parsed_table_entries
