import enum
import gzip


class SupportedOrganisms(enum.Enum):
    HOMO_SAPIENS = "NCBITaxon:9606"
    MUS_MUSCULUS = "NCBITaxon:10090"
    SARS_COV_2 = "NCBITaxon:2697049"
    ERCC = "NCBITaxon:32630"


class GeneChecker:
    """Handles checking gene ids, retrieves symbols"""

    base_prefix = "backend/common/ontology_files/"
    GENE_FILES = {
        SupportedOrganisms.HOMO_SAPIENS: base_prefix + "genes_homo_sapiens.csv.gz",
        SupportedOrganisms.MUS_MUSCULUS: base_prefix + "genes_mus_musculus.csv.gz",
        SupportedOrganisms.SARS_COV_2: base_prefix + "genes_sars_cov_2.csv.gz",
        SupportedOrganisms.ERCC: base_prefix + "genes_ercc.csv.gz",
    }
    gene_dict = {}  # type: ignore

    def __init__(self):
        if not GeneChecker.gene_dict:
            for file in self.GENE_FILES:
                with gzip.open(self.GENE_FILES[file], "rt") as genes:
                    for gene in genes:
                        gene = gene.rstrip().split(",")
                        gene_id = gene[0]
                        gene_label = gene[1]
                        GeneChecker.gene_dict[gene_label] = gene_id

    def get_id(self, gene_label) -> str:
        """
        Gets ENSEMBL id associated to the gene label

        :param str gene_label: gene symbol

        :rtype str
        :return A gene ENSEMBL id
        """
        if not gene_label:
            raise ValueError("No gene label provided.")
        # in processing, if there are duplicate gene labels,
        # a gene might be gene_label + "_" + gene_id
        if "_" in gene_label:
            gene_label = gene_label.split("_")[1]

        if gene_label not in GeneChecker.gene_dict:
            # could be an id already
            if gene_label in GeneChecker.gene_dict.values():
                return gene_label
            raise ValueError(f"Could not find '{gene_label}' in dictionary.")
        return GeneChecker.gene_dict[gene_label]
