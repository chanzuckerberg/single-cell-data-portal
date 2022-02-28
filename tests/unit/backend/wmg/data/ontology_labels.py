import re
import unittest

from backend.wmg.data import ontology_labels


class OntologyLabelTests(unittest.TestCase):
    def test__load_ontologies(self):
        ontology_term_id_labels = ontology_labels.load_ontologies()

        ontology_prefixes = {term_id.split(':')[0] for term_id in ontology_term_id_labels.keys()}

        # A conservative, high-level check that all ontologies have been loaded, without checking explicit counts
        self.assertEqual(ontology_prefixes,
                         {'HANCESTRO', 'HsapDv', 'PATO', 'EFO', 'NCBITaxon', 'CL', 'MONDO', 'MmusDv', 'UBERON'})

    def test__load_genes(self):
        gene_term_id_labels = ontology_labels.load_genes()

        gene_prefixes = {re.match('([A-Za-z]+)', term_id)[1] for term_id in gene_term_id_labels.keys()}

        # A conservative, high-level check that all ontologies have been loaded, without checking explicit counts
        self.assertEqual(gene_prefixes,
                         {'ENSSAST', 'ENSMUSG', 'ENSG', 'ENST', 'ENSSASG', 'ERCC', 'ENSMUST'})


if __name__ == '__main__':
    unittest.main()
