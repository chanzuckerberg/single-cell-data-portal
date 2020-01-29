from utils.id_generator import hca_accession_generator

from typing import List


class Project:
    def __init__(
        self,
        contributor_names: List,
        lab,
        publication_title,
        publication_doi,
        external_accessions: List = [],
    ):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.contributor_names = contributor_names
        self.lab = lab
        self.publication_title = publication_title
        self.publication_doi = publication_doi
        self.external_accessions = external_accessions
