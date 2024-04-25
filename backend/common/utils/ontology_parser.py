from cellxgene_ontology_guide.ontology_parser import OntologyParser

from backend.wmg.data.constants import WMG_PINNED_SCHEMA_VERSION

ontology_parser = OntologyParser(schema_version=f"v{WMG_PINNED_SCHEMA_VERSION}")
