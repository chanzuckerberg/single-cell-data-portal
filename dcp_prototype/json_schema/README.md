# The Human Cell Atlas Metadata Standards: JSON Schemas

This repo contains the HCA metadata JSON schemas for **Core**, **Type**, and **Module** entities. Below is a brief description of the structure of the HCA JSON schemas. More details about the schema structure can be found in the [structure document](../docs/structure.md).

## Schema structure overview

### User-supplied fields

Metadata fields are divided into three major entities: Core, Type, and Module.

Metadata fields in **Core** entity schemas are those that apply to 100% of all entities and are highly stable. Example Core fields include those for IDs and accessions. Core fields are inherited by Type entities.

Metadata fields in **Type** entity schemas are those that apply to a majority of entities and are moderately stable. Example Type fields include Height and Weight for donors and total estimated cell count for cell suspensions. 

Metadata fields in **Module** entity schemas are those that apply to a small number of entities or may be shared across different entities. Example Module fields include Cause of death for donors and Publication information (which can apply to both projects and biomaterials like cell lines). **Ontology** specifications are a special class of Module entity schemas.

### Ingest-supplied fields

Some metadata about a submitted dataset is provided by the [HCA Ingestion Service](https://github.com/HumanCellAtlas/ingest-central) including information about when a dataset was submitted and how all the Type entities relate to each other. The HCA JSON schemas describing these attributes are located in the `system/` directory. Data contributors do not need to know about these schemas.

## Schema structure rules

The following rules are followed by HCA JSON schemas:

1. Core schemas can not reference Type or Module schemas.
1. All Type schemas must reference their associated Core schema.
1. Module schemas can only reference Ontology schemas.
1. A set of fields that are logically grouped or only applicable in certain situations should be put into their own Module schema.
1. Core, Type, and Module schemas are grouped by their entity category (Biomaterial, File, Process, Project, Protocol).

## Directory structure of schemas

```
core/
    biomaterial/
        biomaterial_core.json 
    file/
        file_core.json
    process/
        process_core.json
    project/
        project_core.json
    protocol/ 
        protocol_core.json
module/
    biomaterial/
        cell_morphology.json
        ...
    ontology/ 
        biological_macromolecule_ontology.json
        ...
    process/ 
        purchased_reagents.json
        ...                
    project/
        contact.json
        ...
    protocol/ 
        channel.json
        ...  
type/
    biomaterial/
        cell_line.json
        ...
    file/
        analysis_file.json
        ...
    process/
        analysis/ 
            analysis_process.json
        ...
    project/ 
        project.json
    protocol/ 
        analysis/
            analysis_protocol.json
        ...
system/
    links.json
    provenance.json 
```
