"""
6/12/20 calvinnhieu

Dirty script to load sample matrix (H5AD, Remix) into DB.
This dataset will be used in the short-term to:
1. Verify the DB Schema
2. Provide testing data to unblock API development

The sample dataset can be found in the
`cellxgene-data-wrangling-prod` bucket
in the `single-cell-prod` account.

This script will be deprecated when file upload is implemented in Corpora.
"""

import json
import os
import scanpy as sc
import sys
import uuid

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, pkg_root)

from corpora.common.corpora_orm import (
    DBSessionMaker,
    ProjectStatus,
    ProcessingState,
    ValidationState,
    ProjectLinkType,
    DatasetArtifactType,
    DatasetArtifactFileType,
    DbUser,
    DbProject,
    DbProjectDataset,
    DbProjectLink,
    DbDataset,
    DbDatasetArtifact,
    DbContributor,
    DbDatasetContributor,
)

session = DBSessionMaker().session()

filename = "Single_cell_gene_expression_profiling_of_SARS_CoV_2_infected_human_cell_lines_Calu_3-28-remixed.h5ad"
adata = sc.read(filename)

# admin user
user_id = str(uuid.uuid4())
user_name = "Admin"
email = "admin@example.com"

user = DbUser(id=user_id, name=user_name, email=email)
print(user)
session.add(user)
session.commit()

# project
project_id = str(uuid.uuid4())
status = ProjectStatus.LIVE.name
owner = user_id
project_name = adata.uns["project_title"]
description = ""  # no field
s3_bucket = f"s3://corpora-data-dev/{project_id}"
tc_uri = ""
needs_attestation = False
processing_state = ProcessingState.NA.name
validation_state = ValidationState.NOT_VALIDATED.name

project = DbProject(
    id=project_id,
    status=status,
    owner=owner,
    name=project_name,
    description=description,
    s3_bucket=s3_bucket,
    tc_uri=tc_uri,
    needs_attestation=needs_attestation,
    processing_state=processing_state,
    validation_state=validation_state,
)
print(project)
session.add(project)
session.commit()

# dataset
dataset_id = str(uuid.uuid4())
revision = 1
dataset_name = adata.uns["title"]
organism = adata.uns["organism"]
organism_ontology = adata.uns["organism_ontology"]
tissue = adata.obs["tissue"].dtype.categories[0]
tissue_ontology = adata.obs["tissue_ontology"].dtype.categories[0]
assay = adata.obs["assay"].dtype.categories[0]
assay_ontology = adata.obs["assay_ontology"].dtype.categories[0]
disease = adata.obs["disease"].dtype.categories[0]
disease_ontology = adata.obs["disease_ontology"].dtype.categories[0]
sex = "NA"  # Not available. Required?
ethnicity = "NA"  # Not available. Required?
ethnicity_ontology = "NA"  # Not available. Required?
source_data_location = "NA"  # Not available. Required?
preprint_doi = adata.uns["preprint_doi"]
publication_doi = "NA"  # Not available. Required?

dataset = DbDataset(
    id=dataset_id,
    revision=revision,
    name=dataset_name,
    organism=organism,
    organism_ontology=organism_ontology,
    tissue=tissue,
    tissue_ontology=tissue_ontology,
    assay=assay,
    assay_ontology=assay_ontology,
    disease=disease,
    disease_ontology=disease_ontology,
    sex=sex,
    ethnicity=ethnicity,
    ethnicity_ontology=ethnicity_ontology,
    source_data_location=source_data_location,
    preprint_doi=preprint_doi,
    publication_doi=publication_doi,
)
print(dataset)
session.add(dataset)
session.commit()

# project dataset
pd_id = str(uuid.uuid4())
pd_project_id = project_id
pd_project_status = ProjectStatus.LIVE.name
pd_dataset_id = dataset_id

project_dataset = DbProjectDataset(
    id=pd_id, project_id=pd_project_id, project_status=pd_project_status, dataset_id=pd_dataset_id
)
print(project_dataset)
session.add(project_dataset)
session.commit()

# project links
# RAW DATA
project_links = []
for link in adata.uns["project_raw_data_links"]:
    project_link_id = str(uuid.uuid4())
    pl_p_id = project_id
    pl_p_status = ProjectStatus.LIVE.name
    link_url = link
    link_type = ProjectLinkType.RAW_DATA.name

    project_link = DbProjectLink(
        id=project_link_id, project_id=pl_p_id, project_status=pl_p_status, link_url=link_url, link_type=link_type
    )
    print(project_link)
    project_links.append(project_link)
# OTHER
for link in adata.uns["project_other_links"]:
    project_link_id = str(uuid.uuid4())
    pl_p_id = project_id
    pl_p_status = ProjectStatus.LIVE.name
    link_url = link
    link_type = ProjectLinkType.OTHER.name

    project_link = DbProjectLink(
        id=project_link_id, project_id=pl_p_id, project_status=pl_p_status, link_url=link_url, link_type=link_type
    )
    print(project_link)
    project_links.append(project_link)
session.add_all(project_links)
session.commit()

# dataset artifact
dataset_artifact_id = str(uuid.uuid4())
da_d_id = dataset_id
dataset_filename = filename
dataset_filetype = DatasetArtifactFileType.H5AD.name
dataset_type = DatasetArtifactType.REMIX.name
user_submitted = True
s3_uri = f"s3://corpora-api-dev/matrix/remix/{filename}"

dataset_artifact = DbDatasetArtifact(
    id=dataset_artifact_id,
    dataset_id=da_d_id,
    filename=dataset_filename,
    filetype=dataset_filetype,
    type=dataset_type,
    user_submitted=user_submitted,
    s3_uri=s3_uri,
)
print(dataset_artifact)
session.add(dataset_artifact)
session.commit()

# deployment directory. DNE at this stage.

# contributor
# dataset contributor
db_contributors = []
dcs = []
for contributor_str in adata.uns["contributors"]:
    contributor = json.loads(contributor_str.replace("'", '"'))
    contributor_id = str(uuid.uuid4())
    contributor_name = contributor["name"]
    institution = contributor["institution"]
    contributor_email = contributor["email"]

    db_contributor = DbContributor(
        id=contributor_id, name=contributor_name, institution=institution, email=contributor_email
    )
    print(db_contributor)
    db_contributors.append(db_contributor)

    dataset_contributor_id = str(uuid.uuid4())
    dc_c_id = contributor_id
    dc_d_id = dataset_id

    dataset_contributor = DbDatasetContributor(id=dataset_contributor_id, contributor_id=dc_c_id, dataset_id=dc_d_id)
    print(dataset_contributor)
    dcs.append(dataset_contributor)
session.add_all(db_contributors)
session.commit()
session.add_all(dcs)
session.commit()
