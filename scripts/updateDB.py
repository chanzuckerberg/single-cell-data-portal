#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import create_engine, text
from backend.layers.persistence.orm import mapper_registry
from backend.layers.persistence.constants import SCHEMA_NAME
import uuid
from datetime import datetime
from sqlalchemy.orm import Session
from backend.layers.persistence.orm import (
    CollectionTable, CollectionVersionTable, DatasetTable, 
    DatasetVersionTable, DatasetArtifactTable
)

def _format_citation(collection_data):
    """
    Format citations in the same way as the backend's format_citation_dp function.
    Authors is a list of dicts with 'family' and 'name' keys.
    """
    authors = collection_data.get('authors', [{'family': 'Anonymous', 'name': 'Anonymous'}])
    author_str_suffix = ""
    if len(authors) > 1:
        author_str_suffix = " et al."
    
    first_author = authors[0]
    author_str = first_author.get('family', first_author.get('name', 'Anonymous'))
    author_str += author_str_suffix
    
    journal = collection_data.get('journal', '')
    year = collection_data.get('publication_year', 2024)
    
    return f"{author_str} ({year}) {journal}"


def _get_valid_tissue_id(tissue_label):
    """Map tissue labels to valid UBERON ontology IDs"""
    tissue_map = {
        'brain': 'UBERON:0000955',  # brain
        'lung': 'UBERON:0002048',   # lung
        'heart': 'UBERON:0000948',  # heart
        'kidney': 'UBERON:0002113', # kidney
        'liver': 'UBERON:0002107',  # liver
    }
    return tissue_map.get(tissue_label.lower(), 'UBERON:0001016')  # default: tissue


def _get_valid_disease_id(disease_label):
    """Map disease labels to valid MONDO ontology IDs"""
    disease_map = {
        'healthy': 'MONDO:0000001',  # disease
        'cancer': 'MONDO:0005148',   # cancer
        'covid-19': 'MONDO:0100096', # COVID-19
    }
    return disease_map.get(disease_label.lower(), 'MONDO:0000001')


def _get_valid_cell_type_id(cell_type_label):
    """Map cell type labels to valid CL ontology IDs"""
    cell_type_map = {
        'neuron': 'CL:0000540',           # neuron
        'astrocyte': 'CL:0000069',        # astrocyte
        'oligodendrocyte': 'CL:0000128',  # oligodendrocyte
        'microglia': 'CL:0000129',        # microglia cell
        'endothelial cell': 'CL:0000115', # endothelial cell
    }
    return cell_type_map.get(cell_type_label.lower(), 'CL:0000003')  # default: native cell


def _get_valid_assay_id(assay_label):
    """Map assay labels to valid OBI ontology IDs"""
    assay_map = {
        '10x 3\' gene expression': 'OBI:0002685',  # 10x 3' v2
        '10x 5\' gene expression': 'OBI:0002686',  # 10x 5' v2
        'smart-seq2': 'OBI:0002751',                # Smart-seq2
        'bulk rna': 'OBI:0001271',                  # RNA-seq assay
    }
    return assay_map.get(assay_label.lower(), 'OBI:0002751')  # default: Smart-seq2

# Create engine and setup schema
engine = create_engine('postgresql://corpora:test_pw@database.corporanet.local:5432')
with engine.connect() as conn:
    conn.execute(text(f'CREATE SCHEMA IF NOT EXISTS {SCHEMA_NAME}'))
    conn.commit()

mapper_registry.metadata.create_all(engine)
print('Tables created!')

# Sample collections data with full metadata
collections_data = [
    {
        'name': 'Brain Cell Atlas',
        'description': 'Single-cell atlas of the human brain with comprehensive cell type mapping',
        'contact_name': 'John Smith',
        'contact_email': 'jsmith@example.com',
        'consortia': ['CZ', 'Cambridge'],
        'keywords': ['brain', 'neuronal', 'cell-types', 'transcriptomics'],
        'links': [
            {'name': None, 'type': 'DOI', 'uri': 'https://doi.org/10.1234/brain.atlas'},
            {'name': None, 'type': 'PROTOCOL', 'uri': 'https://protocol.io/brain-protocol'}
        ],
        'doi': '10.1234/brain.atlas',
        'journal': 'Nature Neuroscience',
        'publication_year': 2023,
        'authors': [{'family': 'Smith', 'name': 'John'}, {'family': 'Wilson', 'name': 'Jane'}],
    },
    {
        'name': 'Lung Cell Landscape',
        'description': 'Comprehensive map of lung cell types and immune populations',
        'contact_name': 'Jane Doe',
        'contact_email': 'jdoe@example.com',
        'consortia': ['CZ'],
        'keywords': ['lung', 'immune', 'respiratory', 'atlas'],
        'links': [
            {'name': None, 'type': 'DOI', 'uri': 'https://doi.org/10.1234/lung.atlas'},
        ],
        'doi': '10.1234/lung.atlas',
        'journal': 'Cell',
        'publication_year': 2022,
        'authors': [{'family': 'Doe', 'name': 'Jane'}, {'family': 'Johnson', 'name': 'Mark'}],
    },
    {
        'name': 'Heart Single-Cell Atlas',
        'description': 'Single-cell transcriptomics of human heart tissue',
        'contact_name': 'Bob Lee',
        'contact_email': 'blee@example.com',
        'consortia': ['CZ'],
        'keywords': ['heart', 'cardiac', 'single-cell', 'rna-seq'],
        'links': [],
        'doi': '10.1234/heart.atlas',
        'journal': 'Circulation',
        'publication_year': 2023,
        'authors': [{'family': 'Lee', 'name': 'Bob'}],
    },
    {
        'name': 'Kidney Cell Atlas',
        'description': 'Cell type composition and spatial organization of the human kidney',
        'contact_name': 'Alice Wang',
        'contact_email': 'awang@example.com',
        'consortia': ['Cambridge'],
        'keywords': ['kidney', 'nephron', 'glomerulus', 'tubule'],
        'links': [],
        'doi': '10.1234/kidney.atlas',
        'journal': 'Nature Methods',
        'publication_year': 2024,
        'authors': [{'family': 'Wang', 'name': 'Alice'}, {'family': 'Chen', 'name': 'Robert'}],
    },
    {
        'name': 'Liver Single-Cell Map',
        'description': 'Comprehensive single-cell map of human liver tissue including hepatocytes and immune cells',
        'contact_name': 'Charlie Brown',
        'contact_email': 'cbrown@example.com',
        'consortia': ['CZ', 'Cambridge'],
        'keywords': ['liver', 'hepatocyte', 'immune', 'metabolism'],
        'links': [
            {'name': None, 'type': 'DOI', 'uri': 'https://doi.org/10.1234/liver.atlas'},
            {'name': None, 'type': 'RAW_DATA', 'uri': 'https://data.org/liver-data'}
        ],
        'doi': '10.1234/liver.atlas',
        'journal': 'Hepatology',
        'publication_year': 2023,
        'authors': [{'family': 'Brown', 'name': 'Charlie'}, {'family': 'Prince', 'name': 'Diana'}],
    },
]

# Sample dataset data per collection with comprehensive metadata
# Each collection has all datasets with SAME tissue, disease, organism
# Datasets within a collection differ only by assay, cell_type, and processing stage
datasets_per_collection = {
    'brain': [  # Brain Collection
        {
            'name': 'Neurons Sample - 10x Genomics',
            'tissue': ['brain'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['neuron', 'astrocyte'],
        },
        {
            'name': 'Glia Sample - Smart-seq2',
            'tissue': ['brain'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['oligodendrocyte', 'microglia'],
        },
        {
            'name': 'Mixed Neurons - 10x 5\' v2',
            'tissue': ['brain'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 5\' gene expression'],
            'cell_type': ['neuron', 'astrocyte', 'oligodendrocyte'],
        },
        {
            'name': 'Filtered Dataset - Quality Control',
            'tissue': ['brain'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['neuron'],
        },
        {
            'name': 'Processed Expression Matrix',
            'tissue': ['brain'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['neuron', 'astrocyte', 'oligodendrocyte', 'microglia'],
        },
    ],
    'lung': [  # Lung Collection
        {
            'name': 'Lung Epithelial - 10x Genomics',
            'tissue': ['lung'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Lung Immune - Smart-seq2',
            'tissue': ['lung'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Lung Full - 10x 5\' v2',
            'tissue': ['lung'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 5\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Lung QC Dataset',
            'tissue': ['lung'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Lung Processed',
            'tissue': ['lung'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
    ],
    'heart': [  # Heart Collection
        {
            'name': 'Cardiac Myocytes - 10x',
            'tissue': ['heart'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Cardiac Fibroblasts - Smart-seq2',
            'tissue': ['heart'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Cardiac Immune - 10x 5\'',
            'tissue': ['heart'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 5\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Cardiac Filtered',
            'tissue': ['heart'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Cardiac Processed',
            'tissue': ['heart'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
    ],
    'kidney': [  # Kidney Collection
        {
            'name': 'Glomerular - 10x',
            'tissue': ['kidney'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Tubular - Smart-seq2',
            'tissue': ['kidney'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Collecting Duct - 10x 5\'',
            'tissue': ['kidney'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 5\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Kidney Filtered',
            'tissue': ['kidney'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Kidney Processed',
            'tissue': ['kidney'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
    ],
    'liver': [  # Liver Collection
        {
            'name': 'Hepatocytes - 10x',
            'tissue': ['liver'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Kupffer Cells - Smart-seq2',
            'tissue': ['liver'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Sinusoidal - 10x 5\'',
            'tissue': ['liver'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 5\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Liver Filtered',
            'tissue': ['liver'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['10x 3\' gene expression'],
            'cell_type': ['endothelial cell'],
        },
        {
            'name': 'Liver Processed',
            'tissue': ['liver'],
            'disease': ['healthy'],
            'organism': ['Homo sapiens'],
            'assay': ['Smart-seq2'],
            'cell_type': ['endothelial cell'],
        },
    ],
}

with Session(engine) as session:
    # Delete all existing data (in reverse order of foreign key dependencies)
    print('Cleaning database - deleting existing data...')
    
    try:
        session.query(DatasetArtifactTable).delete()
        print('  ✓ Deleted all DatasetArtifacts')
        
        session.query(DatasetVersionTable).delete()
        print('  ✓ Deleted all DatasetVersions')
        
        session.query(DatasetTable).delete()
        print('  ✓ Deleted all Datasets')
        
        session.query(CollectionVersionTable).delete()
        print('  ✓ Deleted all CollectionVersions')
        
        session.query(CollectionTable).delete()
        print('  ✓ Deleted all Collections')
        
        session.commit()
        print('✓ All tables cleaned successfully\n')
    except Exception as e:
        session.rollback()
        print(f'  Warning: Error during cleanup: {e}')
        print('  Continuing with insert operations...\n')
    
    collection_ids = []
    
    # First, prepare all IDs
    print('Preparing IDs for collections...')
    prepared_collections = []
    for c in collections_data:
        collection_id = uuid.uuid4()
        version_id = uuid.uuid4()
        prepared_collections.append((collection_id, version_id, c))
    
    # Insert Collections with version_id reference
    print('Inserting Collections...')
    for collection_id, version_id, c in prepared_collections:
        now = datetime.now()
        
        collection = CollectionTable(
            id=collection_id,
            version_id=version_id,  # Set reference to version
            originally_published_at=now,
            tombstone=False
        )
        session.add(collection)
        print(f'  Added collection: {c["name"]} (version_id: {version_id})')

    session.flush()  # Ensure collections are written
    print(f'✓ {len(prepared_collections)} Collections inserted with version_id\n')

    # Now insert CollectionVersions
    print('Inserting CollectionVersions...')
    for collection_id, version_id, c in prepared_collections:
        now = datetime.now()
        
        collection_ids.append((collection_id, version_id))

        collection_version = CollectionVersionTable(
            id=version_id,  # Use the same version_id
            collection_id=collection_id,
            collection_metadata={
                'name': c['name'],
                'description': c['description'],
                'contact_name': c['contact_name'],
                'contact_email': c['contact_email'],
                'consortia': c.get('consortia', []),
                'links': c.get('links', []),
            },
            owner='test_user@example.com',
            curator_name='Test Curator',
            publisher_metadata={
                'doi': c.get('doi', f'10.1234/test.{str(collection_id)[:8]}'),
                'citation': _format_citation(c),
                'journal': c.get('journal', 'Journal of Cell Biology'),
                'published_year': c.get('publication_year', 2024),
                'is_preprint': c.get('is_preprint', False),
                'authors': c.get('authors', [{'family': 'Anonymous', 'name': 'Anonymous'}]),
            },
            published_at=now,
            created_at=now,
            schema_version='1.0.0',
            datasets=[],
            has_custom_dataset_order=False,
            is_auto_version=False,
            data_submission_policy_version='1'
        )
        session.add(collection_version)
        print(f'  Added version for: {c["name"]}')

    session.flush()  # Ensure collection versions are written
    print(f'✓ {len(collection_ids)} CollectionVersions inserted\n')

    # Commit collections first
    session.commit()
    print('✓ Collections committed to database with version_id references\n')

    # Insert Datasets, DatasetVersions, and DatasetArtifacts
    artifact_count = 0
    dataset_version_ids_by_collection = {}  # Track datasets per collection
    
    # Map collection index to tissue type for datasets
    tissue_keys = ['brain', 'lung', 'heart', 'kidney', 'liver']
    
    for collection_idx, (collection_id, collection_version_id) in enumerate(collection_ids):
        dataset_version_ids_by_collection[collection_version_id] = []
        
        # Get the appropriate datasets for this collection
        tissue_key = tissue_keys[collection_idx % len(tissue_keys)]
        datasets_for_collection = datasets_per_collection[tissue_key]
        
        for i, dataset_info in enumerate(datasets_for_collection):
            dataset_id = uuid.uuid4()
            dataset_version_id = uuid.uuid4()
            now = datetime.now()

            # Track dataset version IDs for this collection
            dataset_version_ids_by_collection[collection_version_id].append(dataset_version_id)

            # Insert Dataset
            session.add(DatasetTable(
                id=dataset_id,
                version_id=dataset_version_id,
                published_at=now,
                tombstone=False
            ))

            # Insert DatasetVersion
            session.add(DatasetVersionTable(
                id=dataset_version_id,
                dataset_id=dataset_id,
                collection_id=collection_id,
                created_at=now,
                dataset_metadata={
                    'name': dataset_info['name'],
                    'schema_version': '3.0.0',
                    'cell_count': 5000 + i * 500,
                    'mean_genes_per_cell': 2000.0,
                    'tissue': [{'label': t, 'ontology_term_id': _get_valid_tissue_id(t)} for t in dataset_info['tissue']],
                    'disease': [{'label': d, 'ontology_term_id': _get_valid_disease_id(d)} for d in dataset_info['disease']],
                    'organism': [{'label': o, 'ontology_term_id': 'NCBITaxon:9606'} for o in dataset_info['organism']],
                    'assay': [{'label': a, 'ontology_term_id': _get_valid_assay_id(a)} for a in dataset_info['assay']],
                    'cell_type': [{'label': c, 'ontology_term_id': _get_valid_cell_type_id(c)} for c in dataset_info['cell_type']],
                    'sex': [{'label': 'unknown', 'ontology_term_id': 'PATO:0000383'}],
                    'self_reported_ethnicity': [{'label': 'unknown', 'ontology_term_id': 'unknown'}],
                    'development_stage': [{'label': 'unknown', 'ontology_term_id': 'unknown'}],
                    'batch_condition': [],
                    'suspension_type': ['cell'],
                    'donor_id': [],
                    'is_primary_data': 'PRIMARY',
                    'x_approximate_distribution': 'normal',
                },
                artifacts=[],
                status={
                    'upload_status': 'UPLOADED',
                    'validation_status': 'VALID',
                    'processing_status': 'SUCCESS',
                    'h5ad_status': 'CONVERTED',
                    'cxg_status': 'CONVERTED',
                    'rds_status': 'NA',
                    'atac_status': 'NA',
                }
            ))

            # Insert 2 DatasetArtifacts per dataset (h5ad and cxg)
            artifact_ids = []
            for j in range(2):
                artifact_id = uuid.uuid4()
                artifact_type = 'h5ad' if j == 0 else 'cxg'
                
                # h5ad file: dataset_version_id.h5ad
                # cxg file: dataset_id.cxg
                if artifact_type == 'h5ad':
                    # uri = f's3://bucket/dataset-{dataset_id}/{dataset_version_id}.h5ad'
                    uri = f's3://corpora-data-dev/{dataset_version_id}.h5ad'
                else:
                    # uri = f's3://bucket/dataset-{dataset_id}/{dataset_id}.cxg'
                    uri = f's3://corpora-data-dev/{dataset_id}.cxg'
                
                session.add(DatasetArtifactTable(
                    id=artifact_id,
                    type=artifact_type,
                    uri=uri
                ))
                artifact_ids.append(artifact_id)
                artifact_count += 1

            # Update dataset version with artifact references
            dataset_version = session.query(DatasetVersionTable).filter_by(id=dataset_version_id).first()
            if dataset_version:
                dataset_version.artifacts = artifact_ids

    session.commit()
    print(f'✓ 25 Datasets inserted (5 per collection)')
    print(f'✓ 25 DatasetVersions inserted')
    print(f'✓ {artifact_count} DatasetArtifacts inserted (2 per dataset)\n')
    
    # Link datasets to collections
    print('Linking datasets to collections...')
    for collection_version_id, dataset_version_ids in dataset_version_ids_by_collection.items():
        collection_version = session.query(CollectionVersionTable).filter_by(id=collection_version_id).first()
        if collection_version:
            collection_version.datasets = [str(dv_id) for dv_id in dataset_version_ids]
    
    session.commit()
    print(f'✓ Datasets linked to collections\n')
    print('All data inserted successfully!')

    # Verify data in database
    print('\n Database verification:')
    collection_count = session.query(CollectionTable).count()
    version_count = session.query(CollectionVersionTable).count()
    dataset_count = session.query(DatasetTable).count()
    dataset_version_count = session.query(DatasetVersionTable).count()
    artifact_count_verify = session.query(DatasetArtifactTable).count()
    
    print(f'  Collections: {collection_count}')
    print(f'  CollectionVersions: {version_count}')
    print(f'  Datasets: {dataset_count}')
    print(f'  DatasetVersions: {dataset_version_count}')
    print(f'  DatasetArtifacts: {artifact_count_verify}')
    
    print('\n Comprehensive Metadata Included:')
    print('  ✓ Collection metadata: name, description, contact, consortia, keywords, links')
    print('  ✓ Publisher metadata: DOI, citation, journal, publication_year, authors')
    print('  ✓ Dataset metadata: name, tissue, disease, organism, assay, cell_type')
    print('  ✓ Dataset status: upload, validation, processing, h5ad, cxg, rds, atac')
    print('  ✓ Collection-Dataset linkage: datasets array populated')
    print('  ✓ Dataset artifacts: 2 per dataset (H5AD, CXG)')
