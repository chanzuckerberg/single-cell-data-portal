/* eslint-disable sonarjs/no-duplicate-string */

/* Model of response from /datasets/index API endpoint. */
import { DatasetResponse } from "src/common/queries/filter";

const datasetsIndex = [
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "fcb51725-66cf-458e-b194-b234f31073ba",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d5de5327-4c5c-4e3b-ba62-b10feccf476d",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "Smart-seq",
        ontology_term_id: "EFO:0008930",
      },
    ],
    cell_count: 6171,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "b0d4f376-5345-4815-96d6-c964efde7f91",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: SMARTer_nuclei_MOp",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    collection_id: "e7c47289-22fe-440e-ba24-a4434ac3f379",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "b8306f9d-357a-4cc1-8e15-e109b6241a68",
    name: "Massively multiplex chemical transcriptomics at single-cell resolution - A549",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "pulmonary alveolus epithelium (cell culture)",
        ontology_term_id: "UBERON:0004821 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "49ebd634-c566-485d-ae7c-80ecbb7aa12e",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "e9018b0b-19bc-4f54-be8c-6ad2b67b3f39",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "c85d531e-78bd-4b57-b086-9ce93113135a",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "cc6aecfd-a7ca-467b-8623-237eac8cd3b3",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "43208b33-a29c-4cec-a8dc-a2de5ea10dc0",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "22c5b222-01e6-43c3-b030-2648ae75052f",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "e9c3b028-ebd0-40f2-ac49-eb10254d9f74",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "625654a0-048c-4d15-960b-7ef8f0d19fd9",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "ffdd5b73-6387-4676-803f-5484bb2221b9",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "8f98098d-5488-4ced-8fe8-49fc25f9d0bf",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "59ce12a2-ac1c-45bc-93c2-1ffe1e1bd541",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "bce65912-4a99-48b4-9245-e1f8a0dc626b",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "6fb1a1b9-2566-4a60-8754-bd0191812171",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d18cc36c-c5a2-46fb-b678-43bf4966d9bf",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "7d46cf94-2605-489e-91b4-018cd1324a3a",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "119c926d-3f36-452e-8fda-4cdfb416cb5a",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    collection_id: "e7c47289-22fe-440e-ba24-a4434ac3f379",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "bb7a543d-59df-4a61-9a9f-709fc96365c5",
    name: "Massively multiplex chemical transcriptomics at single-cell resolution - K562",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "bone marrow (cell culture)",
        ontology_term_id: "UBERON:0002371 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    cell_count: 146752,
    collection_id: "009275f5-58dc-4531-9f7c-c682fc13bc73",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "c20979ba-3a5b-48b0-9a22-d6452959383a",
    name: "Massively multiplex chemical transcriptomics at single-cell resolution - K562",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "bone marrow (cell culture)",
        ontology_term_id: "UBERON:0002371 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "ae6e9d0e-7227-4f41-99a3-5cc7b41318bb",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "ae2076fa-4b42-4e79-a0ee-520d45655ced",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "63cdd1d4-93f1-42c2-81dd-3d9c8952bf6b",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "fb907a18-e7e1-4e4b-89ec-8e411fcbaa0c",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 4355,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "43208b33-a29c-4cec-a8dc-a2de5ea10dc0",
    development_stage: [
      {
        label: "human adult stage",
        ontology_term_id: "HsapDv:0000087",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Finnish",
        ontology_term_id: "HANCESTRO:0321",
      },
    ],
    id: "d071ce77-f0c3-41e0-aaa8-e549bd1615f3",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - human",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
    ],
    tissue: [
      {
        label: "cervical lymph node",
        ontology_term_id: "UBERON:0002429",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "4bd543b0-219b-42ac-88e1-72cf41c76791",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "276f87b0-5ed9-42bb-b751-d6cf85500a99",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 24213,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "0966df2d-edda-455a-a451-283f65a538c3",
    name: "Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 3-species integration excitory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "a8040c02-551b-444b-a289-0f10375799b6",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "6550ccf9-b1da-4a3e-9286-32578c11472f",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "Smart-seq",
        ontology_term_id: "EFO:0008930",
      },
    ],
    cell_count: 6288,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "10850a40-eebf-421a-b14e-0d5a5c824f17",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: SMARTer_cells_MOp",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "1a5951e2-357a-456c-bd5b-c44f1ae187c3",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "7a2c1961-43c9-4ce8-88d3-a31fc6254c83",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "2a3f5bcd-0b15-4141-9580-667ae2c20d4c",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "104a5888-49bd-4a22-bb05-67d6284bdc59",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "3f9f55c7-fa28-49e0-bc66-cee6965792f4",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "a180f9d2-ddf0-4957-af0e-2ed6af338b6e",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "8338440a-b274-4af1-914e-30fc6832045d",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "72737b57-86fc-4935-98fe-aca9bc401141",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "54c0cc5b-abe2-4b61-8bc6-54f26d1d30ed",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "ed89fbc5-12f4-4ee7-9eb4-0486715abeb5",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "134ab6a6-055e-4292-b2d6-f34ead9c5caf",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: SFG microglia",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "superior frontal gyrus",
        ontology_term_id: "UBERON:0002661",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 29486,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "1cec9ca7-e7a5-42a7-a428-c98da9c4686c",
    name: "Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 3-species integration inhibitory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "2eaec02d-5930-4057-86ba-099de6c34a29",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "4c85be72-6be5-4307-8b2c-727aaf7508f2",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "317ef624-3029-4e96-93b4-3a1415d0681f",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "35a831e3-ba46-45e3-bafc-0ae2fc69b5c9",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "e22fb133-9d12-45dd-9f1c-f0a8cba9ac9d",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "0427ad33-b934-448f-976f-b3ac7835a3a8",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "771fed72-f2a2-46dc-89b1-d5602eaa08e7",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "62018f0a-2786-48e7-adf3-7df543f90774",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "fe14751e-cd63-4dcd-8f70-24412e48c197",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "3b80cdb6-9aa4-42a5-bab9-8a2ea7f41b9e",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "a1292b98-2c77-40cb-99b9-d8db14613271",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "f31a78c9-0815-4dc6-9c53-9ad0b19d72cc",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "4cfeb880-4e60-41ce-84c6-767cc1d8d306",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "4e08a021-f780-467a-8e2a-0a9b93c88ef3",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "9086a1e5-2df8-4459-b6b8-a31cfafa0041",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d0cc4256-6f06-42e9-b66a-2a88094d67d8",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "6e8998fa-15a0-4ed6-b8c1-728882c358b5",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "bf6866f8-86ae-43f9-9a47-8ac4b8e0faa6",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "92eb5814-dd52-4cde-a3ac-5ee8d2e3ce68",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "1b9bac0b-2ddb-459b-8ece-83ef89156cef",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "1182fbbd-042d-4d0a-ae3c-395b32fc2470",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "ced13899-a5ec-4474-bb71-3fc2f9d8a027",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "c6cc0598-a50a-4fbb-9d26-61652f954cdd",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "0708e1e6-aeff-4fe3-b3c4-d2d14ec1ed53",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "4bd543b0-219b-42ac-88e1-72cf41c76791",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "b7105fe9-df36-436d-8bce-c19d23004e8b",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "aaf8aa07-2e3a-476b-b980-09f372464f2d",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "84ea63f7-26f8-4dc0-ba37-9acfe72c49ee",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 10739,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "1e77e1f0-9970-4631-9150-b6d39c93c88b",
    name: "Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 3-species integration non-nuerons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "b97c4133-5431-4f21-abfe-583cd052d2a2",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "46a0075f-786f-42e4-902f-01a8ff47002c",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 29050,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "1e9700e2-68ce-4f4f-8758-24e1ebb372aa",
    name: "Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 4-species integration excitory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "Visium Spatial Gene Expression",
        ontology_term_id: "EFO:0010961",
      },
    ],
    cell_count: 344,
    cell_type: [
      {
        label: "enteric smooth muscle cell",
        ontology_term_id: "CL:0002504",
      },
      {
        label: "enterocyte",
        ontology_term_id: "CL:0000584",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
      {
        label: "group 3 innate lymphoid cell",
        ontology_term_id: "CL:0001071",
      },
      {
        label: "gut absorptive cell",
        ontology_term_id: "CL:0000677",
      },
      {
        label: "gut endothelial cell",
        ontology_term_id: "CL:0000131",
      },
      {
        label: "interneuron",
        ontology_term_id: "CL:0000099",
      },
      {
        label: "interstitial cell of Cajal",
        ontology_term_id: "CL:0002088",
      },
      {
        label: "intestinal crypt stem cell of large intestine",
        ontology_term_id: "CL:0009016",
      },
      {
        label: "intestinal crypt stem cell of small intestine",
        ontology_term_id: "CL:0009017",
      },
      {
        label: "intestine goblet cell",
        ontology_term_id: "CL:0019031",
      },
      {
        label: "leukocyte",
        ontology_term_id: "CL:0000738",
      },
      {
        label: "monocyte",
        ontology_term_id: "CL:0000576",
      },
      {
        label: "motor neuron",
        ontology_term_id: "CL:0000100",
      },
      {
        label: "myofibroblast cell",
        ontology_term_id: "CL:0000186",
      },
      {
        label: "neuroendocrine cell",
        ontology_term_id: "CL:0000165",
      },
      {
        label: "pericyte cell",
        ontology_term_id: "CL:0000669",
      },
      {
        label: "progenitor cell",
        ontology_term_id: "CL:0011026",
      },
      {
        label: "secretory cell",
        ontology_term_id: "CL:0000151",
      },
    ],
    collection_id: "a7fb73ef-78e5-4014-bbab-10e0c49f2dd7",
    development_stage: [
      {
        label: "fetal stage",
        ontology_term_id: "HsapDv:0000037",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "66ba589f-3392-44b1-a20d-02995636d189",
    is_primary_data: "PRIMARY",
    name: "Spatiotemporal analysis of human intestinal development at single-cell resolution: Fetal A7",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "intestine",
        ontology_term_id: "UBERON:0000160",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 14783,
    cell_type: [
      {
        label: "central memory CD4-positive, alpha-beta T cell",
        ontology_term_id: "CL:0000904",
      },
      {
        label: "classical monocyte",
        ontology_term_id: "CL:0000860",
      },
      {
        label: "conventional dendritic cell",
        ontology_term_id: "CL:0000990",
      },
      {
        label: "effector memory CD8-positive, alpha-beta T cell",
        ontology_term_id: "CL:0000913",
      },
      {
        label: "hematopoietic stem cell",
        ontology_term_id: "CL:0000037",
      },
      {
        label: "immature neutrophil",
        ontology_term_id: "CL:0000776",
      },
      {
        label: "lymphocyte",
        ontology_term_id: "CL:0000542",
      },
      {
        label: "memory B cell",
        ontology_term_id: "CL:0000787",
      },
      {
        label: "naive B cell",
        ontology_term_id: "CL:0000788",
      },
      {
        label: "naive thymus-derived CD4-positive, alpha-beta T cell",
        ontology_term_id: "CL:0000895",
      },
      {
        label: "naive thymus-derived CD8-positive, alpha-beta T cell",
        ontology_term_id: "CL:0000900",
      },
      {
        label: "natural killer cell",
        ontology_term_id: "CL:0000623",
      },
      {
        label: "neutrophil",
        ontology_term_id: "CL:0000775",
      },
      {
        label: "non-classical monocyte",
        ontology_term_id: "CL:0000875",
      },
      {
        label: "plasmablast",
        ontology_term_id: "CL:0000980",
      },
      {
        label: "plasmacytoid dendritic cell, human",
        ontology_term_id: "CL:0001058",
      },
      {
        label: "platelet",
        ontology_term_id: "CL:0000233",
      },
      {
        label: "transitional stage B cell",
        ontology_term_id: "CL:0000818",
      },
    ],
    collection_id: "a7fb73ef-78e5-4014-bbab-10e0c49f2dd7",
    development_stage: [
      {
        label: "39-year-old human stage",
        ontology_term_id: "HsapDv:0000133",
      },
      {
        label: "78-year-old human stage",
        ontology_term_id: "HsapDv:0000172",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "cd024f39-2bae-45fc-951b-941d390ca7bf",
    is_primary_data: "PRIMARY",
    name: "Individual Single-Cell RNA-seq PBMC Data from Guo et al.",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
    ],
    tissue: [
      {
        label: "blood",
        ontology_term_id: "UBERON:0000178",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "Smart-seq",
        ontology_term_id: "EFO:0008930",
      },
    ],
    cell_count: 406187,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "317c52fa-5dbc-44cb-86b2-8e34c400c18f",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 71183,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "37e6e58a-cbc8-4599-bf90-69fe2ef7eaff",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: 10X_cells_v3_AIBS",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 40166,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "43fed16e-460a-44ab-8c3d-dd591dd7ab62",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: 10X_nuclei_v3_AIBS",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 611034,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "adult",
        ontology_term_id: "EFO:0001272",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "af46ed05-026f-48f3-9327-d763ccd617b4",
    name: "A transcriptomic atlas of the mouse cerebellum",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "cerebellum",
        ontology_term_id: "UBERON:0002037",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "facc3764-25de-4039-af61-cf14778e00b1",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "cb5d3451-5f71-47bd-aa4c-f97feda2172b",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "7acf8645-e54e-4d40-9c51-d3a2c5f055c2",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "6932ebc5-621c-471d-954f-4d0051de05c2",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "5d8784c3-f8f4-443b-81a5-328a6f96f8a9",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d1807757-5932-46ba-aa6a-503eca8f59f6",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "95981e6b-ce16-486b-ac32-b550d39f8d40",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "50c355d9-4233-44e4-ab41-f5f65cef9dbd",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "400cc5bc-b96d-40d4-ae78-4a1a186ed62e",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "fb2969de-a2ca-4be7-bf58-7bb0d8115460",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "60eb7771-cad1-471c-a0fd-abefcdf707df",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "e33a8d67-aba6-4f2c-9887-800f0a352bcc",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "4ea95e8d-ad26-40ab-b37d-a188034549ac",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d5b69c9f-9a6d-4318-bcf4-d8fcecd45fe9",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "21ec65c6-7148-4360-b79d-8756398753d1",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "321b0b64-c918-4a6a-a5b9-f6d7268d12cc",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "cb86e616-80ec-4fc4-8c66-b8e458be5b02",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "5f85602a-f4a2-466d-bc08-ac5fcfbae1a3",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "335089bf-4aee-41cf-abf0-e57c61d75518",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "6689c60b-bff7-40c9-9ce9-365125afa8c0",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "1442277d-d3a2-4ad7-aa7e-855ec9bf2648",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "bf33daf6-2dcb-49b7-9d65-4cdb53442d2c",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "4c093fc9-d8c6-4240-ba77-45dfcb6f202a",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "e9c19207-ec4a-40a1-aa2c-d2c81d14777a",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "32c35b29-f65a-494a-b417-dd88a93fa687",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "c9b8c5db-7dac-4046-9d5d-66af0f0c9071",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "cd81c7d4-825c-4c0d-8293-3c353d176f9a",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "0f95707b-99ca-47bd-8eec-3ce6e2357757",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "568199d1-24bf-42b9-b6ef-82117e42c5f3",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "12c56474-90f2-4234-94a1-9ccbccd306bf",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "47c5b7c2-e786-40b5-888b-ed5bdee7a1c8",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: superior frontal gyrus",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "superior frontal gyrus",
        ontology_term_id: "UBERON:0002661",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "5f788166-c98c-4705-aed9-d0a229b33a2c",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "9d7f09a5-41d7-460b-8112-f37902646ae8",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "4de33b39-a464-4596-8552-29d52dd3de26",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "33e5da1f-8c80-49d3-b71b-c78be14f68c5",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "f913a951-b1d8-40ba-9e20-916b2e7dd391",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "9e066534-22dd-4ffd-8879-6b09a274deaf",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "d5a41cf2-6c65-4a6d-affb-8d8dd86c2fe5",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "018a8dd0-6160-4446-b39e-92aa02d7e1c7",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "47db8f76-3a1c-4acf-9ba2-32a57c240a70",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: EC microglia",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    cell_count: 143015,
    collection_id: "f6eb8c76-1e3d-4ee8-bd1e-8f645919d9bc",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "4af94bb8-eb36-4d62-b2af-6dc32e65c03a",
    name: "test dataset for upload",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "pulmonary alveolus epithelium (cell culture)",
        ontology_term_id: "UBERON:0004821 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 2638,
    collection_id: "9303f3d8-3199-4421-b45b-8581785a4237",
    development_stage: [
      {
        label: "neurula stage",
        ontology_term_id: "HsapDv:0000012",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "503a8479-5a4b-45fb-a8d1-d5b2350141e3",
    name: "Test PBMC 3k dataset",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "mouthpart",
        ontology_term_id: "UBERON:6004520",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "e66b3812-cdeb-4a20-b343-d80377fd2c35",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "5c8bb3a7-9961-471f-85a6-b5c4ff5977db",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "4f57a163-467e-45d6-ab12-6b44467d26bc",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "7e98a1f7-4b7c-4b6e-80a4-ed3d50df7863",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "3147526a-d574-40a9-986b-335fe06e7a3b",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "6980bec5-e552-488b-ac8c-db20af07711f",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "6feb0aae-88ef-487a-ade0-952b10704b24",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "125aa8c6-2877-417b-937f-6fd0ca1cca45",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "8318eb93-7963-4a05-9f0f-01de6ae67d20",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "200bd77e-7b72-4286-8cad-5ffeac3a9a8f",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "899550af-0918-43e8-bfd9-e7cdd3a6e090",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "e7055be4-35d9-4136-9943-72047f86f19a",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "8be9ea96-68cb-407a-8466-1deefedb19e4",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "5f0986bb-025f-4d1f-8a2a-2f05fa9606f5",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "2d1a245b-496b-42d2-a6e1-ea123e762c26",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "7f5f580b-02e7-41b2-9cec-08b330b09100",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "Visium Spatial Gene Expression",
        ontology_term_id: "EFO:0010961",
      },
    ],
    cell_count: 709,
    cell_type: [
      {
        label: "enteric smooth muscle cell",
        ontology_term_id: "CL:0002504",
      },
      {
        label: "enterocyte",
        ontology_term_id: "CL:0000584",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
      {
        label: "glial cell",
        ontology_term_id: "CL:0000125",
      },
      {
        label: "gut absorptive cell",
        ontology_term_id: "CL:0000677",
      },
      {
        label: "gut endothelial cell",
        ontology_term_id: "CL:0000131",
      },
      {
        label: "inhibitory motor neuron",
        ontology_term_id: "CL:0008015",
      },
      {
        label: "interstitial cell of Cajal",
        ontology_term_id: "CL:0002088",
      },
      {
        label: "intestinal crypt stem cell of large intestine",
        ontology_term_id: "CL:0009016",
      },
      {
        label: "intestinal crypt stem cell of small intestine",
        ontology_term_id: "CL:0009017",
      },
      {
        label: "intestine goblet cell",
        ontology_term_id: "CL:0019031",
      },
      {
        label: "leukocyte",
        ontology_term_id: "CL:0000738",
      },
      {
        label: "mesothelial cell",
        ontology_term_id: "CL:0000077",
      },
      {
        label: "motor neuron",
        ontology_term_id: "CL:0000100",
      },
      {
        label: "myofibroblast cell",
        ontology_term_id: "CL:0000186",
      },
      {
        label: "neural cell",
        ontology_term_id: "CL:0002319",
      },
      {
        label: "pericyte cell",
        ontology_term_id: "CL:0000669",
      },
      {
        label: "progenitor cell",
        ontology_term_id: "CL:0011026",
      },
      {
        label: "secretory cell",
        ontology_term_id: "CL:0000151",
      },
      {
        label: "transit amplifying cell of colon",
        ontology_term_id: "CL:0009011",
      },
    ],
    collection_id: "8318eb93-7963-4a05-9f0f-01de6ae67d20",
    development_stage: [
      {
        label: "fetal stage",
        ontology_term_id: "HsapDv:0000037",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d22cb620-9310-4e7f-929d-d5af10e584ed",
    is_primary_data: "PRIMARY",
    name: "Spatiotemporal analysis of human intestinal development at single-cell resolution: Fetal A8",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "intestine",
        ontology_term_id: "UBERON:0000160",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "b999c227-f71d-4832-879c-ac24dccf09d3",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "edfe043c-20b1-43fe-a778-7297ed303fe5",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "805fe349-4bbd-44c6-980a-2dca03c95ff2",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "46ecd9e2-f716-4efb-a13a-3ad5533c99d1",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "Smart-seq2 protocol",
        ontology_term_id: "EFO:0008442",
      },
    ],
    collection_id: "18e6337b-8c36-4977-9ac8-ca8aa5494481",
    development_stage: [
      {
        label: "75-year-old human stage",
        ontology_term_id: "HsapDv:0000169",
      },
      {
        label: "46-year-old human stage",
        ontology_term_id: "HsapDv:0000140",
      },
      {
        label: "51-year-old human stage",
        ontology_term_id: "HsapDv:0000145",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "53ef8427-ad49-4c4c-8cf4-1d78e8611da4",
    name: "Krasnow Lab Human Lung Cell Atlas, Smart-seq2",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "lung",
        ontology_term_id: "UBERON:0002048",
      },
      {
        label: "blood",
        ontology_term_id: "UBERON:0000178",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "55d11a7b-55d4-4224-be0f-0d216346a38d",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "a89a2f8c-f231-4625-bc1e-4d24cdf90c39",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "1276094c-9b71-482e-95b6-82d2a25ac16b",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "a14bc3c5-e0a6-40a6-b35d-729ff6ffc48a",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "546b3999-6950-4da7-bbe5-40bf5408340b",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: EC excitatory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "f5f9366c-e240-44c8-9eb3-537591a46a8c",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "c0759980-bbbb-4113-b2ea-53056b47658d",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "9d0fb12b-879f-4897-bb9e-d67b474780fb",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "815d6a90-fce3-4648-abfd-11fd4b709ccd",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "b9ed3583-ef3b-488f-acbf-a5ebda3720d9",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "1ef6f45e-3fd9-4ba9-8c4c-0f185d6e8cdd",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "Drop-seq",
        ontology_term_id: "EFO:0008722",
      },
    ],
    collection_id: "673637cf-dcb7-45e1-bb88-72a27c50c8ca",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "severe acute respiratory syndrome",
        ontology_term_id: "MONDO:0005091",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "5ae1469a-f74e-4141-ae89-c6004482312d",
    name: "Single-cell gene expression profiling of SARS-CoV-2 infected human cell lines - Calu-3",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "lung epithelium (cell culture)",
        ontology_term_id: "UBERON:0000115 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    cell_count: 143015,
    collection_id: "6eaeca39-c187-4f10-86af-50572ff6c333",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "5b484cd3-d1f9-456e-9502-a6bddc92caef",
    name: "",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "pulmonary alveolus epithelium (cell culture)",
        ontology_term_id: "UBERON:0004821 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "5dbb5dc2-c174-489f-9106-0439c9a4a8b7",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: EC oligodendrocyte",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "28a014b5-8aed-483d-88cf-b97c0f31e01a",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "3575431c-ca5d-415e-9be9-30f10d9d7dcf",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 122641,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "64491279-920b-46a8-8899-6b1e01a4c25a",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: 10X_cells_v2_AIBS",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "0514143e-0746-4cdd-9d1e-60f7a54d8b06",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "da9d7292-a65d-4a12-9be6-fad7ccde515b",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "68003376-d1a8-4b1c-9ce9-5ce0080f6e6b",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: SFG inhibitory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "superior frontal gyrus",
        ontology_term_id: "UBERON:0002661",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x Genomics",
        ontology_term_id: "EFO:0008995",
      },
      {
        label: "CITE-Seq",
        ontology_term_id: "EFO:0009294",
      },
      {
        label: "Seq-Well",
        ontology_term_id: "EFO:0008919",
      },
    ],
    cell_count: 239696,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "human adult stage",
        ontology_term_id: "HsapDv:0000087",
      },
      {
        label: "human early adulthood stage",
        ontology_term_id: "HsapDv:0000088",
      },
      {
        label: "human late adulthood stage",
        ontology_term_id: "HsapDv:0000091",
      },
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "6ae3b815-e7a4-430a-b935-3d952ce3cb99",
    name: "COVID-19 Immune Altas: Integration of 5 public COVID-19 PBMC single-cell datasets",
    organism: [
      {
        label: "homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "blood",
        ontology_term_id: "UBERON:0000178",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    cell_count: 146752,
    collection_id: "009275f5-58dc-4531-9f7c-c682fc13bc73",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "6d5eadb8-8359-4fa6-a092-a001e1e92f7d",
    name: "Massively multiplex chemical transcriptomics at single-cell resolution - K562",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "bone marrow (cell culture)",
        ontology_term_id: "UBERON:0002371 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "microwell-seq",
        ontology_term_id: "",
      },
    ],
    collection_id: "4bd1b31f-d5ea-41b4-8003-0cb61c386ef8",
    development_stage: [
      {
        label: "embryonic human stage",
        ontology_term_id: "HsapDv:0000002",
      },
      {
        label: "fetal stage",
        ontology_term_id: "HsapDv:0000037",
      },
      {
        label: "human adult stage",
        ontology_term_id: "HsapDv:0000087",
      },
      {
        label: "newborn human stage",
        ontology_term_id: "HsapDv:0000082",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Han Chinese",
        ontology_term_id: "HANCESTRO:0027",
      },
    ],
    id: "6d76b779-f811-4606-9380-4140506d7f2d",
    name: "Construction of a human cell landscape at single-cell level",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "bone marrow",
        ontology_term_id: "UBERON:0002371",
      },
      {
        label: "cerebellum",
        ontology_term_id: "UBERON:0002037",
      },
      {
        label: "lung",
        ontology_term_id: "UBERON:0002048",
      },
      {
        label: "uterine cervix",
        ontology_term_id: "UBERON:0000002",
      },
      {
        label: "blood",
        ontology_term_id: "UBERON:0000178",
      },
      {
        label: "adrenal gland",
        ontology_term_id: "UBERON:0002369",
      },
      {
        label: "eye",
        ontology_term_id: "UBERON:0000970",
      },
      {
        label: "vault of skull",
        ontology_term_id: "UBERON:0004339",
      },
      {
        label: "omentum",
        ontology_term_id: "UBERON:0003688",
      },
      {
        label: "transverse colon",
        ontology_term_id: "UBERON:0001157",
      },
      {
        label: "artery",
        ontology_term_id: "UBERON:0001637",
      },
      {
        label: "duodenum",
        ontology_term_id: "UBERON:0002114",
      },
      {
        label: "intestine",
        ontology_term_id: "UBERON:0000160",
      },
      {
        label: "trachea",
        ontology_term_id: "UBERON:0003126",
      },
      {
        label: "stomach",
        ontology_term_id: "UBERON:0000945",
      },
      {
        label: "vermiform appendix",
        ontology_term_id: "UBERON:0001154",
      },
      {
        label: "chorionic villus",
        ontology_term_id: "UBERON:0007106",
      },
      {
        label: "HESC",
        ontology_term_id: "",
      },
      {
        label: "umbilical cord blood",
        ontology_term_id: "UBERON:0012168",
      },
      {
        label: "prostate gland",
        ontology_term_id: "UBERON:0002367",
      },
      {
        label: "skin of body",
        ontology_term_id: "UBERON:0002097",
      },
      {
        label: "thymus",
        ontology_term_id: "UBERON:0002370",
      },
      {
        label: "spinal cord",
        ontology_term_id: "UBERON:0002240",
      },
      {
        label: "rectum",
        ontology_term_id: "UBERON:0001052",
      },
      {
        label: "ileum",
        ontology_term_id: "UBERON:0002116",
      },
      {
        label: "testis",
        ontology_term_id: "UBERON:0000473",
      },
      {
        label: "gall bladder",
        ontology_term_id: "UBERON:0002110",
      },
      {
        label: "ureter",
        ontology_term_id: "UBERON:0000056",
      },
      {
        label: "ovary",
        ontology_term_id: "UBERON:0000992",
      },
      {
        label: "jejunum",
        ontology_term_id: "UBERON:0002115",
      },
      {
        label: "sigmoid colon",
        ontology_term_id: "UBERON:0001159",
      },
      {
        label: "fallopian tube",
        ontology_term_id: "UBERON:0003889",
      },
      {
        label: "rib",
        ontology_term_id: "UBERON:0002228",
      },
      {
        label: "ascending colon",
        ontology_term_id: "UBERON:0001156",
      },
      {
        label: "placenta",
        ontology_term_id: "UBERON:0001987",
      },
      {
        label: "uterus",
        ontology_term_id: "UBERON:0000995",
      },
      {
        label: "esophagus",
        ontology_term_id: "UBERON:0001043",
      },
      {
        label: "brain",
        ontology_term_id: "UBERON:0000955",
      },
      {
        label: "spleen",
        ontology_term_id: "UBERON:0002106",
      },
      {
        label: "pancreas",
        ontology_term_id: "UBERON:0001264",
      },
      {
        label: "pleura",
        ontology_term_id: "UBERON:0000977",
      },
      {
        label: "thyroid gland",
        ontology_term_id: "UBERON:0002046",
      },
      {
        label: "adipose tissue",
        ontology_term_id: "UBERON:0001013",
      },
      {
        label: "temporal lobe",
        ontology_term_id: "UBERON:0001871",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
      {
        label: "heart",
        ontology_term_id: "UBERON:0000948",
      },
      {
        label: "muscle organ",
        ontology_term_id: "UBERON:0001630",
      },
      {
        label: "bladder organ",
        ontology_term_id: "UBERON:0018707",
      },
      {
        label: "liver",
        ontology_term_id: "UBERON:0002107",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "3477c0b0-cfa1-4fa6-8961-f55f74a48a5b",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "46f1a16f-7d4d-4dee-b885-337b6e0f6347",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "18e6337b-8c36-4977-9ac8-ca8aa5494481",
    development_stage: [
      {
        label: "75-year-old human stage",
        ontology_term_id: "HsapDv:0000169",
      },
      {
        label: "46-year-old human stage",
        ontology_term_id: "HsapDv:0000140",
      },
      {
        label: "51-year-old human stage",
        ontology_term_id: "HsapDv:0000145",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "7307297e-0b8c-4894-8f41-ead34ba90fa5",
    name: "Krasnow Lab Human Lung Cell Atlas, 10X",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "lung",
        ontology_term_id: "UBERON:0002048",
      },
      {
        label: "blood",
        ontology_term_id: "UBERON:0000178",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "7e31ce20-3c77-4e1e-af46-f168177e2855",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: EC astrocytes",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 159738,
    collection_id: "cd48bc53-7021-4a75-9685-54b374e825de",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "7e742c78-038c-41bb-952c-e0645dc16164",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: 10X_nuclei_v3_Broad",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 2638,
    collection_id: "6015e80d-0203-453f-a033-bb59bd439e1d",
    development_stage: [
      {
        label: "neurula stage",
        ontology_term_id: "HsapDv:0000012",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "836e7bf4-2e11-48bc-b773-1f6b65cdfa27",
    name: "Test PBMC 3k dataset",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "mouthpart",
        ontology_term_id: "UBERON:6004520",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "8786a595-ac9a-4e5f-b1a2-1717596dc333",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: EC inhibitory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "d9265653-0ff7-40b1-a738-130f1a73a990",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "1f3922ed-a65d-4966-9c90-eed2aec39091",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "26c14cd9-e85e-4f29-b633-0810e6f358f2",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "ce85512c-469b-4c44-84b1-a56f678d7831",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "7acf8645-e54e-4d40-9c51-d3a2c5f055c2",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "1b3986c6-b3da-482b-b3a7-6a77d514e345",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "d3df11ad-7af3-4662-8bf1-8aa6a841e026",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "29bf810d-0a96-4d55-9f2d-93c449b9bfc2",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "e92c1462-dfe6-4616-8a20-8879c6c80c6b",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "f0e1289f-5ce9-41df-a4e1-a5ec81cf5ad7",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "8df3280a-6f90-4b86-ae69-babd0cf6297a",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: SFG astrocytes",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "superior frontal gyrus",
        ontology_term_id: "UBERON:0002661",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    collection_id: "a55cac37-629f-4d28-9b3f-4ece7a7ac174",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "8e3bf114-5491-4855-99a1-83770a3f75fa",
    name: "Single-cell longitudinal analysis of SARS-CoV-2 infection in human bronchial epithelial cells",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "epithelium of bronchus (cell culture)",
        ontology_term_id: "UBERON:0002031 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "Smart-seq",
        ontology_term_id: "EFO:0008930",
      },
    ],
    cell_count: 6288,
    collection_id: "cd4c7ed9-1b0e-4d65-8e77-a0d2e15b1623",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "91e5af67-5f91-4ed7-a20c-b4454b3c25e5",
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: SMARTer_cells_MOp",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 2638,
    collection_id: "4826a9e2-68b4-4017-ad91-ba56605c9396",
    development_stage: [
      {
        label: "neurula stage",
        ontology_term_id: "HsapDv:0000012",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "93e22f02-5017-4e79-9b4c-d4600f9de930",
    name: "Test PBMC 3k dataset",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "mouthpart",
        ontology_term_id: "UBERON:6004520",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5500,
    collection_id: "009275f5-58dc-4531-9f7c-c682fc13bc73",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "949b764f-d334-40b8-a317-d10715718574",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer\u2019s Disease: EC astrocytes",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "e92c1462-dfe6-4616-8a20-8879c6c80c6b",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "906b206b-7a3f-4c91-b960-f5654392bfa1",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "97e7185f-815c-4b4b-a53d-55c73d215a0a",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: SFG oligodendrocyte",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "superior frontal gyrus",
        ontology_term_id: "UBERON:0002661",
      },
    ],
  },
  {
    assay: [
      {
        label: "Drop-seq",
        ontology_term_id: "EFO:0008722",
      },
    ],
    collection_id: "673637cf-dcb7-45e1-bb88-72a27c50c8ca",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "severe acute respiratory syndrome",
        ontology_term_id: "MONDO:0005091",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "a330bb4f-7495-4773-939c-9a9e600dd96b",
    name: "Single-cell gene expression profiling of SARS-CoV-2 infected human cell lines - H1299",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "lung epithelium (cell culture)",
        ontology_term_id: "UBERON:0000115 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "bbc86167-c02c-478a-9881-2dc91f46a741",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "c7d504c5-59c0-4eb1-b30a-41fccc7a7f06",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "c6b2a804-e24d-4b0a-be95-d615ff872983",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: SFG excitatory neurons",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "superior frontal gyrus",
        ontology_term_id: "UBERON:0002661",
      },
    ],
  },
  {
    assay: [
      {
        label: "Seq-Well",
        ontology_term_id: "EFO:0008919",
      },
    ],
    collection_id: "94e005bf-1d80-4938-a442-6226f0141f74",
    development_stage: [
      {
        label: "human adult stage",
        ontology_term_id: "HsapDv:0000087",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "ca9e42d9-de52-4e7b-9379-352b27f46411",
    name: "Single-cell atlas of peripheral immune response to SARS-CoV-2 infection",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
      {
        label: "female",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "blood",
        ontology_term_id: "UBERON:0000178",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    collection_id: "e7c47289-22fe-440e-ba24-a4434ac3f379",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "d485ea79-7e05-4406-874a-a854892f8722",
    name: "Massively multiplex chemical transcriptomics at single-cell resolution - MCF7",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "mammary gland (cell culture)",
        ontology_term_id: "UBERON:0001911 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    collection_id: "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "Alzheimer disease",
        ontology_term_id: "MONDO:0004975",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "daad0b15-b096-4c4e-a73d-b2d629e7c28a",
    name: "Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease: caudal entorhinal cortex",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "male",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "entorhinal cortex",
        ontology_term_id: "UBERON:0002728",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    cell_count: 146752,
    collection_id: "009275f5-58dc-4531-9f7c-c682fc13bc73",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "effac18f-51b8-4441-a761-016e3aff7c46",
    name: "Massively multiplex chemical transcriptomics at single-cell resolution - K562",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "bone marrow (cell culture)",
        ontology_term_id: "UBERON:0002371 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-plex",
        ontology_term_id: "",
      },
    ],
    cell_count: 143015,
    collection_id: "6eaeca39-c187-4f10-86af-50572ff6c333",
    development_stage: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "",
      },
    ],
    id: "ffa76ea4-ce17-4759-a45d-21170dc708b8",
    name: "",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "unknown",
        sex_ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "pulmonary alveolus epithelium (cell culture)",
        ontology_term_id: "UBERON:0004821 (cell culture)",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "a5ea6c99-f75a-4299-a04e-fe9eb174c4ec",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "ec89a163-77bf-4e08-ae41-335b457f5eeb",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "8d500e11-a671-4e01-a3d2-6ee9ceec014e",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "8ba756fe-0432-4505-8b85-d499548c4339",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "a5ea6c99-f75a-4299-a04e-fe9eb174c4ec",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "80b61fa6-56a9-42b5-9903-7b8108349f36",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "54d83db8-363d-4161-adc6-7b43ca8c88b5",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "e6b3ffbc-e8b6-4f1f-a974-4abdfdbb47eb",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "a6e48cd2-002a-492f-a0aa-58b9e5c819e2",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "c46bf249-c5f4-480c-a960-e22fa431987a",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "ead713c6-f089-4f36-aaec-1f1ba7ff696e",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "b515c8de-dfe8-43a0-a353-c5d60cb1b0fe",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 5625,
    cell_type: [
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
    ],
    collection_id: "a6e48cd2-002a-492f-a0aa-58b9e5c819e2",
    development_stage: [
      {
        label: "early adult stage",
        ontology_term_id: "MmusDv:0000061",
      },
    ],
    disease: [
      {
        label: "lymphadenitis",
        ontology_term_id: "MONDO:0002052",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "na",
        ontology_term_id: "na",
      },
    ],
    id: "4e40ea4a-4412-4a1c-b760-e252c2b82f5d",
    is_primary_data: "PRIMARY",
    name: "A Single-Cell Transcriptional Roadmap of the Mouse and Human Lymph Node Lymphatic Vasculature - mouse",
    organism: [
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "peripheral lymph node",
        ontology_term_id: "UBERON:0003968",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "75657991-b14b-4ff5-81f7-708dc84c4905",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "de38bf9e-48f1-4aa0-9a79-f8acce6792a5",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "d97692cd-7c05-4c1f-87dd-76a1b6f311d3",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "f14d86ab-c3d6-4ba8-92dd-dcefe18e1dcf",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "cd77ad9f-3f37-4504-b492-790fb3d08164",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "a39947ca-ac3f-4f2c-88bd-e6ac72d638cb",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "8a89918a-54e7-4387-be8f-9ff84b1d2b6c",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "d5fddd48-4c7e-46fa-b08c-236e52d86502",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "911c059e-ac19-4092-b2c1-ac327c44d102",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "9bd417e0-be74-4737-a506-93de8de52276",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "e82ab1bf-be5b-40e5-8f69-b12a9fd186d7",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "c3712d37-b328-4f2d-b319-933158e58f8f",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "42bd7631-c5c7-4673-99ab-bf5b1c07f097",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "23644339-28dc-40e0-bb82-89821abb95fb",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "9a216899-fb8f-4a81-8402-df099c18b81b",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "1f5caefe-76f0-4ac4-a2e5-c7e02c5ad194",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "b5d78e8c-cb50-43cc-93fc-aa98f7a2db3d",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "1a921541-f645-4f38-a345-41dc501c6556",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "sci-RNA-seq",
        ontology_term_id: "EFO:0010550",
      },
    ],
    cell_count: 4062980,
    cell_type: [
      {
        label: "Muller cell",
        ontology_term_id: "CL:0011107",
      },
      {
        label: "Purkinje cell",
        ontology_term_id: "CL:0000121",
      },
      {
        label: "Schwann cell",
        ontology_term_id: "CL:0002573",
      },
      {
        label: "acinar cell",
        ontology_term_id: "CL:0000622",
      },
      {
        label: "amacrine cell",
        ontology_term_id: "CL:0000561",
      },
      {
        label: "astrocyte",
        ontology_term_id: "CL:0000127",
      },
      {
        label: "bipolar neuron",
        ontology_term_id: "CL:0000103",
      },
      {
        label: "blood vessel endothelial cell",
        ontology_term_id: "CL:0000071",
      },
      {
        label: "brush cell",
        ontology_term_id: "CL:0002204",
      },
      {
        label: "cardiac muscle cell",
        ontology_term_id: "CL:0000746",
      },
      {
        label: "cell of skeletal muscle",
        ontology_term_id: "CL:0000188",
      },
      {
        label: "chorionic trophoblast cell",
        ontology_term_id: "CL:0011101",
      },
      {
        label: "chromaffin cell",
        ontology_term_id: "CL:0000166",
      },
      {
        label: "ciliated epithelial cell",
        ontology_term_id: "CL:0000067",
      },
      {
        label: "corneal epithelial cell",
        ontology_term_id: "CL:0000575",
      },
      {
        label: "cortical cell of adrenal gland",
        ontology_term_id: "CL:0002097",
      },
      {
        label: "dermis microvascular lymphatic vessel endothelial cell",
        ontology_term_id: "CL:2000041",
      },
      {
        label: "endocardial cell",
        ontology_term_id: "CL:0002350",
      },
      {
        label: "endocrine cell",
        ontology_term_id: "CL:0000163",
      },
      {
        label: "endothelial cell of lymphatic vessel",
        ontology_term_id: "CL:0002138",
      },
      {
        label: "endothelial cell of vascular tree",
        ontology_term_id: "CL:0002139",
      },
      {
        label: "enteric neuron",
        ontology_term_id: "CL:0007011",
      },
      {
        label: "epicardial adipocyte",
        ontology_term_id: "CL:1000309",
      },
      {
        label: "epithelial cell of lower respiratory tract",
        ontology_term_id: "CL:0002632",
      },
      {
        label: "epithelial cell of thymus",
        ontology_term_id: "CL:0002293",
      },
      {
        label: "erythroblast",
        ontology_term_id: "CL:0000765",
      },
      {
        label: "excitatory neuron",
        ontology_term_id: "CL:0008030",
      },
      {
        label: "extravillous trophoblast",
        ontology_term_id: "CL:0008036",
      },
      {
        label: "ganglion interneuron",
        ontology_term_id: "CL:0000397",
      },
      {
        label: "glial cell",
        ontology_term_id: "CL:0000125",
      },
      {
        label: "goblet cell",
        ontology_term_id: "CL:0000160",
      },
      {
        label: "granule cell",
        ontology_term_id: "CL:0000120",
      },
      {
        label: "hematopoietic cell",
        ontology_term_id: "CL:0000988",
      },
      {
        label: "hematopoietic stem cell",
        ontology_term_id: "CL:0000037",
      },
      {
        label: "hepatic stellate cell",
        ontology_term_id: "CL:0000632",
      },
      {
        label: "hepatoblast",
        ontology_term_id: "CL:0005026",
      },
      {
        label: "inhibitory interneuron",
        ontology_term_id: "CL:0000498",
      },
      {
        label: "inhibitory neuron",
        ontology_term_id: "CL:0008029",
      },
      {
        label: "innate lymphoid cell",
        ontology_term_id: "CL:0001065",
      },
      {
        label: "intestinal epithelial cell",
        ontology_term_id: "CL:0002563",
      },
      {
        label: "lens fiber cell",
        ontology_term_id: "CL:0011004",
      },
      {
        label: "macroglial cell",
        ontology_term_id: "CL:0000126",
      },
      {
        label: "megakaryocyte",
        ontology_term_id: "CL:0000556",
      },
      {
        label: "mesangial cell",
        ontology_term_id: "CL:0000650",
      },
      {
        label: "mesothelial cell",
        ontology_term_id: "CL:0000077",
      },
      {
        label: "microglial cell",
        ontology_term_id: "CL:0000129",
      },
      {
        label: "myeloid cell",
        ontology_term_id: "CL:0000763",
      },
      {
        label: "native cell",
        ontology_term_id: "CL:0000003",
      },
      {
        label: "neuroendocrine cell",
        ontology_term_id: "CL:0000165",
      },
      {
        label: "neuron",
        ontology_term_id: "CL:0000540",
      },
      {
        label: "neuron associated cell (sensu Vertebrata)",
        ontology_term_id: "CL:0000123",
      },
      {
        label: "neuronal brush cell",
        ontology_term_id: "CL:0000555",
      },
      {
        label: "oligodendrocyte",
        ontology_term_id: "CL:0000128",
      },
      {
        label: "pancreatic ductal cell",
        ontology_term_id: "CL:0002079",
      },
      {
        label: "photoreceptor cell",
        ontology_term_id: "CL:0000210",
      },
      {
        label: "professional antigen presenting cell",
        ontology_term_id: "CL:0000145",
      },
      {
        label: "regular atrial cardiac myocyte",
        ontology_term_id: "CL:0002129",
      },
      {
        label: "retina horizontal cell",
        ontology_term_id: "CL:0000745",
      },
      {
        label: "retinal ganglion cell",
        ontology_term_id: "CL:0000740",
      },
      {
        label: "retinal pigment epithelial cell",
        ontology_term_id: "CL:0002586",
      },
      {
        label: "skeletal muscle satellite cell",
        ontology_term_id: "CL:0000594",
      },
      {
        label: "smooth muscle cell",
        ontology_term_id: "CL:0000192",
      },
      {
        label: "squamous epithelial cell",
        ontology_term_id: "CL:0000076",
      },
      {
        label: "stellate neuron",
        ontology_term_id: "CL:0000122",
      },
      {
        label: "stromal cell",
        ontology_term_id: "CL:0000499",
      },
      {
        label: "sympathetic neuron",
        ontology_term_id: "CL:0011103",
      },
      {
        label: "syncytiotrophoblast cell",
        ontology_term_id: "CL:0000525",
      },
      {
        label: "taste receptor cell",
        ontology_term_id: "CL:0000209",
      },
      {
        label: "thymocyte",
        ontology_term_id: "CL:0000893",
      },
      {
        label: "trophoblast giant cell",
        ontology_term_id: "CL:0002488",
      },
      {
        label: "ventricular cardiac muscle cell",
        ontology_term_id: "CL:2000046",
      },
      {
        label: "visceromotor neuron",
        ontology_term_id: "CL:0005025",
      },
    ],
    collection_id: "144c0b61-438c-4dd2-b999-d41a137520ce",
    development_stage: [
      {
        label: "10th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000047",
      },
      {
        label: "12th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000049",
      },
      {
        label: "13th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000050",
      },
      {
        label: "14th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000051",
      },
      {
        label: "15th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000052",
      },
      {
        label: "16th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000053",
      },
      {
        label: "17th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000054",
      },
      {
        label: "18th week post-fertilization human stage",
        ontology_term_id: "HsapDv:0000055",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
      {
        label: "trisomy 18",
        ontology_term_id: "MONDO:0018071",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "9834221e-4c9f-41f0-bd88-f903b3cc489b",
    is_primary_data: "PRIMARY",
    name: "Survey of human embryonic development",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
    ],
    tissue: [
      {
        label: "adrenal gland",
        ontology_term_id: "UBERON:0002369",
      },
      {
        label: "cerebellum",
        ontology_term_id: "UBERON:0002037",
      },
      {
        label: "eye",
        ontology_term_id: "UBERON:0000970",
      },
      {
        label: "heart",
        ontology_term_id: "UBERON:0000948",
      },
      {
        label: "intestine",
        ontology_term_id: "UBERON:0000160",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
      {
        label: "liver",
        ontology_term_id: "UBERON:0002107",
      },
      {
        label: "lung",
        ontology_term_id: "UBERON:0002048",
      },
      {
        label: "muscle organ",
        ontology_term_id: "UBERON:0001630",
      },
      {
        label: "pancreas",
        ontology_term_id: "UBERON:0001264",
      },
      {
        label: "placenta",
        ontology_term_id: "UBERON:0001987",
      },
      {
        label: "spleen",
        ontology_term_id: "UBERON:0002106",
      },
      {
        label: "stomach",
        ontology_term_id: "UBERON:0000945",
      },
      {
        label: "telencephalon",
        ontology_term_id: "UBERON:0001893",
      },
      {
        label: "thymus",
        ontology_term_id: "UBERON:0002370",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "40bcddaa-801b-42ce-ae2d-53676ac02fed",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "6d679b7f-9c1c-4512-bcc1-85b1bbdf0531",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "b7db3b7f-8df4-4544-95c1-bdabfc50c85d",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "6e1135dd-aa1a-4fbb-8537-fd057ff0eadb",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "a82cd03d-65c4-47bd-9ff4-a1bd8fb5e24f",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "292fc973-4347-4bca-b6f7-fd3e264a4d7f",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "8af11b02-744e-4300-bdbd-9c8b05713749",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "609b672e-3af3-407f-ba7e-ed692ec8a618",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "71ccf9ad-07c7-4849-b212-753c84609b64",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "873cc5c7-8bac-486b-902a-cdc5c7a5c1ce",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "aae666d6-8b52-4b30-b076-942f36645da4",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "6d3c5c9d-ccf9-4a9c-b3ae-efc029f5335d",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
  {
    assay: [
      {
        label: "10x 3' v2",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3",
        ontology_term_id: "EFO:0009922",
      },
      {
        label: "10x 3' v3 (sci-plex)",
        ontology_term_id: "EFO:0009922 (sci-plex)",
      },
    ],
    cell_count: 320,
    cell_type: [
      {
        label: "endothelial cell",
        ontology_term_id: "CL:0000115",
      },
      {
        label: "fibroblast",
        ontology_term_id: "CL:0000057",
      },
    ],
    collection_id: "abf2ea2a-7bb8-44e5-9f0a-9ad65eb9eac1",
    development_stage: [
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "Theiler stage 01",
        ontology_term_id: "MmusDv:0000003",
      },
      {
        label: "prime adult stage",
        ontology_term_id: "UBERON:0018241",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    disease: [
      {
        label: "COVID-19",
        ontology_term_id: "MONDO:0100096",
      },
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "Yoruban",
        ontology_term_id: "HANCESTRO:0575",
      },
      {
        label: "na",
        ontology_term_id: "na",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    id: "e8bfe846-187e-460f-aca9-a9a76819a0f7",
    is_primary_data: "BOTH",
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Callithrix jacchus",
        ontology_term_id: "NCBITaxon:9483",
      },
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      {
        label: "Mus musculus",
        ontology_term_id: "NCBITaxon:10090",
      },
    ],
    schema_version: "2.0.0",
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
      {
        label: "unknown",
        ontology_term_id: "unknown",
      },
    ],
    tissue: [
      {
        label: "brain (organoid)",
        ontology_term_id: "UBERON:0000955 (organoid)",
      },
      {
        label: "epithelial cell",
        ontology_term_id: "CL:0000066",
      },
      {
        label: "fibroblast (cell culture)",
        ontology_term_id: "CL:0000057 (cell culture)",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "kidney",
        ontology_term_id: "UBERON:0002113",
      },
    ],
  },
];

// Adding unknown to handle schema and future filterable fields (development_stage, ethnicity etc)
export default datasetsIndex as unknown as DatasetResponse[];
