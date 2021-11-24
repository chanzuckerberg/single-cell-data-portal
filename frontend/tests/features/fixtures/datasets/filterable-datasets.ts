/* eslint-disable sonarjs/no-duplicate-string */

/**
 * Model of response from /datasets API endpoint.
 */

// App dependencies
import { FilterableDataset, IS_PRIMARY_DATA } from "src/common/entities";

// Data init
const filterableDatasets = [
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
    ],
    cell_count: 245389,
    cell_type: [
      {
        label: "amacrine cell",
        ontology_term_id: "CL:0000561",
      },
    ],
    disease: [
      {
        label: "type 2 diabetes mellitus",
        ontology_term_id: "MONDO:0005148",
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
    id: "003ab05f-a7a4-4dc5-8a5f-084bee4d9da7",
    is_primary_data: IS_PRIMARY_DATA.PRIMARY,
    name: "A single-cell transcriptomic atlas characterizes ageing tissues in the mouse",
    organism: [{ label: "Mus musculus", ontology_term_id: "NCBITaxon:10090" }],
    sex: [
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
    ],
    tissue: [
      {
        label: "adipose tissue",
        ontology_term_id: "UBERON:0001013",
      },
      {
        label: "mesenteric artery",
        ontology_term_id: "UBERON:0005616",
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
    cell_count: 11243,
    cell_type: [
      {
        label: "retinal ganglion cell",
        ontology_term_id: "CL:0000740",
      },
    ],
    disease: [
      {
        label: "type 2 diabetes mellitus",
        ontology_term_id: "MONDO:0005148",
      },
    ],
    ethnicity: [
      {
        label: "African American or Afro-Caribbean",
        ontology_term_id: "HANCESTRO:0016",
      },
      {
        label: "European",
        ontology_term_id: "HANCESTRO:0005",
      },
    ],
    id: "ddded24c-f3fb-4742-9811-f0a23bfa2a10",
    is_primary_data: IS_PRIMARY_DATA.PRIMARY,
    name: "scRNA-seq data analysis of endothelium-enriched mesenteric arterial tissues from human donors",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    sex: [
      {
        label: "female",
        ontology_term_id: "PATO:0000383",
      },
    ],
    tissue: [
      {
        label: "mesenteric artery",
        ontology_term_id: "UBERON:0005616",
      },
    ],
  },
  {
    assay: [
      {
        label: "Smart-seq2",
        ontology_term_id: "EFO:0008931",
      },
    ],
    cell_count: 3589,
    cell_type: [
      {
        label: "amacrine cell",
        ontology_term_id: "CL:0000561",
      },
    ],
    disease: [
      {
        label: "glioblastoma (disease)",
        ontology_term_id: "MONDO:0018177",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "8026b1f2-a29f-420a-88a0-3700a92236ef",
    is_primary_data: IS_PRIMARY_DATA.PRIMARY,
    name: "Infiltrating Neoplastic Cells Human Glioblastoma",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
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
        label: "cerebral cortex",
        ontology_term_id: "UBERON:0000956",
      },
    ],
  },
  {
    assay: [
      {
        label: "scATAC-seq",
        ontology_term_id: "EFO:0010891",
      },
    ],
    cell_count: 2000,
    cell_type: [
      {
        label: "respiratory goblet cell",
        ontology_term_id: "CL:0002370",
      },
      {
        label: "retina horizontal cell",
        ontology_term_id: "CL:0000745",
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
    id: "2250182e-f496-4fba-9461-dc6f282ae5a7",
    is_primary_data: IS_PRIMARY_DATA.PRIMARY,
    name: "Glutamatergic neurons \u2014 An Atlas of Gene Regulatory Elements in Adult Mouse Cerebrum.",
    organism: [{ label: "Mus musculus", ontology_term_id: "NCBITaxon:10090" }],
    sex: [
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
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
        label: "scATAC-seq",
        ontology_term_id: "EFO:0010891",
      },
    ],
    cell_count: 2000,
    cell_type: [
      {
        label: "respiratory goblet cell",
        ontology_term_id: "CL:0002370",
      },
      {
        label: "retina horizontal cell",
        ontology_term_id: "CL:0000745",
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
    id: "c2500381-f55a-4793-a2a2-2d72c9fecb1e",
    is_primary_data: IS_PRIMARY_DATA.SECONDARY,
    name: "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types: SMARTer_cells_MOp",
    organism: [{ label: "Mus musculus", ontology_term_id: "NCBITaxon:10090" }],
    sex: [
      {
        label: "male",
        ontology_term_id: "PATO:0000384",
      },
    ],
    tissue: [
      {
        label: "primary motor cortex",
        ontology_term_id: "UBERON:0001384",
      },
    ],
  },
] as FilterableDataset[];

export default filterableDatasets as unknown as FilterableDataset[];
