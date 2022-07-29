import { ChangeEvent, Dispatch, SetStateAction } from "react";
import { CellValue, Row } from "react-table";
import { EVENTS } from "src/common/analytics/events";
import { Collection, Ontology, PublisherMetadata } from "src/common/entities";
import { FilterCategoryKey } from "src/common/hooks/useCategoryFilter";

/**
 * Configuration model of category.
 */
export interface CategoryConfig {
  analyticsEvent?: EVENTS;
  categoryKey: CATEGORY_KEY; // TODO(cc) remove
  categoryType: CATEGORY_FILTER_TYPE; // TODO(cc) rename to filter type?
  filterCategoryKey: FILTER_CATEGORY_KEY;
  filterKey: FilterKey; // Key in result set row values to filter on.
  label: string;
  multiselect: boolean; // True if category can have multiple values selected.
  pinnedCategoryValues?: CATEGORY_VALUE_KEY[];
  tooltip?: string;
}

/**
 * Possible set of keys to filter over.
 */
export type FilterKey = keyof CollectionRow | keyof DatasetRow;

/**
 * Configuration model of ontology category.
 */
export interface OntologyCategoryConfig extends CategoryConfig {
  isLabelVisible: boolean; // True if ontology label is to be displayed (e.g. "Other Organisms")
  isSearchable: boolean; // True if ontology category search box is displayed
  isZerosVisible: boolean; // True if zero values are to be displayed
  ontologyTermSet: OntologyTermSet;
}

/**
 * Possible types of category filters, matches React Table's filter types.
 */
export enum CATEGORY_FILTER_TYPE {
  "BETWEEN" = "between",
  "INCLUDES_SOME" = "includesSome",
}

/**
 * Filterable metadata keys.
 * @deprecated - Ideally we want to remove this and use CATEGORY_CONFIGS (and type generated from it) instead.
 * TODO(cc) remove completely (check analytics)
 */
export enum CATEGORY_KEY {
  "ASSAY" = "assay",
  "CELL_COUNT" = "cell_count",
  "CELL_TYPE" = "cell_type",
  "CELL_TYPE_ANCESTORS" = "cell_type_ancestors",
  "DEVELOPMENT_STAGE_ANCESTORS" = "development_stage_ancestors",
  "DISEASE" = "disease",
  "ETHNICITY" = "ethnicity",
  "MEAN_GENES_PER_CELL" = "mean_genes_per_cell",
  "ORGANISM" = "organism",
  "PUBLICATION_AUTHORS" = "publicationAuthors",
  "PUBLICATION_DATE_VALUES" = "publicationDateValues",
  "SEX" = "sex",
  "TISSUE" = "tissue",
  "TISSUE_ANCESTORS" = "tissue_ancestors",
}

/**
 * Set of all categories to be included in filter.
 */
export enum FILTER_CATEGORY_KEY {
  "ASSAY" = "ASSAY",
  "CELL_COUNT" = "CELL_COUNT",
  "CELL_TYPE_DEPRECATED" = "CELL_TYPE_DEPRECATED",
  "CELL_TYPE" = "CELL_TYPE",
  "DEVELOPMENT_STAGE" = "DEVELOPMENT_STAGE",
  "DISEASE" = "DISEASE",
  "ETHNICITY" = "ETHNICITY",
  "GENE_COUNT" = "GENE_COUNT",
  "ORGANISM" = "ORGANISM",
  "PUBLICATION_AUTHORS" = "PUBLICATION_AUTHORS",
  "PUBLICATION_DATE_VALUES" = "PUBLICATION_DATE_VALUES",
  "SEX" = "SEX",
  "TISSUE_DEPRECATED" = "TISSUE_DEPRECATED", // TODO(cc) remove with #2569.
  "TISSUE" = "TISSUE",
}

/**
 * Category value keys.
 */
export enum CATEGORY_VALUE_KEY {
  NORMAL = "normal",
}

/**
 * Filterable values of datasets and collections.
 */
export interface Categories {
  assay: Ontology[];
  cell_type: Ontology[];
  cell_type_ancestors: string[];
  disease: Ontology[];
  development_stage_ancestors: string[];
  ethnicity: Ontology[];
  organism: Ontology[];
  sex: Ontology[];
  tissue: Ontology[]; // TODO(cc) remove with #2569.
  tissue_ancestors: string[];
}

/**
 * Keys of Categories that are ontology arrays.
 */
export type CategoriesKeyOfTypeOntologyArray = {
  [K in keyof Categories]: Categories[K] extends Ontology[] ? K : never;
}[keyof Categories];

/**
 * Join of collection, dataset and aggregated dataset information, optimized for filtering collections (that is,
 * datasets grouped by collection.
 */
export interface CollectionRow extends Categories, PublisherMetadataCategories {
  id: string;
  name: string;
  published_at: number;
  publisher_metadata?: PublisherMetadata; // Collection publication metadata
  recency: number; // Used by sort
  revised_at?: number;
  summaryCitation: string;
}

/*
 * Ontology view keys used to lookup ontology nodes for display . For example, for the ontology ID "HsapDv:0000045",
 * the key is "HsapDv".
 */
export enum ONTOLOGY_VIEW_KEY {
  "CL" = "CL",
  "HsapDv" = "HsapDv",
  "MmusDv" = "MmusDv",
  "UBERON" = "UBERON",
}

/**
 * Labels for displaying ontology views. Currently only used by development stage filter.
 */
export enum ONTOLOGY_VIEW_LABEL {
  "HsapDv" = "Homo Sapiens",
  "MmusDv" = "Mus Musculus",
  "UBERON" = "Other Organisms",
}

/**
 * Cell types to be included for display in cell type ontology tree.
 */
/* eslint-disable sort-keys -- disabling key order for readability. */
export const CELL_TYPE_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.CL]: [
    { label: "bladder cell", ontology_term_id: "CL:1001319" },
    {
      label: "brush cell",
      ontology_term_id: "CL:0002204",
      children: [
        {
          label: "brush cell of trachebronchial tree",
          ontology_term_id: "CL:0002075",
          children: [
            {
              label: "brush cell of bronchus",
              ontology_term_id: "CL:0002208",
            },
            {
              label: "brush cell of trachea",
              ontology_term_id: "CL:0002207",
            },
          ],
        },
      ],
    },
    {
      label: "cardiac muscle cell",
      ontology_term_id: "CL:0000746",
      children: [
        {
          label: "regular cardiac myocyte",
          ontology_term_id: "CL:0002098",
          children: [
            {
              label: "regular atrial cardiac myocyte",
              ontology_term_id: "CL:0002129",
            },
            {
              label: "regular ventricular cardiac myocyte",
              ontology_term_id: "CL:0002131",
            },
          ],
        },
        {
          label: "ventricular cardiac muscle cell",
          ontology_term_id: "CL:2000046",
        },
      ],
    },
    {
      label: "cell in vitro",
      ontology_term_id: "CL:0001034",
      children: [
        {
          label: "experimentally modified cell in vitro",
          ontology_term_id: "CL:0000578",
          children: [
            {
              label: "cultured cell",
              ontology_term_id: "CL:0000010",
            },
          ],
        },
      ],
    },
    {
      label: "cell of skeletal muscle",
      ontology_term_id: "CL:0000188",
      children: [
        {
          label: "skeletal muscle fiber",
          ontology_term_id: "CL:0008002",
          children: [
            {
              label: "fast muscle cell",
              ontology_term_id: "CL:0000190",
            },
            {
              label: "slow muscle cell",
              ontology_term_id: "CL:0000189",
            },
            {
              label: "tongue muscle cell",
              ontology_term_id: "CL:0002673",
            },
          ],
        },
      ],
    },
    {
      label: "chorionic trophoblast cell",
      ontology_term_id: "CL:0011101",
    },
    {
      label: "ciliated cell",
      ontology_term_id: "CL:0000064",
      children: [
        {
          label: "ciliated epithelial cell",
          ontology_term_id: "CL:0000067",
          children: [
            {
              label: "multi-ciliated epithelial cell",
              ontology_term_id: "CL:0005012",
            },
          ],
        },
        {
          label: "lung ciliated cell",
          ontology_term_id: "CL:1000271",
        },
      ],
    },
    {
      label: "conjunctival epithelial cell",
      ontology_term_id: "CL:1000432",
    },
    {
      label: "connective tissue cell",
      ontology_term_id: "CL:0002320",
      children: [
        {
          label: "adipose microvascular endothelial cell",
          ontology_term_id: "CL:2000072",
        },
        { label: "adventitial cell", ontology_term_id: "CL:0002503" },
        {
          label: "fat cell",
          ontology_term_id: "CL:0000136",
          children: [
            {
              label: "epicardial adipocyte",
              ontology_term_id: "CL:1000309",
            },
            {
              label: "subcutaneous fat cell",
              ontology_term_id: "CL:0002521",
            },
          ],
        },
        {
          label: "fibroblast",
          ontology_term_id: "CL:0000057",
          children: [
            {
              label: "bronchus fibroblast of lung",
              ontology_term_id: "CL:2000093",
            },
            {
              label: "fibroblast of breast",
              ontology_term_id: "CL:4006000",
            },
            {
              label: "fibroblast of cardiac tissue",
              ontology_term_id: "CL:0002548",
            },
            {
              label:
                "fibroblast of connective tissue of glandular part of prostate",
              ontology_term_id: "CL:1000305",
            },
            {
              label:
                "fibroblast of connective tissue of nonglandular part of prostate",
              ontology_term_id: "CL:1000304",
            },
            {
              label: "fibroblast of connective tissue of prostate",
              ontology_term_id: "CL:1000299",
            },
            {
              label: "fibroblast of lung",
              ontology_term_id: "CL:0002553",
            },
            { label: "keratocyte", ontology_term_id: "CL:0002363" },
            {
              label: "kidney interstitial fibroblast",
              ontology_term_id: "CL:1000692",
            },
            {
              label: "pancreatic stellate cell",
              ontology_term_id: "CL:0002410",
            },
            {
              label: "pulmonary interstitial fibroblast",
              ontology_term_id: "CL:0002241",
            },
            {
              label: "reticular cell",
              ontology_term_id: "CL:0000432",
            },
            {
              label: "skin fibroblast",
              ontology_term_id: "CL:0002620",
            },
            {
              label: "vascular leptomeningeal cell",
              ontology_term_id: "CL:4023051",
            },
          ],
        },
        {
          label: "pericyte",
          ontology_term_id: "CL:0000669",
          children: [
            {
              label: "brain pericyte",
              ontology_term_id: "CL:2000043",
            },
            {
              label: "mesangial cell",
              ontology_term_id: "CL:0000650",
            },
            {
              label: "renal interstitial pericyte",
              ontology_term_id: "CL:1001318",
            },
          ],
        },
        {
          label: "stromal cell",
          ontology_term_id: "CL:0000499",
          children: [
            { label: "fibrocyte", ontology_term_id: "CL:0000135" },
            {
              label: "hepatic stellate cell",
              ontology_term_id: "CL:0000632",
            },
            {
              label: "leptomeningeal cell",
              ontology_term_id: "CL:0000708",
            },
            {
              label: "prostate stromal cell",
              ontology_term_id: "CL:0002622",
            },
            {
              label: "stromal cell of ovary",
              ontology_term_id: "CL:0002132",
            },
            { label: "tendon cell", ontology_term_id: "CL:0000388" },
          ],
        },
      ],
    },
    { label: "contractile cell", ontology_term_id: "CL:0000183" },
    {
      label: "corneal epithelial cell",
      ontology_term_id: "CL:0000575",
    },
    {
      label: "duct epithelial cell",
      ontology_term_id: "CL:0000068",
      children: [
        {
          label: "basal epithelial cell of prostatic duct",
          ontology_term_id: "CL:0002236",
        },
        {
          label: "cholangiocyte",
          ontology_term_id: "CL:1000488",
          children: [
            {
              label: "intrahepatic cholangiocyte",
              ontology_term_id: "CL:0002538",
            },
          ],
        },
        {
          label: "epithelial cell of lacrimal sac",
          ontology_term_id: "CL:1000436",
        },
        {
          label: "kidney collecting duct epithelial cell",
          ontology_term_id: "CL:1000454",
          children: [
            {
              label: "kidney collecting duct intercalated cell",
              ontology_term_id: "CL:1001432",
            },
            {
              label: "kidney collecting duct principal cell",
              ontology_term_id: "CL:1001431",
            },
            {
              label: "renal alpha-intercalated cell",
              ontology_term_id: "CL:0005011",
            },
            {
              label: "renal beta-intercalated cell",
              ontology_term_id: "CL:0002201",
            },
          ],
        },
        {
          label: "luminal epithelial cell of mammary gland",
          ontology_term_id: "CL:0002326",
        },
      ],
    },
    {
      label: "endothelial cell",
      ontology_term_id: "CL:0000115",
      children: [
        {
          label: "cardiac endothelial cell",
          ontology_term_id: "CL:0010008",
          children: [
            {
              label: "endocardial cell",
              ontology_term_id: "CL:0002350",
            },
          ],
        },
        {
          label: "cerebral cortex endothelial cell",
          ontology_term_id: "CL:1001602",
        },
        {
          label: "endothelial cell of hepatic sinusoid",
          ontology_term_id: "CL:1000398",
          children: [
            {
              label: "endothelial cell of pericentral hepatic sinusoid",
              ontology_term_id: "CL:0019022",
            },
            {
              label: "endothelial cell of periportal hepatic sinusoid",
              ontology_term_id: "CL:0019021",
            },
          ],
        },
        {
          label: "endothelial cell of vascular tree",
          ontology_term_id: "CL:0002139",
          children: [
            {
              label: "aortic endothelial cell",
              ontology_term_id: "CL:0002544",
            },
            {
              label: "blood vessel endothelial cell",
              ontology_term_id: "CL:0000071",
            },
            {
              label: "capillary endothelial cell",
              ontology_term_id: "CL:0002144",
            },
            {
              label: "dermis microvascular lymphatic vessel endothelial cell",
              ontology_term_id: "CL:2000041",
            },
            {
              label: "endothelial cell of artery",
              ontology_term_id: "CL:1000413",
            },
            {
              label: "endothelial cell of coronary artery",
              ontology_term_id: "CL:2000018",
            },
            {
              label: "endothelial cell of lymphatic vessel",
              ontology_term_id: "CL:0002138",
            },
            {
              label: "endothelial cell of umbilical vein",
              ontology_term_id: "CL:0002618",
            },
            {
              label: "endothelial stalk cell",
              ontology_term_id: "CL:0002671",
            },
            {
              label: "kidney capillary endothelial cell",
              ontology_term_id: "CL:1000892",
            },
            {
              label: "lung endothelial cell",
              ontology_term_id: "CL:1001567",
            },
            {
              label: "lung microvascular endothelial cell",
              ontology_term_id: "CL:2000016",
            },
            {
              label: "peritubular capillary endothelial cell",
              ontology_term_id: "CL:1001033",
            },
            {
              label: "prostate gland microvascular endothelial cell",
              ontology_term_id: "CL:2000059",
            },
            {
              label: "pulmonary artery endothelial cell",
              ontology_term_id: "CL:1001568",
            },
            {
              label: "retinal blood vessel endothelial cell",
              ontology_term_id: "CL:0002585",
            },
            {
              label: "vasa recta ascending limb cell",
              ontology_term_id: "CL:1001131",
            },
            {
              label: "vasa recta descending limb cell",
              ontology_term_id: "CL:1001285",
            },
            {
              label: "vein endothelial cell",
              ontology_term_id: "CL:0002543",
            },
          ],
        },
        {
          label: "glomerular endothelial cell",
          ontology_term_id: "CL:0002188",
        },
        {
          label: "gut endothelial cell",
          ontology_term_id: "CL:0000131",
        },
      ],
    },
    { label: "ependymal cell", ontology_term_id: "CL:0000065" },
    {
      label: "epidermal cell",
      ontology_term_id: "CL:0000362",
      children: [
        {
          label: "taste receptor cell",
          ontology_term_id: "CL:0000209",
        },
      ],
    },
    {
      label: "epithelial cell",
      ontology_term_id: "CL:0000066",
      children: [
        {
          label: "epithelial cell of lung",
          ontology_term_id: "CL:0000082",
          children: [
            {
              label: "epithelial cell of alveolus of lung",
              ontology_term_id: "CL:0010003",
            },
            { label: "pneumocyte", ontology_term_id: "CL:0000322" },
            {
              label: "type I pneumocyte",
              ontology_term_id: "CL:0002062",
            },
          ],
        },
        {
          label: "squamous epithelial cell",
          ontology_term_id: "CL:0000076",
        },
      ],
    },
    {
      label: "epithelial cell of exocrine pancreas",
      ontology_term_id: "CL:1001433",
      children: [
        {
          label: "pancreatic ductal cell",
          ontology_term_id: "CL:0002079",
        },
      ],
    },
    {
      label: "epithelial cell of lower respiratory tract",
      ontology_term_id: "CL:0002632",
      children: [
        {
          label: "basal epithelial cell of tracheobronchial tree",
          ontology_term_id: "CL:0002329",
        },
        {
          label: "ciliated cell of the bronchus",
          ontology_term_id: "CL:0002332",
        },
        {
          label: "ciliated columnar cell of tracheobronchial tree",
          ontology_term_id: "CL:0002145",
        },
      ],
    },
    {
      label: "epithelial cell of prostate",
      ontology_term_id: "CL:0002231",
      children: [
        {
          label: "basal cell of prostate epithelium",
          ontology_term_id: "CL:0002341",
        },
        {
          label: "luminal cell of prostate epithelium",
          ontology_term_id: "CL:0002340",
        },
      ],
    },
    {
      label: "epithelial cell of proximal tubule",
      ontology_term_id: "CL:0002306",
      children: [
        {
          label: "kidney proximal convoluted tubule epithelial cell",
          ontology_term_id: "CL:1000838",
        },
        {
          label: "kidney proximal straight tubule epithelial cell",
          ontology_term_id: "CL:1000839",
        },
      ],
    },
    {
      label: "epithelial cell of sweat gland",
      ontology_term_id: "CL:1000448",
    },
    {
      label: "epithelial cell of thymus",
      ontology_term_id: "CL:0002293",
      children: [
        {
          label: "medullary thymic epithelial cell",
          ontology_term_id: "CL:0002365",
        },
      ],
    },
    {
      label: "epithelial cell of urethra",
      ontology_term_id: "CL:1000296",
      children: [
        {
          label: "urethra urothelial cell",
          ontology_term_id: "CL:1001430",
        },
      ],
    },
    {
      label: "epithelial cell of uterus",
      ontology_term_id: "CL:0002149",
    },
    {
      label: "extravillous trophoblast",
      ontology_term_id: "CL:0008036",
    },
    { label: "fenestrated cell", ontology_term_id: "CL:0000666" },
    {
      label: "follicular dendritic cell",
      ontology_term_id: "CL:0000442",
    },
    { label: "gut absorptive cell", ontology_term_id: "CL:0000677" },
    {
      label: "hematopoietic cell",
      ontology_term_id: "CL:0000988",
      children: [
        { label: "blood cell", ontology_term_id: "CL:0000081" },
        {
          label: "hematopoietic precursor cell",
          ontology_term_id: "CL:0008001",
          children: [
            {
              label: "common myeloid progenitor",
              ontology_term_id: "CL:0000049",
            },
          ],
        },
        {
          label: "leukocyte",
          ontology_term_id: "CL:0000738",
          children: [
            { label: "B cell", ontology_term_id: "CL:0000236" },
            {
              label: "CD16-negative, CD56-bright natural killer cell, human",
              ontology_term_id: "CL:0000938",
            },
            {
              label: "CD16-positive, CD56-dim natural killer cell, human",
              ontology_term_id: "CL:0000939",
            },
            {
              label: "CD4-positive helper T cell",
              ontology_term_id: "CL:0000492",
            },
            {
              label:
                "CD4-positive, CD25-positive, alpha-beta regulatory T cell",
              ontology_term_id: "CL:0000792",
            },
            {
              label: "CD4-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000624",
            },
            {
              label: "CD4-positive, alpha-beta cytotoxic T cell",
              ontology_term_id: "CL:0000934",
            },
            {
              label: "CD4-positive, alpha-beta memory T cell",
              ontology_term_id: "CL:0000897",
            },
            {
              label: "CD4-positive, alpha-beta memory T cell, CD45RO-positive",
              ontology_term_id: "CL:0001204",
            },
            {
              label: "CD4-positive, alpha-beta thymocyte",
              ontology_term_id: "CL:0000810",
            },
            {
              label:
                "CD8-alpha alpha positive, gamma-delta intraepithelial T cell",
              ontology_term_id: "CL:0000802",
            },
            {
              label: "CD8-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000625",
            },
            {
              label:
                "CD8-positive, alpha-beta cytokine secreting effector T cell",
              ontology_term_id: "CL:0000908",
            },
            {
              label: "CD8-positive, alpha-beta cytotoxic T cell",
              ontology_term_id: "CL:0000794",
            },
            {
              label: "CD8-positive, alpha-beta memory T cell",
              ontology_term_id: "CL:0000909",
            },
            {
              label: "CD8-positive, alpha-beta memory T cell, CD45RO-positive",
              ontology_term_id: "CL:0001203",
            },
            {
              label: "CD8-positive, alpha-beta thymocyte",
              ontology_term_id: "CL:0000811",
            },
            {
              label: "DN3 thymocyte",
              ontology_term_id: "CL:0000807",
            },
            {
              label: "DN4 thymocyte",
              ontology_term_id: "CL:0000808",
            },
            { label: "ILC1, human", ontology_term_id: "CL:0001077" },
            {
              label: "IgA plasma cell",
              ontology_term_id: "CL:0000987",
            },
            {
              label: "IgA plasmablast",
              ontology_term_id: "CL:0000984",
            },
            {
              label: "IgG memory B cell",
              ontology_term_id: "CL:0000979",
            },
            {
              label: "IgG plasma cell",
              ontology_term_id: "CL:0000985",
            },
            {
              label: "IgG plasmablast",
              ontology_term_id: "CL:0000982",
            },
            {
              label: "IgG-negative class switched memory B cell",
              ontology_term_id: "CL:0002117",
            },
            {
              label: "IgM plasma cell",
              ontology_term_id: "CL:0000986",
            },
            { label: "T cell", ontology_term_id: "CL:0000084" },
            {
              label: "T follicular helper cell",
              ontology_term_id: "CL:0002038",
            },
            {
              label: "T-helper 1 cell",
              ontology_term_id: "CL:0000545",
            },
            {
              label: "T-helper 17 cell",
              ontology_term_id: "CL:0000899",
            },
            {
              label: "T-helper 2 cell",
              ontology_term_id: "CL:0000546",
            },
            {
              label: "T-helper 22 cell",
              ontology_term_id: "CL:0001042",
            },
            {
              label: "activated CD4-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000896",
            },
            {
              label: "activated CD4-positive, alpha-beta T cell, human",
              ontology_term_id: "CL:0001043",
            },
            {
              label: "activated CD8-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000906",
            },
            {
              label: "activated CD8-positive, alpha-beta T cell, human",
              ontology_term_id: "CL:0001049",
            },
            {
              label: "alpha-beta T cell",
              ontology_term_id: "CL:0000789",
            },
            {
              label: "central memory CD4-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000904",
            },
            {
              label: "central memory CD8-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000907",
            },
            {
              label: "class switched memory B cell",
              ontology_term_id: "CL:0000972",
            },
            {
              label: "conventional dendritic cell",
              ontology_term_id: "CL:0000990",
            },
            {
              label: "dendritic cell",
              ontology_term_id: "CL:0000451",
            },
            {
              label: "dendritic cell, human",
              ontology_term_id: "CL:0001056",
            },
            {
              label: "double negative T regulatory cell",
              ontology_term_id: "CL:0011024",
            },
            {
              label: "double negative thymocyte",
              ontology_term_id: "CL:0002489",
            },
            {
              label: "double-positive, alpha-beta thymocyte",
              ontology_term_id: "CL:0000809",
            },
            {
              label: "effector CD4-positive, alpha-beta T cell",
              ontology_term_id: "CL:0001044",
            },
            {
              label: "effector CD8-positive, alpha-beta T cell",
              ontology_term_id: "CL:0001050",
            },
            {
              label: "effector memory CD4-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000905",
            },
            {
              label: "effector memory CD8-positive, alpha-beta T cell",
              ontology_term_id: "CL:0000913",
            },
            {
              label:
                "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
              ontology_term_id: "CL:0001062",
            },
            {
              label: "exhausted T cell",
              ontology_term_id: "CL:0011025",
            },
            {
              label: "follicular B cell",
              ontology_term_id: "CL:0000843",
            },
            {
              label: "gamma-delta T cell",
              ontology_term_id: "CL:0000798",
            },
            {
              label: "germinal center B cell",
              ontology_term_id: "CL:0000844",
            },
            {
              label: "group 2 innate lymphoid cell, human",
              ontology_term_id: "CL:0001081",
            },
            {
              label: "group 3 innate lymphoid cell",
              ontology_term_id: "CL:0001071",
            },
            {
              label: "immature B cell",
              ontology_term_id: "CL:0000816",
            },
            {
              label: "immature NK T cell",
              ontology_term_id: "CL:0000914",
            },
            {
              label: "immature T cell",
              ontology_term_id: "CL:0002420",
            },
            {
              label: "immature alpha-beta T cell",
              ontology_term_id: "CL:0000790",
            },
            {
              label: "immature natural killer cell",
              ontology_term_id: "CL:0000823",
            },
            {
              label: "innate lymphoid cell",
              ontology_term_id: "CL:0001065",
            },
            {
              label: "late pro-B cell",
              ontology_term_id: "CL:0002048",
            },
            {
              label: "liver dendritic cell",
              ontology_term_id: "CL:2000055",
            },
            { label: "lymphocyte", ontology_term_id: "CL:0000542" },
            {
              label: "mature B cell",
              ontology_term_id: "CL:0000785",
            },
            {
              label: "mature NK T cell",
              ontology_term_id: "CL:0000814",
            },
            {
              label: "mature alpha-beta T cell",
              ontology_term_id: "CL:0000791",
            },
            {
              label: "mature gamma-delta T cell",
              ontology_term_id: "CL:0000800",
            },
            {
              label: "memory B cell",
              ontology_term_id: "CL:0000787",
            },
            {
              label: "memory T cell",
              ontology_term_id: "CL:0000813",
            },
            {
              label: "mononuclear cell",
              ontology_term_id: "CL:0000842",
            },
            {
              label: "mucosal invariant T cell",
              ontology_term_id: "CL:0000940",
            },
            { label: "naive B cell", ontology_term_id: "CL:0000788" },
            { label: "naive T cell", ontology_term_id: "CL:0000898" },
            {
              label: "naive regulatory T cell",
              ontology_term_id: "CL:0002677",
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
              label: "peripheral blood mononuclear cell",
              ontology_term_id: "CL:2000001",
            },
            { label: "plasma cell", ontology_term_id: "CL:0000786" },
            { label: "plasmablast", ontology_term_id: "CL:0000980" },
            {
              label: "plasmacytoid dendritic cell",
              ontology_term_id: "CL:0000784",
            },
            {
              label: "plasmacytoid dendritic cell, human",
              ontology_term_id: "CL:0001058",
            },
            {
              label: "precursor B cell",
              ontology_term_id: "CL:0000817",
            },
            {
              label: "professional antigen presenting cell",
              ontology_term_id: "CL:0000145",
            },
            {
              label: "regulatory T cell",
              ontology_term_id: "CL:0000815",
            },
            { label: "thymocyte", ontology_term_id: "CL:0000893" },
            {
              label: "transitional stage B cell",
              ontology_term_id: "CL:0000818",
            },
            {
              label: "type I NK T cell",
              ontology_term_id: "CL:0000921",
            },
            {
              label: "unswitched memory B cell",
              ontology_term_id: "CL:0000970",
            },
          ],
        },
        {
          label: "myeloid cell",
          ontology_term_id: "CL:0000763",
          children: [
            {
              label: "CD141-positive myeloid dendritic cell",
              ontology_term_id: "CL:0002394",
            },
            {
              label: "CD1c-positive myeloid dendritic cell",
              ontology_term_id: "CL:0002399",
            },
            { label: "Kupffer cell", ontology_term_id: "CL:0000091" },
            {
              label: "Langerhans cell",
              ontology_term_id: "CL:0000453",
            },
            {
              label: "alternatively activated macrophage",
              ontology_term_id: "CL:0000890",
            },
            {
              label: "alveolar macrophage",
              ontology_term_id: "CL:0000583",
            },
            { label: "basophil", ontology_term_id: "CL:0000767" },
            {
              label: "colon macrophage",
              ontology_term_id: "CL:0009038",
            },
            {
              label: "elicited macrophage",
              ontology_term_id: "CL:0000861",
            },
            {
              label: "enucleate erythrocyte",
              ontology_term_id: "CL:0000595",
            },
            {
              label: "enucleated reticulocyte",
              ontology_term_id: "CL:0002422",
            },
            {
              label: "epidermal Langerhans cell",
              ontology_term_id: "CL:0002457",
            },
            { label: "erythrocyte", ontology_term_id: "CL:0000232" },
            {
              label: "erythroid lineage cell",
              ontology_term_id: "CL:0000764",
            },
            { label: "granulocyte", ontology_term_id: "CL:0000094" },
            {
              label: "immature neutrophil",
              ontology_term_id: "CL:0000776",
            },
            {
              label: "inflammatory macrophage",
              ontology_term_id: "CL:0000863",
            },
            {
              label: "kidney resident macrophage",
              ontology_term_id: "CL:1000698",
            },
            {
              label: "lung macrophage",
              ontology_term_id: "CL:1001603",
            },
            { label: "macrophage", ontology_term_id: "CL:0000235" },
            {
              label: "mature conventional dendritic cell",
              ontology_term_id: "CL:0000841",
            },
            {
              label: "mature microglial cell",
              ontology_term_id: "CL:0002629",
            },
            {
              label: "megakaryocyte",
              ontology_term_id: "CL:0000556",
            },
            {
              label: "microglial cell",
              ontology_term_id: "CL:0000129",
            },
            {
              label: "multinucleated giant cell",
              ontology_term_id: "CL:0000647",
            },
            {
              label: "myeloid dendritic cell",
              ontology_term_id: "CL:0000782",
            },
            {
              label: "myeloid dendritic cell, human",
              ontology_term_id: "CL:0001057",
            },
            {
              label: "myeloid leukocyte",
              ontology_term_id: "CL:0000766",
            },
            { label: "neutrophil", ontology_term_id: "CL:0000775" },
            {
              label: "proerythroblast",
              ontology_term_id: "CL:0000547",
            },
            {
              label: "tissue-resident macrophage",
              ontology_term_id: "CL:0000864",
            },
          ],
        },
      ],
    },
    {
      label: "hepatocyte",
      ontology_term_id: "CL:0000182",
      children: [
        {
          label: "centrilobular region hepatocyte",
          ontology_term_id: "CL:0019029",
        },
        {
          label: "midzonal region hepatocyte",
          ontology_term_id: "CL:0019028",
        },
        {
          label: "periportal region hepatocyte",
          ontology_term_id: "CL:0019026",
        },
      ],
    },
    { label: "inflammatory cell", ontology_term_id: "CL:0009002" },
    {
      label: "interstitial cell of Cajal",
      ontology_term_id: "CL:0002088",
    },
    {
      label: "intestinal epithelial cell",
      ontology_term_id: "CL:0002563",
      children: [
        { label: "enterocyte", ontology_term_id: "CL:0000584" },
        {
          label: "enterocyte of epithelium of small intestine",
          ontology_term_id: "CL:1000334",
          children: [
            {
              label: "enterocyte of epithelium proper of ileum",
              ontology_term_id: "CL:1000342",
            },
          ],
        },
        {
          label: "epithelial cell of large intestine",
          ontology_term_id: "CL:0002253",
          children: [
            {
              label: "brush cell of epithelium proper of large intestine",
              ontology_term_id: "CL:0002203",
            },
            {
              label: "enterocyte of epithelium of large intestine",
              ontology_term_id: "CL:0002071",
            },
          ],
        },
        {
          label: "intestinal tuft cell",
          ontology_term_id: "CL:0019032",
        },
      ],
    },
    {
      label: "ionocyte",
      ontology_term_id: "CL:0005006",
      children: [
        {
          label: "pulmonary ionocyte",
          ontology_term_id: "CL:0017000",
        },
      ],
    },
    {
      label: "kidney cell",
      ontology_term_id: "CL:1000497",
      children: [
        {
          label: "kidney collecting duct cell",
          ontology_term_id: "CL:1001225",
          children: [
            {
              label: "kidney inner medulla collecting duct epithelial cell",
              ontology_term_id: "CL:1000547",
            },
          ],
        },
        {
          label: "kidney epithelial cell",
          ontology_term_id: "CL:0002518",
          children: [
            {
              label: "epithelial cell of glomerular capsule",
              ontology_term_id: "CL:1000450",
            },
            {
              label: "epithelial cell of nephron",
              ontology_term_id: "CL:1000449",
            },
            {
              label: "kidney connecting tubule epithelial cell",
              ontology_term_id: "CL:1000768",
            },
            {
              label: "kidney cortex artery cell",
              ontology_term_id: "CL:1001045",
            },
            {
              label: "kidney loop of Henle epithelial cell",
              ontology_term_id: "CL:1000909",
            },
            {
              label:
                "kidney loop of Henle thin descending limb epithelial cell",
              ontology_term_id: "CL:1001111",
            },
            {
              label: "nephron tubule epithelial cell",
              ontology_term_id: "CL:1000494",
            },
            {
              label: "parietal epithelial cell",
              ontology_term_id: "CL:1000452",
            },
          ],
        },
        {
          label: "kidney interstitial cell",
          ontology_term_id: "CL:1000500",
        },
        {
          label: "papillary tips cell",
          ontology_term_id: "CL:1000597",
        },
      ],
    },
    {
      label: "kidney distal convoluted tubule epithelial cell",
      ontology_term_id: "CL:1000849",
    },
    {
      label: "kidney loop of Henle ascending limb epithelial cell",
      ontology_term_id: "CL:1001016",
      children: [
        {
          label: "kidney loop of Henle thick ascending limb epithelial cell",
          ontology_term_id: "CL:1001106",
        },
        {
          label: "kidney loop of Henle thin ascending limb epithelial cell",
          ontology_term_id: "CL:1001107",
        },
      ],
    },
    { label: "lens fiber cell", ontology_term_id: "CL:0011004" },
    {
      label: "macula densa epithelial cell",
      ontology_term_id: "CL:1000850",
    },
    {
      label: "melanocyte",
      ontology_term_id: "CL:0000148",
      children: [
        {
          label: "melanocyte of skin",
          ontology_term_id: "CL:1000458",
        },
      ],
    },
    { label: "mesenchymal cell", ontology_term_id: "CL:0008019" },
    {
      label: "mesothelial cell",
      ontology_term_id: "CL:0000077",
      children: [
        {
          label: "mesothelial cell of pleura",
          ontology_term_id: "CL:1000491",
          children: [
            {
              label: "mesothelial cell of visceral pleura",
              ontology_term_id: "CL:1000493",
            },
          ],
        },
      ],
    },
    {
      label: "mononuclear phagocyte",
      ontology_term_id: "CL:0000113",
    },
    {
      label: "muscle cell",
      ontology_term_id: "CL:0000187",
      children: [
        {
          label: "smooth muscle cell",
          ontology_term_id: "CL:0000192",
          children: [
            {
              label: "aortic smooth muscle cell",
              ontology_term_id: "CL:0002539",
            },
            {
              label: "blood vessel smooth muscle cell",
              ontology_term_id: "CL:0019018",
            },
            {
              label: "bronchial smooth muscle cell",
              ontology_term_id: "CL:0002598",
            },
            {
              label: "enteric smooth muscle cell",
              ontology_term_id: "CL:0002504",
            },
            {
              label: "myometrial cell",
              ontology_term_id: "CL:0002366",
            },
            {
              label: "smooth muscle cell of large intestine",
              ontology_term_id: "CL:1000279",
            },
            {
              label: "smooth muscle cell of prostate",
              ontology_term_id: "CL:1000487",
            },
            {
              label: "smooth muscle cell of small intestine",
              ontology_term_id: "CL:1000275",
            },
            {
              label: "smooth muscle cell of the pulmonary artery",
              ontology_term_id: "CL:0002591",
            },
            {
              label: "smooth muscle cell of trachea",
              ontology_term_id: "CL:0002600",
            },
            {
              label: "smooth muscle fiber of ileum",
              ontology_term_id: "CL:1000278",
            },
            {
              label: "tracheobronchial smooth muscle cell",
              ontology_term_id: "CL:0019019",
            },
            {
              label: "uterine smooth muscle cell",
              ontology_term_id: "CL:0002601",
            },
            {
              label: "vascular associated smooth muscle cell",
              ontology_term_id: "CL:0000359",
            },
          ],
        },
      ],
    },
    { label: "myoepithelial cell", ontology_term_id: "CL:0000185" },
    {
      label: "myoepithelial cell of mammary gland",
      ontology_term_id: "CL:0002324",
    },
    { label: "myofibroblast cell", ontology_term_id: "CL:0000186" },
    {
      label: "neoplastic cell",
      ontology_term_id: "CL:0001063",
      children: [{ label: "malignant cell", ontology_term_id: "CL:0001064" }],
    },
    {
      label: "neural cell",
      ontology_term_id: "CL:0002319",
      children: [
        { label: "amacrine cell", ontology_term_id: "CL:0000561" },
        {
          label: "glial cell",
          ontology_term_id: "CL:0000125",
          children: [
            {
              label: "Bergmann glial cell",
              ontology_term_id: "CL:0000644",
            },
            { label: "Mueller cell", ontology_term_id: "CL:0000636" },
            { label: "Schwann cell", ontology_term_id: "CL:0002573" },
            { label: "astrocyte", ontology_term_id: "CL:0000127" },
            {
              label: "astrocyte of the cerebellum",
              ontology_term_id: "CL:0002603",
            },
            {
              label: "astrocyte of the cerebral cortex",
              ontology_term_id: "CL:0002605",
            },
            {
              label: "forebrain radial glial cell",
              ontology_term_id: "CL:0013000",
            },
            {
              label: "macroglial cell",
              ontology_term_id: "CL:0000126",
            },
            {
              label: "mature astrocyte",
              ontology_term_id: "CL:0002627",
            },
            {
              label: "oligodendrocyte",
              ontology_term_id: "CL:0000128",
            },
            {
              label: "radial glial cell",
              ontology_term_id: "CL:0000681",
            },
          ],
        },
        {
          label: "neuron",
          ontology_term_id: "CL:0000540",
          children: [
            {
              label: "Cajal-Retzius cell",
              ontology_term_id: "CL:0000695",
            },
            {
              label: "Purkinje cell",
              ontology_term_id: "CL:0000121",
            },
            {
              label: "bipolar neuron",
              ontology_term_id: "CL:0000103",
            },
            {
              label: "cardiac neuron",
              ontology_term_id: "CL:0010022",
            },
            {
              label: "cerebellar granule cell",
              ontology_term_id: "CL:0001031",
            },
            {
              label: "enteric neuron",
              ontology_term_id: "CL:0007011",
            },
            {
              label: "eye photoreceptor cell",
              ontology_term_id: "CL:0000287",
            },
            {
              label: "ganglion interneuron",
              ontology_term_id: "CL:0000397",
            },
            { label: "granule cell", ontology_term_id: "CL:0000120" },
            {
              label: "inhibitory interneuron",
              ontology_term_id: "CL:0000498",
            },
            {
              label: "inhibitory motor neuron",
              ontology_term_id: "CL:0008015",
            },
            { label: "interneuron", ontology_term_id: "CL:0000099" },
            {
              label: "medium spiny neuron",
              ontology_term_id: "CL:1001474",
            },
            { label: "motor neuron", ontology_term_id: "CL:0000100" },
            {
              label: "neuronal brush cell",
              ontology_term_id: "CL:0000555",
            },
            {
              label: "photoreceptor cell",
              ontology_term_id: "CL:0000210",
            },
            {
              label: "retinal cone cell",
              ontology_term_id: "CL:0000573",
            },
            {
              label: "retinal ganglion cell",
              ontology_term_id: "CL:0000740",
            },
            {
              label: "retinal rod cell",
              ontology_term_id: "CL:0000604",
            },
            {
              label: "stellate neuron",
              ontology_term_id: "CL:0000122",
            },
            {
              label: "sympathetic neuron",
              ontology_term_id: "CL:0011103",
            },
            {
              label: "visceromotor neuron",
              ontology_term_id: "CL:0005025",
            },
          ],
        },
        {
          label: "neuron associated cell (sensu Vertebrata)",
          ontology_term_id: "CL:0000123",
        },
        {
          label: "retina horizontal cell",
          ontology_term_id: "CL:0000745",
        },
        {
          label: "retinal pigment epithelial cell",
          ontology_term_id: "CL:0002586",
        },
      ],
    },
    {
      label: "neural progenitor cell",
      ontology_term_id: "CL:0011020",
    },
    { label: "phagocyte", ontology_term_id: "CL:0000234" },
    {
      label: "pigmented ciliary epithelial cell",
      ontology_term_id: "CL:0002303",
    },
    {
      label: "podocyte (sensu Diptera)",
      ontology_term_id: "CL:0000391",
    },
    { label: "primordial germ cell", ontology_term_id: "CL:0000670" },
    {
      label: "progenitor cell",
      ontology_term_id: "CL:0011026",
      children: [
        { label: "erythroblast", ontology_term_id: "CL:0000765" },
        {
          label: "granulocyte monocyte progenitor cell",
          ontology_term_id: "CL:0000557",
        },
        {
          label: "lymphoid lineage restricted progenitor cell",
          ontology_term_id: "CL:0000838",
          children: [
            {
              label: "DN1 thymic pro-T cell",
              ontology_term_id: "CL:0000894",
            },
            {
              label: "early pro-B cell",
              ontology_term_id: "CL:0002046",
            },
            { label: "pro-B cell", ontology_term_id: "CL:0000826" },
          ],
        },
        {
          label: "macrophage dendritic cell progenitor",
          ontology_term_id: "CL:0002009",
        },
        {
          label: "megakaryocyte-erythroid progenitor cell",
          ontology_term_id: "CL:0000050",
        },
        {
          label: "monocyte",
          ontology_term_id: "CL:0000576",
          children: [
            {
              label: "CD14-low, CD16-positive monocyte",
              ontology_term_id: "CL:0002396",
            },
            {
              label: "CD14-positive monocyte",
              ontology_term_id: "CL:0001054",
            },
            {
              label: "CD14-positive, CD16-negative classical monocyte",
              ontology_term_id: "CL:0002057",
            },
            {
              label: "CD14-positive, CD16-positive monocyte",
              ontology_term_id: "CL:0002397",
            },
            {
              label: "classical monocyte",
              ontology_term_id: "CL:0000860",
            },
            {
              label: "intermediate monocyte",
              ontology_term_id: "CL:0002393",
            },
            {
              label: "non-classical monocyte",
              ontology_term_id: "CL:0000875",
            },
          ],
        },
        {
          label: "myeloid lineage restricted progenitor cell",
          ontology_term_id: "CL:0000839",
          children: [
            {
              label: "erythroid progenitor cell",
              ontology_term_id: "CL:0000038",
            },
            {
              label: "erythroid progenitor cell, mammalian",
              ontology_term_id: "CL:0001066",
            },
            {
              label: "granulocytopoietic cell",
              ontology_term_id: "CL:0002191",
            },
            { label: "promonocyte", ontology_term_id: "CL:0000559" },
          ],
        },
        {
          label: "oligodendrocyte precursor cell",
          ontology_term_id: "CL:0002453",
        },
      ],
    },
    {
      label: "renal intercalated cell",
      ontology_term_id: "CL:0005010",
    },
    { label: "renal principal cell", ontology_term_id: "CL:0005009" },
    { label: "salivary gland cell", ontology_term_id: "CL:0009005" },
    { label: "sebaceous gland cell", ontology_term_id: "CL:2000021" },
    {
      label: "secretory cell",
      ontology_term_id: "CL:0000151",
      children: [
        {
          label: "GABAergic neuron",
          ontology_term_id: "CL:0000617",
          children: [
            {
              label:
                "caudal ganglionic eminence derived GABAergic cortical interneuron",
              ontology_term_id: "CL:4023070",
            },
            {
              label: "cerebellar Golgi cell",
              ontology_term_id: "CL:0000119",
            },
            {
              label: "cerebral cortex GABAergic interneuron",
              ontology_term_id: "CL:0010011",
            },
            {
              label: "chandelier pvalb GABAergic cortical interneuron",
              ontology_term_id: "CL:4023036",
            },
            {
              label: "lamp5 GABAergic cortical interneuron",
              ontology_term_id: "CL:4023011",
            },
            {
              label: "pvalb GABAergic cortical interneuron",
              ontology_term_id: "CL:4023018",
            },
            {
              label: "sncg GABAergic cortical interneuron",
              ontology_term_id: "CL:4023015",
            },
            {
              label: "sst GABAergic cortical interneuron",
              ontology_term_id: "CL:4023017",
            },
            {
              label: "vip GABAergic cortical interneuron",
              ontology_term_id: "CL:4023016",
            },
          ],
        },
        {
          label: "PP cell",
          ontology_term_id: "CL:0000696",
          children: [
            {
              label: "pancreatic PP cell",
              ontology_term_id: "CL:0002275",
            },
          ],
        },
        { label: "Sertoli cell", ontology_term_id: "CL:0000216" },
        {
          label: "acinar cell of salivary gland",
          ontology_term_id: "CL:0002623",
        },
        { label: "chondrocyte", ontology_term_id: "CL:0000138" },
        { label: "club cell", ontology_term_id: "CL:0000158" },
        {
          label: "endocrine cell",
          ontology_term_id: "CL:0000163",
          children: [
            {
              label: "chromaffin cell",
              ontology_term_id: "CL:0000166",
            },
            {
              label: "cortical cell of adrenal gland",
              ontology_term_id: "CL:0002097",
            },
            {
              label: "enteroendocrine cell",
              ontology_term_id: "CL:0000164",
            },
            {
              label: "enteroendocrine cell of small intestine",
              ontology_term_id: "CL:0009006",
            },
            {
              label: "granulosa cell",
              ontology_term_id: "CL:0000501",
            },
            {
              label: "intestinal enteroendocrine cell",
              ontology_term_id: "CL:1001516",
            },
            {
              label: "neuroendocrine cell",
              ontology_term_id: "CL:0000165",
            },
            {
              label: "pancreatic endocrine cell",
              ontology_term_id: "CL:0008024",
            },
            { label: "theca cell", ontology_term_id: "CL:0000503" },
          ],
        },
        {
          label: "glandular epithelial cell",
          ontology_term_id: "CL:0000150",
          children: [
            { label: "acinar cell", ontology_term_id: "CL:0000622" },
            {
              label: "duodenum glandular cell",
              ontology_term_id: "CL:1001589",
            },
            {
              label: "paneth cell of colon",
              ontology_term_id: "CL:0009009",
            },
            {
              label: "paneth cell of epithelium of small intestine",
              ontology_term_id: "CL:1000343",
            },
            { label: "peptic cell", ontology_term_id: "CL:0000155" },
            {
              label: "seminal vesicle glandular cell",
              ontology_term_id: "CL:1001597",
            },
            {
              label: "thyroid follicular cell",
              ontology_term_id: "CL:0002258",
            },
          ],
        },
        {
          label: "glutamatergic neuron",
          ontology_term_id: "CL:0000679",
          children: [
            {
              label:
                "L2/3-6 intratelencephalic projecting glutamatergic cortical neuron",
              ontology_term_id: "CL:4023040",
            },
            {
              label:
                "L5 extratelencephalic projecting glutamatergic cortical neuron",
              ontology_term_id: "CL:4023041",
            },
            {
              label: "L6b glutamatergic cortical neuron",
              ontology_term_id: "CL:4023038",
            },
            {
              label: "OFF-bipolar cell",
              ontology_term_id: "CL:0000750",
            },
            {
              label: "ON-bipolar cell",
              ontology_term_id: "CL:0000749",
            },
            {
              label: "corticothalamic-projecting glutamatergic cortical neuron",
              ontology_term_id: "CL:4023013",
            },
            {
              label:
                "intratelencephalic-projecting glutamatergic cortical neuron",
              ontology_term_id: "CL:4023008",
            },
            {
              label: "near-projecting glutamatergic cortical neuron",
              ontology_term_id: "CL:4023012",
            },
            {
              label: "retinal bipolar neuron",
              ontology_term_id: "CL:0000748",
            },
            {
              label: "rod bipolar cell",
              ontology_term_id: "CL:0000751",
            },
          ],
        },
        {
          label: "kidney granular cell",
          ontology_term_id: "CL:0000648",
        },
        { label: "lung goblet cell", ontology_term_id: "CL:1000143" },
        {
          label: "lung neuroendocrine cell",
          ontology_term_id: "CL:1000223",
        },
        {
          label: "mammary alveolar cell",
          ontology_term_id: "CL:0002325",
        },
        {
          label: "mast cell",
          ontology_term_id: "CL:0000097",
          children: [
            {
              label: "mucosal type mast cell",
              ontology_term_id: "CL:0000485",
            },
          ],
        },
        {
          label: "mucus secreting cell",
          ontology_term_id: "CL:0000319",
          children: [
            {
              label: "bronchial goblet cell",
              ontology_term_id: "CL:1000312",
            },
            { label: "goblet cell", ontology_term_id: "CL:0000160" },
            {
              label: "ileal goblet cell",
              ontology_term_id: "CL:1000326",
            },
            {
              label: "intestine goblet cell",
              ontology_term_id: "CL:0019031",
            },
            {
              label: "large intestine goblet cell",
              ontology_term_id: "CL:1000320",
            },
            {
              label: "nasal mucosa goblet cell",
              ontology_term_id: "CL:0002480",
            },
            {
              label: "respiratory goblet cell",
              ontology_term_id: "CL:0002370",
            },
            {
              label: "small intestine goblet cell",
              ontology_term_id: "CL:1000495",
            },
            {
              label: "tracheal goblet cell",
              ontology_term_id: "CL:1000329",
            },
            {
              label: "tracheobronchial goblet cell",
              ontology_term_id: "CL:0019003",
            },
          ],
        },
        {
          label: "pancreatic acinar cell",
          ontology_term_id: "CL:0002064",
        },
        {
          label: "pancreatic epsilon cell",
          ontology_term_id: "CL:0005019",
        },
        { label: "platelet", ontology_term_id: "CL:0000233" },
        { label: "podocyte", ontology_term_id: "CL:0000653" },
        {
          label: "serous secreting cell",
          ontology_term_id: "CL:0000313",
          children: [
            {
              label: "serous cell of epithelium of bronchus",
              ontology_term_id: "CL:1000331",
            },
            {
              label: "serous cell of epithelium of trachea",
              ontology_term_id: "CL:1000330",
            },
            {
              label: "tracheobronchial serous cell",
              ontology_term_id: "CL:0019001",
            },
          ],
        },
        {
          label: "type A enteroendocrine cell",
          ontology_term_id: "CL:0002067",
          children: [
            {
              label: "pancreatic A cell",
              ontology_term_id: "CL:0000171",
            },
          ],
        },
        {
          label: "type B pancreatic cell",
          ontology_term_id: "CL:0000169",
        },
        {
          label: "type D enteroendocrine cell",
          ontology_term_id: "CL:0000502",
          children: [
            {
              label: "pancreatic D cell",
              ontology_term_id: "CL:0000173",
            },
          ],
        },
        {
          label: "type II pneumocyte",
          ontology_term_id: "CL:0002063",
        },
      ],
    },
    {
      label: "skeletal muscle satellite cell",
      ontology_term_id: "CL:0000594",
      children: [
        {
          label: "skeletal muscle satellite stem cell",
          ontology_term_id: "CL:0008011",
        },
      ],
    },
    { label: "sperm", ontology_term_id: "CL:0000019" },
    {
      label: "stem cell",
      ontology_term_id: "CL:0000034",
      children: [
        {
          label: "basal cell",
          ontology_term_id: "CL:0000646",
          children: [
            {
              label: "basal cell of epithelium of bronchus",
              ontology_term_id: "CL:1000349",
            },
            {
              label: "basal cell of epithelium of trachea",
              ontology_term_id: "CL:1000348",
            },
            {
              label: "epithelial cell of stratum germinativum of esophagus",
              ontology_term_id: "CL:1000447",
            },
            {
              label: "respiratory basal cell",
              ontology_term_id: "CL:0002633",
            },
          ],
        },
        {
          label: "embryonic stem cell",
          ontology_term_id: "CL:0002322",
        },
        {
          label: "hematopoietic stem cell",
          ontology_term_id: "CL:0000037",
          children: [
            {
              label: "CD34-positive, CD38-negative hematopoietic stem cell",
              ontology_term_id: "CL:0001024",
            },
            {
              label: "cord blood hematopoietic stem cell",
              ontology_term_id: "CL:2000095",
            },
          ],
        },
        { label: "hepatoblast", ontology_term_id: "CL:0005026" },
        {
          label: "intestinal crypt stem cell",
          ontology_term_id: "CL:0002250",
          children: [
            {
              label: "intestinal crypt stem cell of large intestine",
              ontology_term_id: "CL:0009016",
            },
            {
              label: "intestinal crypt stem cell of small intestine",
              ontology_term_id: "CL:0009017",
            },
          ],
        },
        {
          label: "keratinocyte stem cell",
          ontology_term_id: "CL:0002337",
        },
        {
          label: "mesenchymal stem cell",
          ontology_term_id: "CL:0000134",
          children: [
            {
              label: "mesenchymal stem cell of adipose tissue",
              ontology_term_id: "CL:0002570",
            },
          ],
        },
        {
          label: "neuroepithelial stem cell",
          ontology_term_id: "CL:0002259",
        },
        {
          label: "neuronal stem cell",
          ontology_term_id: "CL:0000047",
        },
        {
          label: "stem cell of epidermis",
          ontology_term_id: "CL:1000428",
          children: [
            {
              label: "basal cell of epidermis",
              ontology_term_id: "CL:0002187",
            },
          ],
        },
        {
          label: "vascular lymphangioblast",
          ontology_term_id: "CL:0005022",
        },
      ],
    },
    {
      label: "stratified epithelial cell",
      ontology_term_id: "CL:0000079",
      children: [{ label: "keratinocyte", ontology_term_id: "CL:0000312" }],
    },
    { label: "supporting cell", ontology_term_id: "CL:0000630" },
    {
      label: "surface ectodermal cell",
      ontology_term_id: "CL:0000114",
    },
    {
      label: "syncytiotrophoblast cell",
      ontology_term_id: "CL:0000525",
    },
    {
      label: "transit amplifying cell of colon",
      ontology_term_id: "CL:0009011",
    },
    {
      label: "transit amplifying cell of small intestine",
      ontology_term_id: "CL:0009012",
    },
    {
      label: "trophoblast giant cell",
      ontology_term_id: "CL:0002488",
    },
    {
      label: "urothelial cell",
      ontology_term_id: "CL:0000731",
      children: [
        {
          label: "bladder urothelial cell",
          ontology_term_id: "CL:1001428",
        },
      ],
    },
    { label: "valve cell", ontology_term_id: "CL:0000663" },
  ],
};
/* eslint-enable sort-keys -- disabling key order for readability. */

/**
 * Homo sapiens, Mus musculus and other organisms development stage ontology tree.
 */
/* eslint-disable sort-keys -- disabling key order for readability. */
export const DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.HsapDv]: [
    {
      label: "Prenatal (conception–birth)",
      ontology_term_id: "HsapDv:0000045",
      children: [
        {
          label: "Embryonic human (0–56 days)",
          ontology_term_id: "HsapDv:0000002",
          children: [
            {
              label: "Carnegie (CS1)",
              ontology_term_id: "HsapDv:0000003",
            },
            {
              label: "Cleavage (CS2)",
              ontology_term_id: "HsapDv:0000004",
            },
            {
              label: "Blastula (CS3–5)",
              ontology_term_id: "HsapDv:0000006",
            },
            {
              label: "Gastrula (CS6)",
              ontology_term_id: "HsapDv:0000010",
            },
            {
              label: "Neurula (CS7–8)",
              ontology_term_id: "HsapDv:0000012",
            },
            {
              label: "Organogenesis (CS9–23)",
              ontology_term_id: "HsapDv:0000015",
            },
          ],
        },
        {
          label: "Fetal (>56 days–birth)",
          ontology_term_id: "HsapDv:0000037",
        },
      ],
    },
    {
      label: "Immature (0–12 years)",
      ontology_term_id: "HsapDv:0000080",
      children: [
        {
          label: "Newborn human (0–1 month)",
          ontology_term_id: "HsapDv:0000082",
        },
        {
          label: "Infant (1–23 months)",
          ontology_term_id: "HsapDv:0000083",
        },
        {
          label: "Child (2–12 years)",
          ontology_term_id: "HsapDv:0000081",
        },
      ],
    },
    {
      label: "Mature (13+ years)",
      ontology_term_id: "HsapDv:0000204",
      children: [
        {
          label: "Adolescent (13–19 years)",
          ontology_term_id: "HsapDv:0000086",
        },
        {
          label: "Human adult (19+ years)",
          ontology_term_id: "HsapDv:0000087",
          children: [
            {
              label: "Early adulthood (19–45 years)",
              ontology_term_id: "HsapDv:0000088",
            },
            {
              label: "Late adulthood (45+ years)",
              ontology_term_id: "HsapDv:0000091",
            },
          ],
        },
      ],
    },
  ],
  [ONTOLOGY_VIEW_KEY.MmusDv]: [
    {
      label: "Prenatal",
      ontology_term_id: "MmusDv:0000042",
      children: [
        {
          label: "Embryonic mouse",
          ontology_term_id: "MmusDv:0000002",
          children: [
            {
              label: "Thelier stage 1 (TS1)",
              ontology_term_id: "MmusDv:0000003",
            },
            {
              label: "Cleavage (TS2–3)",
              ontology_term_id: "MmusDv:0000004",
            },
            {
              label: "Blastula (TS4–8)",
              ontology_term_id: "MmusDv:0000007",
            },
            {
              label: "Gastrula (TS9–10)",
              ontology_term_id: "MmusDv:0000013",
            },
            {
              label: "Thelier stage 11 (TS11)",
              ontology_term_id: "MmusDv:0000017",
            },
            {
              label: "Organogenesis (TS11–22)",
              ontology_term_id: "MmusDv:0000018",
            },
          ],
        },
        {
          label: "Fetal (TS23–26)",
          ontology_term_id: "MmusDv:0000031",
        },
      ],
    },
    {
      label: "Post-partum (Birth+)",
      ontology_term_id: "MmusDv:0000092",
      children: [
        {
          label: "Immature (0–6 weeks)",
          ontology_term_id: "MmusDv:0000043",
          children: [
            {
              label: "Thelier stage 27 (0–3 days)",
              ontology_term_id: "MmusDv:0000036",
            },
            {
              label: "Premature (3 days–6 weeks)",
              ontology_term_id: "MmusDv:0000112",
            },
          ],
        },
        {
          label: "Mature (6+ weeks)",
          ontology_term_id: "MmusDv:0000110",
          children: [
            {
              label: "Early adulthood (6 weeks–7 months)",
              ontology_term_id: "MmusDv:0000061",
            },
            {
              label: "Late adulthood (7+ months)",
              ontology_term_id: "MmusDv:0000097",
            },
          ],
        },
      ],
    },
  ],
  [ONTOLOGY_VIEW_KEY.UBERON]: [
    {
      label: "Embryo",
      ontology_term_id: "UBERON:0000068",
      children: [
        {
          label: "Zygote",
          ontology_term_id: "UBERON:0000106",
        },
        {
          label: "Cleavage",
          ontology_term_id: "UBERON:0000107",
        },
        {
          label: "Blastula",
          ontology_term_id: "UBERON:0000108",
        },
        {
          label: "Gastrula",
          ontology_term_id: "UBERON:0000109",
        },
        {
          label: "Neurula",
          ontology_term_id: "UBERON:0000110",
        },
        {
          label: "Organogenesis",
          ontology_term_id: "UBERON:0000111",
        },
        {
          label: "Late embryonic",
          ontology_term_id: "UBERON:0007220",
        },
      ],
    },
    {
      label: "Post embryonic",
      ontology_term_id: "UBERON:0000092",
      children: [
        {
          label: "Larval",
          ontology_term_id: "UBERON:0000069",
        },
        {
          label: "Pupal",
          ontology_term_id: "UBERON:0000070",
        },
        {
          label: "Nursing",
          ontology_term_id: "UBERON:0018685",
        },
        {
          label: "Fully formed",
          ontology_term_id: "UBERON:0000066",
          children: [
            {
              label: "Sexually immature",
              ontology_term_id: "UBERON:0000112",
            },
            {
              label: "Post-juvenile adult",
              ontology_term_id: "UBERON:0000113",
            },
          ],
        },
      ],
    },
  ],
};
/* eslint-enable sort-keys -- disabling key order for readability. */

/**
 * Tissues to be included for display in tissue ontology tree.
 */
/* eslint-disable sort-keys -- disabling key order for readability. */
export const TISSUE_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.UBERON]: [
    {
      label: "adipose tissue",
      ontology_term_id: "UBERON:0001013",
    },
    {
      label: "adrenal gland",
      ontology_term_id: "UBERON:0002369",
    },
    {
      label: "axilla",
      ontology_term_id: "UBERON:0009472",
    },
    {
      label: "bladder organ",
      ontology_term_id: "UBERON:0018707",
    },
    {
      label: "blood",
      ontology_term_id: "UBERON:0000178",
    },
    {
      label: "blood vascular",
      ontology_term_id: "UBERON:0004537",
    },
    {
      label: "bone marrow",
      ontology_term_id: "UBERON:0002371",
    },
    {
      label: "brain",
      ontology_term_id: "UBERON:0000955",
    },
    {
      label: "breast",
      ontology_term_id: "UBERON:0000310",
    },
    {
      label: "bronchial epithelial cell (cell culture)",
      ontology_term_id: "CL:0002328 (cell culture)",
    },
    {
      label: "central nervous system",
      ontology_term_id: "UBERON:0001017",
    },
    {
      label: "circulatory system",
      ontology_term_id: "UBERON:0001009",
    },
    {
      label: "chorionic villus",
      ontology_term_id: "UBERON:0007106",
    },
    {
      label: "colon",
      ontology_term_id: "UBERON:0001155",
    },
    {
      label: "cultured cell (cell culture)",
      ontology_term_id: "CL:0000010 (cell culture)",
    },
    {
      label: "digestive system",
      ontology_term_id: "UBERON:0001007",
    },
    {
      label: "embryo",
      ontology_term_id: "UBERON:0000922",
    },
    {
      label: "embryonic stem cell (cell culture)",
      ontology_term_id: "CL:0002322 (cell culture)",
    },
    {
      label: "endocrine gland",
      ontology_term_id: "UBERON:0002368",
    },
    {
      label: "endothelial cell (cell culture)",
      ontology_term_id: "CL:0000115 (cell culture)",
    },
    {
      label: "epithelial cell of alveolus of lung (cell culture)",
      ontology_term_id: "CL:0010003 (cell culture)",
    },
    {
      label: "epithelial cell of lung (cell culture)",
      ontology_term_id: "CL:0000082 (cell culture)",
    },
    {
      label: "esophagus",
      ontology_term_id: "UBERON:0001043",
    },
    {
      label: "exocrine gland",
      ontology_term_id: "UBERON:0002365",
    },
    {
      label: "eye",
      ontology_term_id: "UBERON:0000970",
    },
    {
      label: "fallopian tube",
      ontology_term_id: "UBERON:0003889",
    },
    {
      label: "gall bladder",
      ontology_term_id: "UBERON:0002110 ",
    },
    {
      label: "heart",
      ontology_term_id: "UBERON:0000948",
    },
    {
      label: "immune system",
      ontology_term_id: "UBERON:0002405",
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
      label: "knee",
      ontology_term_id: "UBERON:0001465",
    },
    {
      label: "lamina propria",
      ontology_term_id: "UBERON:0000030",
    },
    {
      label: "large intestine",
      ontology_term_id: "UBERON:0000059",
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
      label: "lymph node",
      ontology_term_id: "UBERON:0000029",
    },
    {
      label: "lymph vasculature",
      ontology_term_id: "UBERON:0004536",
    },
    {
      label: "mammary gland epithelial cell (cell culture)",
      ontology_term_id: "CL:0002327 (cell culture)",
    },
    {
      label: "mucosa",
      ontology_term_id: "UBERON:0000344",
    },
    {
      label: "musculature",
      ontology_term_id: "UBERON:0001015",
    },
    {
      label: "nose",
      ontology_term_id: "UBERON:0000004",
    },
    {
      label: "omentum",
      ontology_term_id: "UBERON:0003688 ",
    },
    {
      label: "ovary",
      ontology_term_id: "UBERON:0000992",
    },
    {
      label: "pancreas",
      ontology_term_id: "UBERON:0001264",
    },
    {
      label: "peripheral nervous system",
      ontology_term_id: "UBERON:0000010",
    },
    {
      label: "placenta",
      ontology_term_id: "UBERON:0001987",
    },
    {
      label: "pleura",
      ontology_term_id: "UBERON:0000977",
    },
    {
      label: "pleural fluid",
      ontology_term_id: "UBERON:0001087",
    },
    {
      label: "prostate",
      ontology_term_id: "UBERON:0002367",
    },
    {
      label: "renal system",
      ontology_term_id: "UBERON:0001008",
    },
    {
      label: "reproductive system",
      ontology_term_id: "UBERON:0000990",
    },
    {
      label: "respiratory system",
      ontology_term_id: "UBERON:0001004",
    },
    {
      label: "rib",
      ontology_term_id: "UBERON:0002228",
    },
    {
      label: "saliva",
      ontology_term_id: "UBERON:0001836",
    },
    {
      label: "sensory system",
      ontology_term_id: "UBERON:0001032",
    },

    {
      label: "skeletal system",
      ontology_term_id: "UBERON:0001434",
    },
    {
      label: "skin",
      ontology_term_id: "UBERON:0002097",
    },
    {
      label: "skull",
      ontology_term_id: "UBERON:0003129",
    },
    {
      label: "small intestine",
      ontology_term_id: "UBERON:0002108",
    },
    {
      label: "spinal cord",
      ontology_term_id: "UBERON:0002240",
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
      label: "testis",
      ontology_term_id: "UBERON:0000473",
    },
    {
      label: "tongue",
      ontology_term_id: " UBERON:0001723",
    },
    {
      label: "thymus",
      ontology_term_id: "UBERON:0002370",
    },
    {
      label: "ureter",
      ontology_term_id: "UBERON:0000056",
    },
    {
      label: "urethra",
      ontology_term_id: "UBERON:0000057",
    },
    {
      label: "urinary bladder",
      ontology_term_id: "UBERON:0001255",
    },
    {
      label: "uterus",
      ontology_term_id: "UBERON:0000995",
    },
    {
      label: "uterus (organoid)",
      ontology_term_id: "UBERON:0000995 (organoid)",
    },
    {
      label: "vasculature",
      ontology_term_id: "UBERON:0002049",
    },
  ],
};
/* eslint-enable sort-keys -- disabling key order for readability. */

/**
 * Configuration for each category.
 */
const CATEGORY_CONFIGS: (CategoryConfig | OntologyCategoryConfig)[] = [
  {
    analyticsEvent: EVENTS.FILTER_SELECT_ASSAY,
    categoryKey: CATEGORY_KEY.ASSAY,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.ASSAY,
    filterKey: "assay",
    label: "Assay",
    multiselect: true,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_CELL_COUNT,
    categoryKey: CATEGORY_KEY.CELL_COUNT,
    categoryType: CATEGORY_FILTER_TYPE.BETWEEN,
    filterCategoryKey: FILTER_CATEGORY_KEY.CELL_COUNT,
    filterKey: "cell_count",
    label: "Cell Count",
    multiselect: false,
  },
  {
    // TODO(cc) possibly remove with #2569.
    analyticsEvent: EVENTS.FILTER_SELECT_CELL_TYPE,
    categoryKey: CATEGORY_KEY.CELL_TYPE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.CELL_TYPE_DEPRECATED,
    filterKey: "cell_type",
    label: "Cell Type",
    multiselect: true,
  },
  {
    // TODO(cc) add analytics event with #2569.
    categoryKey: CATEGORY_KEY.CELL_TYPE_ANCESTORS,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.CELL_TYPE,
    filterKey: "cell_type_ancestors",
    isLabelVisible: false,
    isSearchable: true,
    isZerosVisible: true,
    label: "Cell Type (Ontology)",
    multiselect: true,
    ontologyTermSet: CELL_TYPE_ONTOLOGY_TERM_SET,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_DEVELOPMENT_STAGE,
    categoryKey: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.DEVELOPMENT_STAGE,
    filterKey: "development_stage_ancestors",
    isLabelVisible: true,
    isSearchable: false,
    isZerosVisible: true,
    label: "Development Stage",
    multiselect: true,
    ontologyTermSet: DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_DISEASE,
    categoryKey: CATEGORY_KEY.DISEASE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.DISEASE,
    filterKey: "disease",
    label: "Disease",
    multiselect: true,
    pinnedCategoryValues: [CATEGORY_VALUE_KEY.NORMAL],
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_ETHNICITY,
    categoryKey: CATEGORY_KEY.ETHNICITY,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.ETHNICITY,
    filterKey: "ethnicity",
    label: "Ethnicity",
    multiselect: true,
    tooltip:
      "Ethnicity only applies to Homo sapiens which is not selected in the Organism filter.",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_GENE_COUNT,
    categoryKey: CATEGORY_KEY.MEAN_GENES_PER_CELL,
    categoryType: CATEGORY_FILTER_TYPE.BETWEEN,
    filterCategoryKey: FILTER_CATEGORY_KEY.GENE_COUNT,
    filterKey: "mean_genes_per_cell",
    label: "Gene Count",
    multiselect: false,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_ORGANISM,
    categoryKey: CATEGORY_KEY.ORGANISM,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.ORGANISM,
    filterKey: "organism",
    label: "Organism",
    multiselect: true,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_AUTHORS,
    categoryKey: CATEGORY_KEY.PUBLICATION_AUTHORS,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.PUBLICATION_AUTHORS,
    filterKey: "publicationAuthors",
    label: "Author",
    multiselect: true,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_PUBLICATION_DATE,
    categoryKey: CATEGORY_KEY.PUBLICATION_DATE_VALUES,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.PUBLICATION_DATE_VALUES,
    filterKey: "publicationDateValues",
    label: "Publication Date",
    multiselect: false,
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_SEX,
    categoryKey: CATEGORY_KEY.SEX,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.SEX,
    filterKey: "sex",
    label: "Sex",
    multiselect: true,
  },
  {
    // TODO(cc) possibly remove with #2569.
    analyticsEvent: EVENTS.FILTER_SELECT_TISSUE,
    categoryKey: CATEGORY_KEY.TISSUE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.TISSUE_DEPRECATED,
    filterKey: "tissue",
    label: "Tissue",
    multiselect: true,
  },
  {
    // TODO(cc) add analytics event with #2569.
    categoryKey: CATEGORY_KEY.TISSUE_ANCESTORS,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    filterCategoryKey: FILTER_CATEGORY_KEY.TISSUE,
    filterKey: "tissue_ancestors",
    isLabelVisible: false,
    isSearchable: true,
    isZerosVisible: false,
    label: "Tissue (Ontology)",
    multiselect: true,
    ontologyTermSet: TISSUE_ONTOLOGY_TERM_SET,
  },
];

/**
 * Category configs keyed by filter category key, for convenience. Using object literal with type
 * KeyedFilterCategoryConfig rather than generic Map to prevent having to null check values.
 */
export const CATEGORY_CONFIGS_BY_FILTER_CATEGORY_KEY: KeyedFilterCategoryConfigs =
  CATEGORY_CONFIGS.reduce(
    (accum: KeyedFilterCategoryConfigs, config: CategoryConfig) => {
      return {
        ...accum,
        [config.filterCategoryKey]: config,
      };
    },
    {} as KeyedFilterCategoryConfigs
  );

/**
 * "value" prop passed to react-table's Cell function
 */
export type CellPropsValue<T> = { value: CellValue<T> };

/**
 * View model of category.
 */
export type CategoryView =
  | OntologyCategoryView
  | RangeCategoryView
  | SelectCategoryView;

/**
 * Category values to be used as keys. For example, "Homo sapiens" or "10X 3' v2 sequencing".
 */
export type CategoryValueKey = string;

/**
 * Join of dataset and collection information, optimized for filtering datasets.
 */
export interface DatasetRow extends Categories, PublisherMetadataCategories {
  cell_count: number | null;
  collection_id: Collection["id"];
  collection_name: Collection["name"];
  explorer_url: string;
  id: string;
  isOverMaxCellCount: boolean;
  mean_genes_per_cell: number | null;
  name: string;
  published_at: number;
  recency: number; // Used by sort
  revised_at?: number;
}

/**
 * Display values of unspecified ethnicity labels.
 */
export enum ETHNICITY_UNSPECIFIED_LABEL {
  "unknown" = "Unknown",
}

/**
 * List of ethnicity ontology labels to exclude from filter functionality.
 */
export const ETHNICITY_DENY_LIST = ["na"];

/**
 * Model of filter category configs keyed by category key. Used instead of generic Map to prevent null checking when
 * grabbing keyed value.
 */
type KeyedFilterCategoryConfigs = {
  [K in FILTER_CATEGORY_KEY]: CategoryConfig;
};

/**
 * Possible set of organism values.
 */
export enum ORGANISM {
  "HOMO_SAPIENS" = "Homo sapiens",
  "MUS_MUSCULUS" = "Mus musculus",
}

/**
 * Function invoked when selected state of a category value is toggled or range is selected.
 */
export type OnFilterFn = (
  categoryKey: FilterCategoryKey, // TODO(cc) rename categoryKey here and throughout.
  selectedValue: CategoryValueKey | Range
) => void;

/**
 * Function invoked when filter category input value is changed.
 */
export type OnUpdateSearchValueFn = (
  changeEvent: ChangeEvent<HTMLInputElement>,
  setSearchValue: SetSearchValueFn
) => void;

/**
 * Min and max values selected in range category. Empty array if no range is specified (e.g. on clear of range).
 */
export type Range = [number, number] | [];

/**
 * Tree view model of ontology view. For example, development stage has three tree views (human, mouse and other)
 * whereas tissue has one tree view.
 */
export interface OntologyCategoryTreeView {
  children: OntologyCategoryTreeNodeView[];
  label?: string;
  selectedViews: OntologyCategoryTreeNodeView[];
}

/**
 * View model of a node in an ontology tree view including partial selected state and possibly any children of this
 * node.
 */
export interface OntologyCategoryTreeNodeView extends SelectCategoryValueView {
  children?: OntologyCategoryTreeNodeView[];
  selectedPartial: boolean; // True if value is a parent node and some children node are selected.
}

/**
 * View model of ontology category.
 */
export interface OntologyCategoryView {
  isDisabled?: boolean;
  isSearchable: boolean;
  isZerosVisible: boolean;
  key: FilterCategoryKey;
  label: string;
  views: OntologyCategoryTreeView[];
  tooltip?: string;
}

/**
 * Node of ontology tree including ontology ID, label and any possible children.
 */
export interface OntologyNode extends Ontology {
  children?: OntologyNode[];
}

/**
 * Ontology tree structures, keyed by view key. This is the allowed set of ontology values, configured per category.
 */
export type OntologyTermSet = { [K in ONTOLOGY_VIEW_KEY]?: OntologyNode[] };

/**
 * View model of range metadata key.
 */
export interface RangeCategoryView {
  isDisabled?: boolean;
  key: FilterCategoryKey;
  label: string;
  max: number;
  min: number;
  selectedMax?: number;
  selectedMin?: number;
  tooltip?: string;
}

/**
 * Possible set of values that publication dates can be binned into.
 */
export const PUBLICATION_DATE_VALUES: number[] = [1, 3, 6, 12, 24, 36];

/**
 * Display values of publicationDateValue labels. Enum must be non-numeric key for reverse lookup.
 */
export enum PUBLICATION_DATE_LABELS {
  "LABEL_1" = "1 month",
  "LABEL_12" = "12 months",
  "LABEL_24" = "2 years",
  "LABEL_3" = "3 months",
  "LABEL_36" = "3 years",
  "LABEL_6" = "6 months",
}

/**
 * Publication-related filterable values of collections and datasets.
 */
export interface PublisherMetadataCategories {
  publicationAuthors?: string[];
  publicationDateValues?: number[]; // Set of date bins that publication date falls within
}

/**
 * "row" prop passed to react-table's Cell function.
 */
export type RowPropsValue<T extends Categories> = { row: Row<T> };

/**
 * View model of metadata value, selected state and count for single or multiselect categories.
 */
export interface SelectCategoryValueView {
  key: CategoryValueKey;
  count: number;
  label: string;
  selected: boolean;
}

/**
 * Metadata values grouped by metadata key, for single or multiselect categories.
 */
export interface SelectCategoryView {
  isDisabled?: boolean;
  key: FilterCategoryKey;
  label: string;
  pinnedValues: SelectCategoryValueView[];
  tooltip?: string;
  unpinnedValues: SelectCategoryValueView[];
  values: SelectCategoryValueView[]; // both pinned and unpinned values
}

/**
 * Function invoked to update state for the filter category input value.
 */
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;
