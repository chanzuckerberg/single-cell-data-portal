import { EVENTS } from "src/common/analytics/events";
import {
  ANALYTICS_PAYLOAD_KEY,
  CATEGORY_FILTER_ID,
  CATEGORY_FILTER_PANEL_ID,
  CATEGORY_VALUE_KEY,
  CategoryFilterConfig,
  KeyedCategoryFilterConfigs,
  ONTOLOGY_VIEW_KEY,
  OntologyDescendants,
  OntologyTermSet,
  PINNED_POSITION,
} from "src/components/common/Filter/common/entities";
import cell_type_descendants_json from "src/components/common/Filter/descendant_mappings/cell_type_descendants.json";
import tissue_descendants_json from "src/components/common/Filter/descendant_mappings/tissue_descendants.json";

/**
 * Homo sapiens, Mus musculus and other organisms development stage ontology tree.
 */
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

/**
 * List of ethnicity ontology labels to exclude from filter functionality.
 */
export const SELF_REPORTED_ETHNICITY_DENY_LIST = ["na"];

/**
 * List of suspension types to exclude from filter functionality.
 */
export const SUSPENSION_TYPE_DENY_LIST = ["na"];

/**
 * String value to append to labels in multi-panel categories if the value appears in more than one panel.
 */
export const LABEL_SUFFIX_NON_SPECIFIC = ", non-specific";

/**
 * Possible set of values that publication dates can be binned into.
 */
export const PUBLICATION_DATE_VALUES: number[] = [1, 3, 6, 12, 24, 36];

/**
 * Cell types to be included for display in cell type cell class ontology tree.
 */
/* eslint-disable sonarjs/no-duplicate-string -- disable no dupes, this value is hand-curated */
export const CELL_TYPE_CELL_CLASS_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.UBERON]: [
    {
      label: "cardiocyte",
      ontology_term_id: "CL:0002494",
    },
    {
      label: "connective tissue cell",
      ontology_term_id: "CL:0002320",
    },
    {
      label: "defensive cell",
      ontology_term_id: "CL:0000473",
    },
    {
      label: "epithelial cell",
      ontology_term_id: "CL:0000066",
    },
    {
      label: "hematopoietic cell",
      ontology_term_id: "CL:0000988",
    },
    {
      label: "neural cell",
      ontology_term_id: "CL:0002319",
    },
    {
      label: "precursor cell",
      ontology_term_id: "CL:0011115",
    },
    {
      label: "secretory cell",
      ontology_term_id: "CL:0000151",
    },
    {
      label: "germ cell line",
      ontology_term_id: "CL:0000039",
    },
    {
      label: "ciliated cell",
      ontology_term_id: "CL:0000064",
    },
    {
      label: "contractile cell",
      ontology_term_id: "CL:0000183",
    },
    {
      label: "cell of skeletal muscle",
      ontology_term_id: "CL:0000188",
    },
    {
      label: "motile cell",
      ontology_term_id: "CL:0000219",
    },
    {
      label: "stuff accumulating cell",
      ontology_term_id: "CL:0000325",
    },
    {
      label: "extraembryonic cell",
      ontology_term_id: "CL:0000349",
    },
    {
      label: "germ cell",
      ontology_term_id: "CL:0000586",
    },
    {
      label: "supporting cell",
      ontology_term_id: "CL:0000630",
    },
    {
      label: "bone cell",
      ontology_term_id: "CL:0001035",
    },
    {
      label: "abnormal cell",
      ontology_term_id: "CL:0001061",
    },
    {
      label: "embryonic cell (metazoa)",
      ontology_term_id: "CL:0002321",
    },
    {
      label: "transit amplifying cell",
      ontology_term_id: "CL:0009010",
    },
    {
      label: "lower urinary tract cell",
      ontology_term_id: "CL:1000600",
    },
    {
      label: "perivascular cell",
      ontology_term_id: "CL:4033054",
    },
  ],
};
/* eslint-enable sonarjs/no-duplicate-string -- disable no dupes, this value is hand-curated */

/**
 * Cell types to be included for display in cell type level 1 ontology tree.
 */
/* eslint-disable sonarjs/no-duplicate-string -- disable no dupes, this value is hand-curated */
export const CELL_TYPE_LEVEL_2_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.UBERON]: [
    {
      label: "CD4-positive, alpha-beta T cell",
      ontology_term_id: "CL:0000624",
    },
    {
      label: "CD8-positive, alpha-beta T cell",
      ontology_term_id: "CL:0000625",
    },
    {
      label: "T cell",
      ontology_term_id: "CL:0000084",
    },
    {
      label: "B cell",
      ontology_term_id: "CL:0000236",
    },
    {
      label: "dendritic cell",
      ontology_term_id: "CL:0000451",
    },
    {
      label: "monocyte",
      ontology_term_id: "CL:0000576",
    },
    {
      label: "macrophage",
      ontology_term_id: "CL:0000235",
    },
    {
      label: "lymphocyte",
      ontology_term_id: "CL:0000542",
    },
    {
      label: "leukocyte",
      ontology_term_id: "CL:0000738",
    },
    {
      label: "myeloid cell",
      ontology_term_id: "CL:0000763",
    },
    {
      label: "hematopoietic precursor cell",
      ontology_term_id: "CL:0008001",
    },
    {
      label: "phagocyte",
      ontology_term_id: "CL:0000234",
    },
    {
      label: "glutamatergic neuron",
      ontology_term_id: "CL:0000679",
    },
    {
      label: "GABAergic neuron",
      ontology_term_id: "CL:0000617",
    },
    {
      label: "interneuron",
      ontology_term_id: "CL:0000099",
    },
    {
      label: "glial cell",
      ontology_term_id: "CL:0000125",
    },
    {
      label: "sensory neuron",
      ontology_term_id: "CL:0000101",
    },
    {
      label: "motor neuron",
      ontology_term_id: "CL:0000100",
    },
    {
      label: "CNS neuron (sensu Vertebrata)",
      ontology_term_id: "CL:0000117",
    },
    {
      label: "neuron",
      ontology_term_id: "CL:0000540",
    },
    {
      label: "pericyte",
      ontology_term_id: "CL:0000669",
    },
    {
      label: "stromal cell",
      ontology_term_id: "CL:0000499",
    },
    {
      label: "fibroblast",
      ontology_term_id: "CL:0000057",
    },
    {
      label: "exocrine cell",
      ontology_term_id: "CL:0000152",
    },
    {
      label: "endocrine cell",
      ontology_term_id: "CL:0000163",
    },
    {
      label: "endothelial cell",
      ontology_term_id: "CL:0000115",
    },
    {
      label: "endo-epithelial cell",
      ontology_term_id: "CL:0002076",
    },
    {
      label: "meso-epithelial cell",
      ontology_term_id: "CL:0002078",
    },
    {
      label: "progenitor cell",
      ontology_term_id: "CL:0011026",
    },
    {
      label: "male germ cell",
      ontology_term_id: "CL:0000015",
    },
    {
      label: "female germ cell",
      ontology_term_id: "CL:0000021",
    },
    {
      label: "stem cell",
      ontology_term_id: "CL:0000034",
    },
    {
      label: "non-terminally differentiated cell",
      ontology_term_id: "CL:0000055",
    },
    {
      label: "duct epithelial cell",
      ontology_term_id: "CL:0000068",
    },
    {
      label: "columnar/cuboidal epithelial cell",
      ontology_term_id: "CL:0000075",
    },
    {
      label: "squamous epithelial cell",
      ontology_term_id: "CL:0000076",
    },
    {
      label: "stratified epithelial cell",
      ontology_term_id: "CL:0000079",
    },
    {
      label: "epithelial cell of lung",
      ontology_term_id: "CL:0000082",
    },
    {
      label: "epithelial cell of pancreas",
      ontology_term_id: "CL:0000083",
    },
    {
      label: "neuron associated cell",
      ontology_term_id: "CL:0000095",
    },
    {
      label: "sensory epithelial cell",
      ontology_term_id: "CL:0000098",
    },
    {
      label: "fat cell",
      ontology_term_id: "CL:0000136",
    },
    {
      label: "pigment cell",
      ontology_term_id: "CL:0000147",
    },
    {
      label: "glandular epithelial cell",
      ontology_term_id: "CL:0000150",
    },
    {
      label: "seromucus secreting cell",
      ontology_term_id: "CL:0000159",
    },
    {
      label: "hepatocyte",
      ontology_term_id: "CL:0000182",
    },
    {
      label: "myofibroblast cell",
      ontology_term_id: "CL:0000186",
    },
    {
      label: "muscle cell",
      ontology_term_id: "CL:0000187",
    },
    {
      label: "ectodermal cell",
      ontology_term_id: "CL:0000221",
    },
    {
      label: "mesodermal cell",
      ontology_term_id: "CL:0000222",
    },
    {
      label: "urothelial cell",
      ontology_term_id: "CL:0000244",
    },
    {
      label: "trophoblast cell",
      ontology_term_id: "CL:0000351",
    },
    {
      label: "enterocyte",
      ontology_term_id: "CL:0000584",
    },
    {
      label: "germ cell",
      ontology_term_id: "CL:0000586",
    },
    {
      label: "primordial germ cell",
      ontology_term_id: "CL:0000670",
    },
    {
      label: "muscle precursor cell",
      ontology_term_id: "CL:0000680",
    },
    {
      label: "neoplastic cell",
      ontology_term_id: "CL:0001063",
    },
    {
      label: "ecto-epithelial cell",
      ontology_term_id: "CL:0002077",
    },
    {
      label: "vertebrate lens cell",
      ontology_term_id: "CL:0002222",
    },
    {
      label: "mammary gland epithelial cell",
     o ontology_term_id: "CL:0002327",
    },
    {
      label: "adventitial cell",
      ontology_term_id: "CL:0002503",
    },
    {
      label: "kidney epithelial cell",
      ontology_term_id: "CL:0002518",
    },
    {
      label: "epithelial cell of cervix",
      ontology_term_id: "CL:0002535",
    },
    {
      label: "epithelial cell of amnion",
      ontology_term_id: "CL:0002536",
    },
    {
      label: "ionocyte",
      ontology_term_id: "CL:0005006",
    },
    {
      label: "mesenchymal cell",
      ontology_term_id: "CL:0008019",
    },
    {
      label: "mural cell",
      ontology_term_id: "CL:0008034",
    },
    {
      label: "transit amplifying cell",
      ontology_term_id: "CL:0009010",
    },
    {
      label: "epithelial cell of urethra",
      ontology_term_id: "CL:1000296",
    },
    {
      label: "kidney cell",
      ontology_term_id: "CL:1000497",
    },
    {
      label: "pituitary gland cell",
      ontology_term_id: "CL:2000004",
    },
    {
      label: "ovarian surface epithelial cell",
      ontology_term_id: "CL:2000064",
    },
    {
      label: "interstitial cell",
      ontology_term_id: "CL:4030031",
    },
      ],
};
/* eslint-enable sonarjs/no-duplicate-string -- disable no dupes, this value is hand-curated */

export const CELL_TYPE_DESCENDANTS: OntologyDescendants =
  cell_type_descendants_json;
export const TISSUE_DESCENDANTS: OntologyDescendants = tissue_descendants_json;

/**
 * Tissues to be included for display in tissue system ontology tree.
 */
/* eslint-disable sonarjs/no-duplicate-string -- disable no dupes, this value is hand-curated */
export const TISSUE_ORGAN_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.UBERON]: [
    {
      label: "adipose tissue",
      ontology_term_id: "UBERON:0001013",
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
      label: "esophagus",
      ontology_term_id: "UBERON:0001043",
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
      ontology_term_id: "UBERON:0002110",
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
      label: "lymph node",
      ontology_term_id: "UBERON:0000029",
    },
    {
      label: "nose",
      ontology_term_id: "UBERON:0000004",
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
      label: "placenta",
      ontology_term_id: "UBERON:0001987",
    },
    {
      label: "skin of body",
      ontology_term_id: "UBERON:0002097",
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
      label: "thymus",
      ontology_term_id: "UBERON:0002370",
    },
    {
      label: "thyroid gland",
      ontology_term_id: "UBERON:0002046",
    },
    {
      label: "tongue",
      ontology_term_id: "UBERON:0001723",
    },
    {
      label: "uterus",
      ontology_term_id: "UBERON:0000995",
    },
  ],
};
/* eslint-enable sonarjs/no-duplicate-string -- disable no dupes, this value is hand-curated */

/**
 * Tissues to be included for display in tissue system ontology tree.
 */
export const TISSUE_SYSTEM_ONTOLOGY_TERM_SET: OntologyTermSet = {
  [ONTOLOGY_VIEW_KEY.UBERON]: [
    {
      label: "central nervous system",
      ontology_term_id: "UBERON:0001017",
    },
    {
      label: "cardiovascular system",
      ontology_term_id: "UBERON:0004535",
    },
    {
      label: "circulatory system",
      ontology_term_id: "UBERON:0001009",
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
      label: "endocrine system",
      ontology_term_id: "UBERON:0000949",
    },
    {
      label: "exocrine system",
      ontology_term_id: "UBERON:0002330",
    },
    {
      label: "hematopoietic system",
      ontology_term_id: "UBERON:0002390",
    },
    {
      label: "immune system",
      ontology_term_id: "UBERON:0002405",
    },
    {
      label: "musculature of body",
      ontology_term_id: "UBERON:0000383",
    },
    {
      label: "nervous system",
      ontology_term_id: "UBERON:0001016",
    },
    {
      label: "peripheral nervous system",
      ontology_term_id: "UBERON:0000010",
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
      label: "sensory system",
      ontology_term_id: "UBERON:0001032",
    },
    {
      label: "skeletal system",
      ontology_term_id: "UBERON:0001434",
    },
  ],
};

/**
 * Generic tooltip displayed when select, ontology or range category is disabled due to no values matching current
 * filter.
 */
export const TOOLTIP_CATEGORY_DISABLED =
  "There are no values that meet the current filter criteria";

/**
 * Configuration of each category filter.
 */
const CATEGORY_FILTER_CONFIGS: CategoryFilterConfig[] = [
  {
    analyticsEvent: EVENTS.FILTER_SELECT_ASSAY,
    categoryFilterId: CATEGORY_FILTER_ID.ASSAY,
    filterOnKey: "assay",
    label: "Assay",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_CELL_COUNT,
    categoryFilterId: CATEGORY_FILTER_ID.CELL_COUNT,
    filterOnKey: "cell_count",
    label: "Cell Count",
    labelKind: "VALUE",
    matchKind: "BETWEEN",
    multiselect: false,
    valueSourceKind: "NONE",
    viewKind: "RANGE",
  },
  {
    categoryFilterId: CATEGORY_FILTER_ID.CELL_TYPE_CALCULATED,
    descendants: CELL_TYPE_DESCENDANTS,
    filterOnKey: "cellTypeCalculated",
    label: "Cell Type",
    labelKind: "LOOKUP_LABEL_BY_TERM_ID",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    panels: [
      {
        analyticsEvent: EVENTS.FILTER_SELECT_CELL_CLASS,
        analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.CELL_CLASS,
        filterValueKind: "INFERRED_EXPLICIT",
        id: CATEGORY_FILTER_PANEL_ID.CELL_TYPE_CELL_CLASS,
        label: "Cell Class",
        searchKind: "SEARCH_SINGLE_SELECT",
        source: CELL_TYPE_CELL_CLASS_ONTOLOGY_TERM_SET,
        sourceKind: "CURATED",
      },
      {
        analyticsEvent: EVENTS.FILTER_SELECT_CELL_SUBCLASS,
        analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.CELL_SUBCLASS,
        filterValueKind: "INFERRED_EXPLICIT",
        id: CATEGORY_FILTER_PANEL_ID.CELL_TYPE_CELL_SUBCLASS,
        label: "Cell Subclass",
        searchKind: "SEARCH_SINGLE_SELECT",
        source: CELL_TYPE_LEVEL_2_ONTOLOGY_TERM_SET,
        sourceKind: "CURATED",
      },
      {
        analyticsEvent: EVENTS.FILTER_SELECT_CELL_TYPE,
        analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.CELL_TYPE,
        filterValueKind: "EXPLICIT_ONLY",
        id: CATEGORY_FILTER_PANEL_ID.CELL_TYPE,
        label: "Cell Type",
        searchKind: "SEARCH_MULTI_SELECT",
        sourceKind: "EXPLICIT_ONLY",
      },
    ],
    valueSourceKind: "NONE",
    viewKind: "MULTI_PANEL",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_CONSORTIA,
    analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.CONSORTIA,
    categoryFilterId: CATEGORY_FILTER_ID.CONSORTIA,
    filterOnKey: "consortia",
    label: "Consortia",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    pinnedCategoryValues: [CATEGORY_VALUE_KEY.NO_CONSORTIUM],
    pinnedPosition: PINNED_POSITION.BOTTOM,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    // analyticsEvent: EVENTS.FILTER_SELECT_CURATOR,
    categoryFilterId: CATEGORY_FILTER_ID.CURATOR_NAME,
    filterOnKey: "curator_name",
    label: "Curator",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_DEVELOPMENT_STAGE,
    categoryFilterId: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
    filterOnKey: "development_stage_ancestors",
    isLabelVisible: true,
    isSearchable: false,
    isZerosVisible: true,
    label: "Development Stage",
    labelKind: "LOOKUP_LABEL_BY_TERM_ID",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    source: DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET,
    valueSourceKind: "CURATED",
    viewKind: "CURATED_ONTOLOGY",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_DISEASE,
    categoryFilterId: CATEGORY_FILTER_ID.DISEASE,
    filterOnKey: "disease",
    label: "Disease",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    pinnedCategoryValues: [CATEGORY_VALUE_KEY.NORMAL],
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_SELF_REPORTED_ETHNICITY,
    categoryFilterId: CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY,
    filterOnKey: "self_reported_ethnicity",
    label: "Self-Reported Ethnicity",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    tooltip:
      "Self-Reported Ethnicity only applies to Homo sapiens which is not selected in the Organism filter.",
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_GENE_COUNT,
    categoryFilterId: CATEGORY_FILTER_ID.GENE_COUNT,
    filterOnKey: "mean_genes_per_cell",
    label: "Gene Count",
    labelKind: "VALUE",
    matchKind: "BETWEEN",
    multiselect: false,
    valueSourceKind: "NONE",
    viewKind: "RANGE",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_ORGANISM,
    categoryFilterId: CATEGORY_FILTER_ID.ORGANISM,
    filterOnKey: "organism",
    label: "Organism",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_PUBLICATION,
    categoryFilterId: CATEGORY_FILTER_ID.PUBLICATION,
    filterOnKey: "summaryCitation",
    label: "Publication",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    pinnedCategoryValues: [CATEGORY_VALUE_KEY.NO_PUBLICATION],
    pinnedPosition: PINNED_POSITION.BOTTOM,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_PUBLICATION_DATE,
    categoryFilterId: CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES,
    filterOnKey: "publicationDateValues",
    label: "Publication Date",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    // analyticsEvent: EVENTS.FILTER_SELECT_STATUS,
    categoryFilterId: CATEGORY_FILTER_ID.STATUS,
    filterOnKey: "status",
    label: "Status",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_SEX,
    categoryFilterId: CATEGORY_FILTER_ID.SEX,
    filterOnKey: "sex",
    label: "Sex",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    analyticsEvent: EVENTS.FILTER_SELECT_SUSPENSION_TYPE,
    analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.SUSPENSION_TYPE,
    categoryFilterId: CATEGORY_FILTER_ID.SUSPENSION_TYPE,
    filterOnKey: "suspension_type",
    label: "Suspension Type",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    valueSourceKind: "NONE",
    viewKind: "SELECT",
  },
  {
    categoryFilterId: CATEGORY_FILTER_ID.TISSUE_CALCULATED,
    descendants: TISSUE_DESCENDANTS,
    filterOnKey: "tissueCalculated",
    label: "Tissue",
    labelKind: "LOOKUP_LABEL_BY_TERM_ID",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    panels: [
      {
        analyticsEvent: EVENTS.FILTER_SELECT_SYSTEM,
        analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.SYSTEM,
        filterValueKind: "INFERRED_EXPLICIT",
        id: CATEGORY_FILTER_PANEL_ID.TISSUE_SYSTEM,
        label: "System",
        searchKind: "SEARCH_SINGLE_SELECT",
        source: TISSUE_SYSTEM_ONTOLOGY_TERM_SET,
        sourceKind: "CURATED",
      },
      {
        analyticsEvent: EVENTS.FILTER_SELECT_ORGAN,
        analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.ORGAN,
        filterValueKind: "INFERRED_EXPLICIT",
        id: CATEGORY_FILTER_PANEL_ID.TISSUE_ORGAN,
        label: "Organ",
        searchKind: "SEARCH_SINGLE_SELECT",
        source: TISSUE_ORGAN_ONTOLOGY_TERM_SET,
        sourceKind: "CURATED",
      },
      {
        analyticsEvent: EVENTS.FILTER_SELECT_TISSUE,
        analyticsPayloadKey: ANALYTICS_PAYLOAD_KEY.TISSUE,
        filterValueKind: "EXPLICIT_ONLY",
        id: CATEGORY_FILTER_PANEL_ID.TISSUE,
        label: "Tissue",
        searchKind: "SEARCH_MULTI_SELECT",
        sourceKind: "EXPLICIT_ONLY",
      },
    ],
    valueSourceKind: "NONE",
    viewKind: "MULTI_PANEL",
  },
];

/**
 * Category filter configs keyed by category filter ID, for convenience. Using object literal with type
 * KeyedCategoryFilterConfigs rather than generic Map to prevent having to null check values.
 */
export const CATEGORY_FILTER_CONFIGS_BY_ID: KeyedCategoryFilterConfigs =
  CATEGORY_FILTER_CONFIGS.reduce(
    (accum: KeyedCategoryFilterConfigs, config: CategoryFilterConfig) => {
      return {
        ...accum,
        [config.categoryFilterId]: config,
      };
    },
    {} as KeyedCategoryFilterConfigs
  );

/**
 * Case insensitive sorter
 */
export const COLLATOR_CASE_INSENSITIVE = new Intl.Collator("en", {
  numeric: true,
  sensitivity: "base",
});
