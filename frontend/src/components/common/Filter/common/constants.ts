import { EVENTS } from "src/common/analytics/events";
import {
  ANALYTICS_PAYLOAD_KEY,
  CategoryFilterConfig,
  CATEGORY_FILTER_ID,
  CATEGORY_FILTER_PANEL_ID,
  CATEGORY_VALUE_KEY,
  KeyedCategoryFilterConfigs,
  OntologyDescendants,
  OntologyTermSet,
  ONTOLOGY_VIEW_KEY,
} from "src/components/common/Filter/common/entities";

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
 * List of ethnicity ontology labels to exclude from filter functionality.
 */
export const ETHNICITY_DENY_LIST = ["na"];

/**
 * String value to append to labels in multi-panel categories if the value appears in more than one panel.
 */
export const LABEL_SUFFIX_NON_SPECIFIC = ", non-specific";

/**
 * Possible set of values that publication dates can be binned into.
 */
export const PUBLICATION_DATE_VALUES: number[] = [1, 3, 6, 12, 24, 36];

/* eslint-disable sort-keys -- disable key order, contents of this file are generated */
/* eslint-disable sonarjs/no-duplicate-string -- disable no dupes, this value is generated */
export const TISSUE_DESCENDANTS: OntologyDescendants = {
  "UBERON:0001017": [
    "UBERON:0002811",
    "UBERON:0034751",
    "UBERON:0009835",
    "UBERON:0002802",
    "UBERON:0001872",
    "UBERON:0000956",
    "UBERON:0001871",
    "UBERON:0004725",
    "UBERON:0002808",
    "UBERON:0008933",
    "UBERON:0002809",
    "UBERON:0002421",
    "UBERON:0001882",
    "UBERON:0008934",
    "UBERON:0002728",
    "UBERON:0002810",
    "UBERON:0001890",
    "UBERON:0002771",
    "UBERON:0002894",
    "UBERON:0001384",
    "UBERON:0016540",
    "UBERON:0023852",
    "UBERON:0000451",
    "UBERON:0004167",
    "UBERON:0003876",
    "UBERON:0009834",
    "UBERON:0002807",
    "UBERON:0002420",
    "UBERON:0001950",
    "UBERON:0002240",
    "UBERON:0001893",
    "UBERON:0016530",
    "UBERON:0016525",
    "UBERON:0002266",
    "UBERON:0002803",
    "UBERON:0002037",
    "UBERON:0016538",
    "UBERON:0005383",
    "UBERON:0002436",
    "UBERON:0002661",
    "UBERON:0000955",
    "UBERON:0006514",
    "UBERON:0002021",
    "UBERON:0001885",
    "UBERON:0001870",
    "UBERON:0007628",
  ],
  "UBERON:0004535": [
    "UBERON:0002078",
    "UBERON:0036288",
    "UBERON:0002098",
    "UBERON:0002082",
    "UBERON:0001624",
    "UBERON:0003517",
    "UBERON:0013756",
    "UBERON:0002084",
    "UBERON:0001637",
    "UBERON:0005616",
    "UBERON:0002080",
    "UBERON:0000948",
    "UBERON:0002081",
    "UBERON:0000947",
    "UBERON:0002049",
    "UBERON:0002079",
    "UBERON:0002094",
    "UBERON:0001621",
  ],
  "UBERON:0001009": [
    "UBERON:0002078",
    "UBERON:0036288",
    "UBERON:0002098",
    "UBERON:0002082",
    "UBERON:0001624",
    "UBERON:0003517",
    "UBERON:0013756",
    "UBERON:0002084",
    "UBERON:0001637",
    "UBERON:0005616",
    "UBERON:0002080",
    "UBERON:0000948",
    "UBERON:0002081",
    "UBERON:0000947",
    "UBERON:0002049",
    "UBERON:0002079",
    "UBERON:0002094",
    "UBERON:0001621",
  ],
  "UBERON:0001007": [
    "UBERON:0010033",
    "UBERON:0002107",
    "UBERON:0002106",
    "UBERON:0001159",
    "UBERON:8410000",
    "UBERON:0002116",
    "UBERON:0002110",
    "UBERON:0001154",
    "UBERON:0001043",
    "UBERON:0000400",
    "UBERON:0002114",
    "UBERON:0000160",
    "UBERON:0001736",
    "UBERON:0001153",
    "UBERON:0001976",
    "UBERON:0002108",
    "UBERON:0000945",
    "UBERON:0001723",
    "UBERON:0002115",
    "UBERON:0010032",
    "UBERON:0004648",
    "UBERON:0000059",
    "UBERON:0001831",
    "UBERON:0001156",
    "UBERON:0001155",
    "UBERON:0001052",
    "UBERON:0001117",
    "UBERON:0001157",
  ],
  "UBERON:0000922": [
    "UBERON:0012168",
    "UBERON:0004026",
    "UBERON:0007106",
    "UBERON:0004024",
  ],
  "UBERON:0000949": [
    "UBERON:0018303",
    "UBERON:0002370",
    "UBERON:0002046",
    "UBERON:0002107",
    "UBERON:0000016",
    "UBERON:0001117",
    "UBERON:0000006",
    "UBERON:0002369",
  ],
  "UBERON:0002330": [
    "UBERON:0001911",
    "UBERON:0002107",
    "UBERON:0000017",
    "UBERON:0001817",
    "UBERON:0001736",
    "UBERON:0001117",
    "UBERON:0001831",
  ],
  "UBERON:0002390": [
    "UBERON:0002371",
    "UBERON:0000178",
    "UBERON:0013756",
    "UBERON:0002106",
    "UBERON:0002370",
    "UBERON:0012168",
  ],
  "UBERON:0002405": [
    "UBERON:0001542",
    "UBERON:0002371",
    "UBERON:0003968",
    "UBERON:0002370",
    "UBERON:0002429",
    "UBERON:0002509",
    "UBERON:0002106",
    "UBERON:0039167",
    "UBERON:0000029",
    "UBERON:0007644",
  ],
  "UBERON:0000383": [
    "UBERON:0003661",
    "UBERON:0001388",
    "UBERON:0001134",
    "UBERON:0001630",
    "UBERON:0002385",
    "UBERON:0008612",
    "UBERON:0002382",
    "UBERON:0002378",
    "UBERON:0004648",
    "UBERON:0001103",
  ],
  "UBERON:0001016": [
    "UBERON:0002811",
    "UBERON:0034751",
    "UBERON:0009835",
    "UBERON:0002802",
    "UBERON:0001872",
    "UBERON:0000956",
    "UBERON:0001871",
    "UBERON:0004725",
    "UBERON:0004024",
    "UBERON:0002808",
    "UBERON:0008933",
    "UBERON:0002809",
    "UBERON:0004026",
    "UBERON:0000966",
    "UBERON:0000966 (organoid)",
    "UBERON:0002421",
    "UBERON:0003902",
    "UBERON:0001882",
    "UBERON:0008934",
    "UBERON:0002728",
    "UBERON:0002810",
    "UBERON:0001890",
    "UBERON:0002771",
    "UBERON:0013682",
    "UBERON:0002894",
    "UBERON:0001384",
    "UBERON:0001786",
    "UBERON:0016540",
    "UBERON:0023852",
    "UBERON:0000451",
    "UBERON:0004167",
    "UBERON:0003876",
    "UBERON:0009834",
    "UBERON:0002807",
    "UBERON:0002420",
    "UBERON:0001950",
    "UBERON:0002240",
    "UBERON:0001893",
    "UBERON:0016530",
    "UBERON:0016525",
    "UBERON:0002266",
    "UBERON:0002822",
    "UBERON:0002803",
    "UBERON:0002037",
    "UBERON:0016538",
    "UBERON:0005383",
    "UBERON:0002436",
    "UBERON:0002661",
    "UBERON:0000955",
    "UBERON:0006514",
    "UBERON:0002021",
    "UBERON:0001885",
    "UBERON:0001870",
    "UBERON:0007628",
  ],
  "UBERON:0001008": [
    "UBERON:0003517",
    "UBERON:0001228",
    "UBERON:0001224",
    "UBERON:0000362",
    "UBERON:0009958",
    "UBERON:0001225",
    "UBERON:0002113",
    "UBERON:0001293",
    "UBERON:0000074",
    "UBERON:0000057",
    "UBERON:0001294",
    "UBERON:0000056",
  ],
  "UBERON:0000990": [
    "UBERON:0016632",
    "UBERON:0001296",
    "UBERON:0003889",
    "UBERON:0000473",
    "UBERON:8410010",
    "UBERON:0000992",
    "UBERON:0002367",
    "UBERON:0000995",
    "UBERON:0000995 (organoid)",
    "UBERON:8410026",
    "UBERON:8410025",
    "UBERON:0001987",
    "UBERON:0012648",
    "UBERON:0001295",
    "UBERON:0000002",
    "UBERON:0003428",
  ],
  "UBERON:0001004": [
    "UBERON:0004802",
    "UBERON:0002185",
    "UBERON:0004903",
    "UBERON:0003126",
    "UBERON:0001005",
    "UBERON:0000115",
    "UBERON:0001707",
    "UBERON:0008953",
    "UBERON:0002048",
    "UBERON:0002048 (organoid)",
    "UBERON:0008954",
    "UBERON:0008946",
    "UBERON:0000977",
    "UBERON:0001103",
  ],
  "UBERON:0001032": [
    "UBERON:0010033",
    "UBERON:0000970",
    "UBERON:0002822",
    "UBERON:0000966",
    "UBERON:0000966 (organoid)",
    "UBERON:0003902",
    "UBERON:0001707",
    "UBERON:0000964",
    "UBERON:0000004",
    "UBERON:0001723",
    "UBERON:0001811",
    "UBERON:0001817",
    "UBERON:0010032",
    "UBERON:0013682",
    "UBERON:0001773",
    "UBERON:0001786",
    "UBERON:0007625",
  ],
  "UBERON:0001434": [
    "UBERON:0004339",
    "UBERON:0002228",
    "UBERON:0002371",
    "UBERON:0013706",
  ],
  "UBERON:0000029": [
    "UBERON:0001542",
    "UBERON:0003968",
    "UBERON:0002429",
    "UBERON:0002509",
    "UBERON:0007644",
    "UBERON:0039167",
  ],
  "UBERON:0002048": [
    "UBERON:0008953",
    "UBERON:0008946",
    "UBERON:0000115",
    "UBERON:0008954",
    "UBERON:0002048 (organoid)",
  ],
  "UBERON:0001043": ["UBERON:0004648", "UBERON:0001976"],
  "UBERON:0003889": ["UBERON:0012648", "UBERON:0016632", "UBERON:8410010"],
  "UBERON:0018707": ["UBERON:0002110", "UBERON:0009958"],
  "UBERON:0000178": ["UBERON:0013756", "UBERON:0012168"],
  "UBERON:0000955": [
    "UBERON:0002811",
    "UBERON:0034751",
    "UBERON:0009835",
    "UBERON:0002802",
    "UBERON:0001872",
    "UBERON:0000956",
    "UBERON:0001871",
    "UBERON:0004725",
    "UBERON:0002808",
    "UBERON:0008933",
    "UBERON:0002809",
    "UBERON:0002421",
    "UBERON:0001882",
    "UBERON:0008934",
    "UBERON:0002728",
    "UBERON:0002810",
    "UBERON:0001890",
    "UBERON:0002771",
    "UBERON:0002894",
    "UBERON:0001384",
    "UBERON:0016540",
    "UBERON:0023852",
    "UBERON:0000451",
    "UBERON:0004167",
    "UBERON:0003876",
    "UBERON:0009834",
    "UBERON:0002807",
    "UBERON:0002420",
    "UBERON:0001950",
    "UBERON:0001893",
    "UBERON:0016530",
    "UBERON:0016525",
    "UBERON:0002266",
    "UBERON:0002803",
    "UBERON:0002037",
    "UBERON:0016538",
    "UBERON:0005383",
    "UBERON:0002436",
    "UBERON:0002661",
    "UBERON:0006514",
    "UBERON:0002021",
    "UBERON:0001885",
    "UBERON:0001870",
    "UBERON:0007628",
  ],
  "UBERON:0000310": ["UBERON:0035328"],
  "UBERON:0000970": [
    "UBERON:0002822",
    "UBERON:0000966",
    "UBERON:0000966 (organoid)",
    "UBERON:0003902",
    "UBERON:0000964",
    "UBERON:0001811",
    "UBERON:0001817",
    "UBERON:0007625",
    "UBERON:0013682",
    "UBERON:0001773",
    "UBERON:0001786",
  ],
  "UBERON:0000948": [
    "UBERON:0002082",
    "UBERON:0002080",
    "UBERON:0036288",
    "UBERON:0002079",
    "UBERON:0002078",
    "UBERON:0002094",
    "UBERON:0001621",
    "UBERON:0002098",
    "UBERON:0002084",
    "UBERON:0002081",
  ],
  "UBERON:0000160": [
    "UBERON:0001156",
    "UBERON:0001155",
    "UBERON:0001052",
    "UBERON:0002116",
    "UBERON:0001154",
    "UBERON:0002108",
    "UBERON:0000400",
    "UBERON:0002115",
    "UBERON:0002114",
    "UBERON:0001153",
    "UBERON:0001157",
    "UBERON:0000059",
    "UBERON:0001159",
    "UBERON:8410000",
  ],
  "UBERON:0002113": [
    "UBERON:0003517",
    "UBERON:0001228",
    "UBERON:0001224",
    "UBERON:0000362",
    "UBERON:0001225",
    "UBERON:0001293",
    "UBERON:0000074",
    "UBERON:0001294",
  ],
  "UBERON:0002107": ["UBERON:0001117"],
  "UBERON:0000004": ["UBERON:0001707"],
  "UBERON:0001264": ["UBERON:0000017", "UBERON:0000006", "UBERON:0000016"],
  "UBERON:0002097": [
    "UBERON:0001511",
    "UBERON:0000014",
    "UBERON:0001416",
    "UBERON:0001868",
  ],
  "UBERON:0001723": ["UBERON:0010033", "UBERON:0010032"],
  "UBERON:0000995": [
    "UBERON:0001296",
    "UBERON:0000002",
    "UBERON:0001295",
    "UBERON:0000995 (organoid)",
  ],
  "UBERON:0001013": [
    "UBERON:0015143",
    "UBERON:0002190",
    "UBERON:0003428",
    "UBERON:0001348",
  ],
};
/* eslint-enable sonarjs/no-duplicate-string -- disable no dupes, this value is generated */
/* eslint-enable sort-keys -- disable key order, contents of this file are generated */

/**
 * Tissues to be included for display in tissue system ontology tree.
 */
/* eslint-disable sort-keys -- disabling key order for readability. */
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
/* eslint-enable sort-keys -- disabling key order for readability. */

/**
 * Tissues to be included for display in tissue system ontology tree.
 */
/* eslint-disable sort-keys -- disabling key order for readability. */
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
/* eslint-enable sort-keys -- disabling key order for readability. */

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
    // TODO(cc) remove with #2569/#3138.
    analyticsEvent: EVENTS.FILTER_SELECT_CELL_TYPE,
    categoryFilterId: CATEGORY_FILTER_ID.CELL_TYPE_DEPRECATED,
    filterOnKey: "cell_type",
    label: "Cell Type",
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
    analyticsEvent: EVENTS.FILTER_SELECT_ETHNICITY,
    categoryFilterId: CATEGORY_FILTER_ID.ETHNICITY,
    filterOnKey: "ethnicity",
    label: "Ethnicity",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
    tooltip:
      "Ethnicity only applies to Homo sapiens which is not selected in the Organism filter.",
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
    analyticsEvent: EVENTS.FILTER_SELECT_AUTHORS,
    categoryFilterId: CATEGORY_FILTER_ID.PUBLICATION_AUTHORS,
    filterOnKey: "publicationAuthors",
    label: "Author",
    labelKind: "VALUE",
    matchKind: "INCLUDES_SOME",
    multiselect: true,
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
