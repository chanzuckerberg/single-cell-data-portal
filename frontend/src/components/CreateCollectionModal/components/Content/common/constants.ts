export enum CONSORTIA {
  ALLEN_INSTITUTE_FOR_BRAIN_SCIENCE = "Allen Institute for Brain Science",
  BRAIN_INITIATIVE = "BRAIN Initiative",
  CZ_BIOHUB = "CZ Biohub",
  CZI_CS = "CZI Cell Science",
  CZI_NDCN = "CZI Neurodegeneration Challenge Network",
  EU_HORIZON_2020 = "European Union’s Horizon 2020",
  GUDMAP = "GenitoUrinary Development Molecular Anatomy Project (GUDMAP)",
  GUT_CELL_ATLAS = "Gut Cell Atlas",
  HUBMAP = "Human BioMolecular Atlas Program (HuBMAP)",
  HCA = "Human Cell Atlas (HCA)",
  HPAP = "Human Pancreas Analysis Program (HPAP)",
  HTAN = "Human Tumor Atlas Network (HTAN)",
  KPMP = "Kidney Precision Medicine Project (KPMP)",
  NEPTUNE = "Nephrotic Syndrome Study Network (NEPTUNE)",
  LUNG_MAP = "LungMAP",
  PCEN = "Pediatric Center of Excellence in Nephrology (PCEN)",
  SEA_AD = "SEA-AD",
  WELCOME_HCA_STRATEGIC_SCIENCE_SUPPORT = "Wellcome HCA Strategic Science Support",
}

export const DEBOUNCE_TIME_MS = 100;

/**
 * Text displayed when BE has identified DOI as invalid.
 */
export const INVALID_DOI_ERROR_MESSAGE =
  "This DOI could not be found. Please correct or remove it.";

export const INVALID_DATASET_STATUS_FOR_DOI_UPDATE =
  "Cannot update DOI while a dataset is processing or awaiting upload.";
