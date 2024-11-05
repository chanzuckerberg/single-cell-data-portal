export const CELL_GUIDE_CARD_CL_DESCRIPTION = "cell-guide-card-cl-description";

export const CELL_GUIDE_CARD_GPT_DESCRIPTION =
  "cell-guide-card-gpt-description";

export const CELL_GUIDE_CARD_GPT_TOOLTIP_LINK =
  "cell-guide-card-gpt-tooltip-link";

export const DESCRIPTION_BREAKPOINT_HEIGHT_PX = 400;

export const CELL_GUIDE_CARD_SYNONYMS = "cell-guide-card-synonyms";

export const CELL_GUIDE_CARD_VALIDATED_DESCRIPTION =
  "cell-guide-card-validated-description";

export const CELL_GUIDE_CARD_DEFAULT_DESCRIPTION_CL =
  "Description not available";

export const getDefaultGptDescription = (cellTypeName: string) =>
  `Description for ${cellTypeName} is not available at the moment, please check back at a later time, or click on the link below for more information.`;
