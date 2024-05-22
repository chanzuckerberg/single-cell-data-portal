import { EVENTS } from "src/common/analytics/events";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";

export const QUERY_GROUP_KEYS = [
  "tissues",
  "cellTypes",
  "publicationCitations",
  "diseases",
  "ethnicities",
  "sexes",
];

export const QUERY_GROUP_LABELS = [
  "Tissue",
  "Cell Type",
  "Publications",
  "Disease",
  "Ethnicity",
  "Sex",
];

export const QUERY_GROUP_KEYS_TO_FILTER_EVENT_MAP: Partial<
  Record<keyof QueryGroup, EVENTS>
> = {
  tissues: EVENTS.DE_SELECT_TISSUE,
  cellTypes: EVENTS.DE_SELECT_CELL_TYPE,
  publicationCitations: EVENTS.DE_SELECT_PUBLICATION,
  diseases: EVENTS.DE_SELECT_DISEASE,
  ethnicities: EVENTS.DE_SELECT_ETHNICITY,
  sexes: EVENTS.DE_SELECT_SEX,
};
