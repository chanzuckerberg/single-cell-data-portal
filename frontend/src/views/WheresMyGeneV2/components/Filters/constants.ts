import { EVENTS } from "src/common/analytics/events";
import { Filters as IFilters } from "src/views/WheresMyGeneV2/common/types";

export const ANALYTICS_MAPPING: {
  [key in keyof IFilters]: { eventName: EVENTS; label: string };
} = {
  datasets: {
    eventName: EVENTS.FILTER_SELECT_DATASET,
    label: "dataset_name",
  },
  diseases: {
    eventName: EVENTS.FILTER_SELECT_DISEASE,
    label: "disease",
  },
  ethnicities: {
    eventName: EVENTS.FILTER_SELECT_SELF_REPORTED_ETHNICITY,
    label: "ethnicity",
  },
  publications: {
    eventName: EVENTS.FILTER_SELECT_PUBLICATION,
    label: "publication",
  },
  sexes: {
    eventName: EVENTS.FILTER_SELECT_SEX,
    label: "gender",
  },
  tissues: {
    eventName: EVENTS.FILTER_SELECT_TISSUE,
    label: "tissue",
  },
};
