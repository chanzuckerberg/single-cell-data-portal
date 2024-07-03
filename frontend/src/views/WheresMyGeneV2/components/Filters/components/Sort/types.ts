import { SORT_BY } from "src/views/WheresMyGeneV2/common/types";

export interface Props {
  areFiltersDisabled: boolean;
}

export type CellTypeOptionType = { id?: SORT_BY; name: string };
