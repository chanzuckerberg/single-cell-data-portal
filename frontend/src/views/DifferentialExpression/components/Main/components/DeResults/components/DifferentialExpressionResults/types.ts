import { QueryGroups } from "src/views/DifferentialExpression/common/store/reducer";
import { DifferentialExpressionRow } from "../../types";
import { Dispatch, SetStateAction } from "react";

export interface Props {
  queryGroups: QueryGroups;
  queryGroupsWithNames: QueryGroups;
  organismId: string;
  sortedAndFilteredResults: DifferentialExpressionRow[];
  nCellsOverlap: number;
  setSearchQuery: Dispatch<SetStateAction<string>>;
  setLfcFilter: Dispatch<SetStateAction<string>>;
  setEffectSizeFilter: Dispatch<SetStateAction<string>>;
  sortDirection: "asc" | "desc";
  setSortDirection: Dispatch<SetStateAction<"asc" | "desc">>;
  errorMessage: string | null;
}
