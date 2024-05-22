import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";

export interface FilterOption {
  name: string;
  id: string;
  unavailable?: boolean;
}

export interface Props {
  queryGroup: QueryGroup;
  isQueryGroup1: boolean;
}
