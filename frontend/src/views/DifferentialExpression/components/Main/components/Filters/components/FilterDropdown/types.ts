import {
  FilterOption,
  QueryGroup,
} from "src/views/DifferentialExpression/common/store/reducer";

export interface Props {
  label: string;
  options: FilterOption[];
  allAvailableOptions: FilterOption[];
  selectedOptionIds: string[];
  handleChange: (options: FilterOption[]) => void;
  isQueryGroup1: boolean;
  queryGroupKey: keyof QueryGroup;
}
