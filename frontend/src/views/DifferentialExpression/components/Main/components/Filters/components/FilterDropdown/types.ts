import { FilterOption } from "../../types";

export interface Props {
  label: string;
  options: FilterOption[];
  allAvailableOptions: FilterOption[];
  selectedOptionIds: string[];
  handleChange: (options: FilterOption[]) => void;
}
