import { Dispatch, SetStateAction } from "react";
import { FilterDimensions } from "src/common/queries/wheresMyGene";

export interface FilterOption {
  name: string;
  label: string;
  id: string;
}

export interface Props {
  isLoading: boolean;
  availableFilters: Partial<FilterDimensions>;
  setAvailableFilters: Dispatch<SetStateAction<Partial<FilterDimensions>>>;
  setIsScaled: Dispatch<SetStateAction<boolean>>;
}
