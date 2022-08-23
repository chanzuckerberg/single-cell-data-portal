// Core dependencies
import { useState } from "react";
import { SetSearchValueFn } from "src/components/common/Filter/common/entities";

type ClearSearchValueFn = () => void;

export interface FilterSearchState {
  clearSearchValueFn: ClearSearchValueFn;
  searchValue: string;
  setSearchValue: SetSearchValueFn;
}

/**
 * Faceted filter search functionality.
 * @returns a function to clear search value, the current search value, and the function to update current search value.
 */
export function useFilterSearch(): FilterSearchState {
  const [searchValue, setSearchValue] = useState<string>("");

  /**
   * Function to clear current search value.
   */
  const clearSearchValueFn = () => {
    setSearchValue("");
  };

  return {
    clearSearchValueFn,
    searchValue,
    setSearchValue,
  };
}
