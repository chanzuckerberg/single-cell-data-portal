// Core dependencies
import { Dispatch, SetStateAction, useState } from "react";

/**
 * Function invoked to update state for the filter search value.
 */
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;

/**
 * Function invoked to clear state for the filter search value.
 */
type ClearSearchValueFn = () => void;

/**
 * Model of the filter search state value and the corresponding functions to update the value.
 */
export interface FilterSearchState {
  clearSearchValueFn: ClearSearchValueFn;
  searchValue: string;
  setSearchValue: SetSearchValueFn;
}

/**
 * Faceted filter search functionality.
 * @returns the current search value, a function to update current search value, and a function to clear search value.
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
