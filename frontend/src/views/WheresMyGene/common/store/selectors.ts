import {
  INITIAL_STATE,
  State,
} from "src/views/WheresMyGene/common/store/reducer";

export function selectHasCustomFiltersOrGenesSelected(state: State): boolean {
  const {
    selectedFilters,
    selectedGenes,
    compare,
    sortBy,
    filteredCellTypes,
    filteredCellTypeIds,
  } = state;

  // check if any arrays in selectedFilters is not empty using .some()
  const hasCustomFilters = Object.values(selectedFilters).some(
    (filter) => filter.length
  );

  return (
    hasCustomFilters ||
    selectedGenes.length > 0 ||
    Boolean(compare) ||
    !isSortByDefault(sortBy) ||
    filteredCellTypes.length > 0 ||
    filteredCellTypeIds.length > 0
  );
}

/**
 * (thuang): Returns true if the sortBy value is the same as the default value.
 */
function isSortByDefault(sortBy: State["sortBy"]): boolean {
  const {
    cellTypes: initialCellTypes,
    genes: initialGenes,
    scaled: initialScaled,
  } = INITIAL_STATE.sortBy;
  const { cellTypes, genes, scaled } = sortBy;

  return (
    cellTypes === initialCellTypes &&
    genes === initialGenes &&
    scaled === initialScaled
  );
}
