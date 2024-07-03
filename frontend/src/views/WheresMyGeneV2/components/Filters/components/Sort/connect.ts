import { track } from "src/common/analytics";
import { CellTypeOptionType, Props } from "./types";
import { EVENTS } from "src/common/analytics/events";
import { selectSortBy } from "src/views/WheresMyGeneV2/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGeneV2/common/types";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import { useContext, useMemo } from "react";
import {
  CELL_TYPE_OPTIONS,
  DEFAULT_INPUT_DROPDOWN_PROPS,
  GENE_OPTIONS,
} from "./constants";
import { AutocompleteValue } from "@mui/base";

export const useConnect = ({ areFiltersDisabled }: Props) => {
  const dispatch = useContext(DispatchContext);
  const { sortBy } = useContext(StateContext);

  const InputDropdownProps = useMemo(
    () => ({
      ...DEFAULT_INPUT_DROPDOWN_PROPS,
      disabled: areFiltersDisabled,
    }),
    [areFiltersDisabled]
  );

  const cellTypeSelectedOption = useMemo(() => {
    return (
      CELL_TYPE_OPTIONS.find((option) => option.id === sortBy.cellTypes) ||
      CELL_TYPE_OPTIONS[0]
    );
  }, [sortBy]);

  const geneSelectedOption = useMemo(() => {
    return (
      GENE_OPTIONS.find((option) => option.id === sortBy.genes) ||
      GENE_OPTIONS[0]
    );
  }, [sortBy]);

  function cellTypesOnChange(
    _: React.SyntheticEvent, // event
    value: AutocompleteValue<CellTypeOptionType, false, false, false>
  ): void {
    if (!dispatch || !value || cellTypeSelectedOption.name === value.name)
      return;

    track(EVENTS.WMG_OPTION_SELECT_CELL_TYPES, {
      sort_cell_types_view_option: value.name,
    });

    dispatch(selectSortBy({ cellTypes: value.id as SORT_BY }));
  }

  function genesOnChange(
    _: React.SyntheticEvent, // event
    value: AutocompleteValue<CellTypeOptionType, false, false, false>
  ): void {
    if (!dispatch || !value || geneSelectedOption.name === value.name) return;

    track(EVENTS.WMG_OPTION_SELECT_SORT_GENES, {
      sort_genes_view_option: value.name,
    });

    dispatch(selectSortBy({ genes: value.id as SORT_BY }));
  }

  return {
    InputDropdownProps,
    cellTypeSelectedOption,
    geneSelectedOption,
    cellTypesOnChange,
    genesOnChange,
  };
};
