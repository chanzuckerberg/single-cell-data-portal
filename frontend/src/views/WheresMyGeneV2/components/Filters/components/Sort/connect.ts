import { track } from "src/common/analytics";
import { Props } from "./types";
import { EVENTS } from "src/common/analytics/events";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { useContext, useMemo } from "react";
import {
  CELL_TYPE_OPTIONS,
  DEFAULT_INPUT_DROPDOWN_PROPS,
  GENE_OPTIONS,
} from "./constants";

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
    value: { id?: SORT_BY; name: string } | null
  ): void {
    if (!dispatch || !value || cellTypeSelectedOption.name === value.name)
      return;

    track(EVENTS.WMG_OPTION_SELECT_CELL_TYPES, {
      sort_cell_types_view_option: value.name,
    });

    dispatch(selectSortBy({ cellTypes: value.id as SORT_BY }));
  }

  function genesOnChange(value: { id?: SORT_BY; name: string } | null): void {
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
