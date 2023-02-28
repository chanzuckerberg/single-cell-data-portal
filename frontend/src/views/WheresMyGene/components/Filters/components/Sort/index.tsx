import { InputDropdownProps as IInputDropdownProps } from "czifui";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { FilterLabel, FilterWrapper, Label, StyledDropdown } from "./style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

const CELL_TYPE_OPTIONS = [
  { id: SORT_BY.CELL_ONTOLOGY, name: "Cell Ontology" },
  { id: SORT_BY.H_CLUSTER, name: "Hierarchical" },
];

const GENE_OPTIONS = [
  { id: SORT_BY.USER_ENTERED, name: "As Entered" },
  { id: SORT_BY.H_CLUSTER, name: "Hierarchical" },
];

interface Props {
  areFiltersDisabled: boolean;
}

export default function Sort({ areFiltersDisabled }: Props): JSX.Element {
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

  return (
    <div>
      <Label>View Options</Label>
      <FilterWrapper>
        <FilterLabel>Sort Cell Types</FilterLabel>
        <StyledDropdown
          data-test-id="cell-type-sort-dropdown"
          onChange={cellTypesOnChange}
          label={cellTypeSelectedOption.name}
          options={CELL_TYPE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={cellTypeSelectedOption}
        />
      </FilterWrapper>
      <FilterWrapper>
        <FilterLabel>Sort Genes</FilterLabel>
        <StyledDropdown
          data-test-id="gene-sort-dropdown"
          onChange={genesOnChange}
          label={geneSelectedOption.name}
          options={GENE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={geneSelectedOption}
        />
      </FilterWrapper>
    </div>
  );

  function cellTypesOnChange(
    value: { id?: SORT_BY; name: string } | null
  ): void {
    if (!dispatch || !value) return;

    track(EVENTS.WMG_OPTION_SELECT_CELL_TYPES, {
      sort_cell_types_view_option: value.name,
    });

    dispatch(selectSortBy({ cellTypes: value.id as SORT_BY }));
  }

  function genesOnChange(value: { id?: SORT_BY; name: string } | null): void {
    if (!dispatch || !value) return;

    track(EVENTS.WMG_OPTION_SELECT_SORT_GENES, {
      sort_genes_view_option: value.name,
    });

    dispatch(selectSortBy({ genes: value.id as SORT_BY }));
  }
}
