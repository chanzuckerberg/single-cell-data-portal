import { InputDropdownProps as IInputDropdownProps } from "@czi-sds/components";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { Wrapper, FilterLabel, StyledDropdown } from "../common/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ViewOptionsWrapper } from "./style";

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
    <ViewOptionsWrapper>
      <Wrapper>
        <FilterLabel>Sort Cell Types</FilterLabel>
        <StyledDropdown
          data-testid="cell-type-sort-dropdown"
          onChange={cellTypesOnChange}
          label={cellTypeSelectedOption.name}
          options={CELL_TYPE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={cellTypeSelectedOption}
        />
      </Wrapper>
      <Wrapper>
        <FilterLabel>Sort Genes</FilterLabel>
        <StyledDropdown
          data-testid="gene-sort-dropdown"
          onChange={genesOnChange}
          label={geneSelectedOption.name}
          options={GENE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={geneSelectedOption}
        />
      </Wrapper>
    </ViewOptionsWrapper>
  );

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
}
