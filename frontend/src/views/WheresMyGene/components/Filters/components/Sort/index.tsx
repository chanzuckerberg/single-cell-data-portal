import { Dropdown, InputDropdownProps as IInputDropdownProps } from "czifui";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { FilterLabel, FilterWrapper, Label } from "./style";

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

  const cellTypeSelectedOptionLabel = useMemo(() => {
    return (
      CELL_TYPE_OPTIONS.find((option) => option.id === sortBy.cellTypes)
        ?.name || CELL_TYPE_OPTIONS[0].name
    );
  }, [sortBy]);

  const geneSelectedOptionLabel = useMemo(() => {
    return (
      GENE_OPTIONS.find((option) => option.id === sortBy.genes)?.name ||
      GENE_OPTIONS[0].name
    );
  }, [sortBy]);

  return (
    <div>
      <Label>View Options</Label>
      <FilterWrapper>
        <FilterLabel>Sort Cell Types</FilterLabel>
        <Dropdown
          data-test-id="cell-type-sort-dropdown"
          onChange={cellTypesOnChange}
          label={cellTypeSelectedOptionLabel}
          options={CELL_TYPE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
        />
      </FilterWrapper>
      <FilterWrapper>
        <FilterLabel>Sort Genes</FilterLabel>
        <Dropdown
          data-test-id="gene-sort-dropdown"
          onChange={genesOnChange}
          label={geneSelectedOptionLabel}
          options={GENE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
        />
      </FilterWrapper>
    </div>
  );

  function cellTypesOnChange(
    value: { id?: SORT_BY; name: string } | null
  ): void {
    if (!dispatch || !value) return;

    dispatch(selectSortBy({ cellTypes: value.id as SORT_BY }));
  }

  function genesOnChange(value: { id?: SORT_BY; name: string } | null): void {
    if (!dispatch || !value) return;

    dispatch(selectSortBy({ genes: value.id as SORT_BY }));
  }
}
