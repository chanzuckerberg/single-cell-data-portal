import { InputDropdownProps as IInputDropdownProps } from "czifui";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { FilterLabel, FilterWrapper, Label, StyledDropdown } from "./style";

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

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
        <FilterLabel>Sort Genes</FilterLabel>
        <StyledDropdown
          data-test-id="gene-sort-dropdown"
          onChange={genesOnChange}
          label={geneSelectedOptionLabel}
          options={GENE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
        />
      </FilterWrapper>
    </div>
  );

  function genesOnChange(value: { id?: SORT_BY; name: string } | null): void {
    if (!dispatch || !value) return;

    dispatch(selectSortBy({ genes: value.id as SORT_BY }));
  }
}
