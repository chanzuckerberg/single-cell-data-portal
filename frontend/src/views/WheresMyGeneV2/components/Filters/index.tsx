import {
  ComplexFilterInputDropdown,
  DefaultMenuSelectOption,
} from "@czi-sds/components";
import { memo } from "react";
import Organism from "./components/Organism";
import Compare from "./components/Compare";
import Sort from "./components/Sort";
import {
  FilterContainer,
  StyledComplexFilter,
  StyledComplexFilterInputDropdown,
  ViewOptionsLabel,
  Wrapper,
} from "./style";
import ColorScale from "./components/ColorScale";
import { ViewOptionsWrapper } from "./components/Sort/style";
import { Props } from "./types";
import { useConnect } from "./connect";

export default memo(function Filters({
  isLoading,
  availableFilters,
  setAvailableFilters,
  setIsScaled,
}: Props): JSX.Element {
  const {
    DropdownMenuProps,
    isHeatmapShown,
    InputDropdownProps,
    selected,
    handle,
    terms,
  } = useConnect({
    availableFilters,
    setAvailableFilters,
  });

  return (
    <Wrapper>
      <FilterContainer>
        <StyledComplexFilter
          multiple
          data-testid="disease-filter"
          search
          label="Disease"
          options={terms.disease as unknown as DefaultMenuSelectOption[]}
          onChange={handle.diseasesChange}
          value={selected.diseases as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />
        <StyledComplexFilter
          multiple
          data-testid="self-reported-ethnicity-filter"
          search
          label="Self-Reported Ethnicity"
          options={
            terms.self_reported_ethnicity as unknown as DefaultMenuSelectOption[]
          }
          onChange={handle.ethnicitiesChange}
          value={selected.ethnicities as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />

        <StyledComplexFilter
          multiple
          data-testid="publication-filter"
          search
          label="Publication"
          options={terms.publication as unknown as DefaultMenuSelectOption[]}
          onChange={handle.publicationsChange}
          value={selected.publications as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />

        <StyledComplexFilter
          multiple
          data-testid="sex-filter"
          search
          label="Sex"
          options={terms.sex as unknown as DefaultMenuSelectOption[]}
          onChange={handle.sexesChange}
          value={selected.sexes as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />
        <StyledComplexFilter
          multiple
          data-testid="tissue-filter"
          search
          label="Tissue"
          options={terms.tissue as unknown as DefaultMenuSelectOption[]}
          onChange={handle.tissuesChange}
          value={selected.tissues as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />
      </FilterContainer>

      <Organism isLoading={isLoading} />

      <Compare areFiltersDisabled={false} />

      <FilterContainer>
        <ViewOptionsLabel>View Options</ViewOptionsLabel>
        <ViewOptionsWrapper>
          <Sort areFiltersDisabled={!isHeatmapShown} />
          <ColorScale setIsScaled={setIsScaled} />
        </ViewOptionsWrapper>
      </FilterContainer>
    </Wrapper>
  );
});
