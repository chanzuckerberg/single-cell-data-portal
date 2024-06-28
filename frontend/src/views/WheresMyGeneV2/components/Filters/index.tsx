import { DefaultMenuSelectOption } from "@czi-sds/components";
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
        <StyledComplexFilter<DefaultMenuSelectOption, true, false, false>
          multiple
          data-testid="disease-filter"
          search
          label="Disease"
          options={terms.disease}
          onChange={handle.diseasesChange}
          value={selected.diseases}
          InputDropdownComponent={StyledComplexFilterInputDropdown}
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />
        <StyledComplexFilter<DefaultMenuSelectOption, true, false, false>
          multiple
          data-testid="self-reported-ethnicity-filter"
          search
          label="Self-Reported Ethnicity"
          options={terms.self_reported_ethnicity}
          onChange={handle.ethnicitiesChange}
          value={selected.ethnicities}
          InputDropdownComponent={StyledComplexFilterInputDropdown}
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />

        <StyledComplexFilter<DefaultMenuSelectOption, true, false, false>
          multiple
          data-testid="publication-filter"
          search
          label="Publication"
          options={terms.publication}
          onChange={handle.publicationsChange}
          value={selected.publications}
          InputDropdownComponent={StyledComplexFilterInputDropdown}
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />

        <StyledComplexFilter<DefaultMenuSelectOption, true, false, false>
          multiple
          data-testid="sex-filter"
          search
          label="Sex"
          options={terms.sex}
          onChange={handle.sexesChange}
          value={selected.sexes}
          InputDropdownComponent={StyledComplexFilterInputDropdown}
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
        />
        <StyledComplexFilter<DefaultMenuSelectOption, true, false, false>
          multiple
          data-testid="tissue-filter"
          search
          label="Tissue"
          options={terms.tissue}
          onChange={handle.tissuesChange}
          value={selected.tissues}
          InputDropdownComponent={StyledComplexFilterInputDropdown}
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
