import styled from "@emotion/styled";
import {
  CommonThemeProps,
  ComplexFilter,
  ComplexFilterInputDropdown,
  ComplexFilterPopper,
  fontBodyS,
  getColors,
  getFontWeights,
  TagFilter,
  Tag,
  fontBodyXs,
} from "czifui";

export const StyledPopper = styled(ComplexFilterPopper)`
  max-width: 750px;
`;

export const StyledComplexFilterInputDropdown = styled(
  ComplexFilterInputDropdown
)`
  padding-left: 0;
  /* (thuang): Sex filter is short and doesn't need the default 64px width */
  min-width: 0;

  .MuiButton-label {
    margin-left: 0;
  }

  &.Mui-disabled {
    border: 0;
  }

  && :nth-last-child(1) {
    display: none !important;
  }

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    return `
      .styled-label {
        // (thuang): Override the default color black to be gray until SDS fixes it
        color: ${colors?.gray[500]} !important;
      }
    `;
  }}
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: row;
  width: 444px;
  flex-wrap: wrap;
  column-gap: 16px;

  // (atarashansky): remove the blue chips from the filters
  div > div > .MuiChip-root {
    display: none;
  }
`;

export const FilterLabel = styled.div`
  ${fontBodyS}

  margin-bottom: 8px;

  ${(props: CommonThemeProps) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.semibold};
    `;
  }}
`;

export const StyledComplexFilter = styled(ComplexFilter)`
  margin-bottom: 16px;
  min-width: 75px;
` as typeof ComplexFilter;

export const StyledTagFilter = styled(TagFilter)`
  ${fontBodyXs}
  height: 28px;
  border-radius: 4px;
`;

export const StyledTag = styled(Tag)`
  ${fontBodyXs}
  height: 28px;
  border-radius: 4px;
`;

export const TagWrapper = styled.div`
  width: 100%;
`;
