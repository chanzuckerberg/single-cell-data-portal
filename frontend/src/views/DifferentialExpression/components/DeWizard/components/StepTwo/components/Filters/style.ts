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

interface WrapperProps extends CommonThemeProps {
  active?: boolean;
}

export const Wrapper = styled.div<WrapperProps>`
  display: flex;
  flex-direction: row;
  width: 444px;
  flex-wrap: wrap;
  column-gap: 16px;
  padding-top: 12px;
  padding-left: 12px;
  padding-right: 12px;
  margin-bottom: 8px;

  // (atarashansky): remove the blue chips from the filters
  div > div > .MuiChip-root {
    display: none;
  }

  ${(props: WrapperProps) => {
    const colors = getColors(props);
    return `
      border: 1px ${props.active ? "solid" : "dashed"} ${colors?.gray[300]};
      border-radius: 6px;
    `;
  }}
`;

export const QueryGroupTitle = styled.div`
  ${fontBodyS}
  margin-bottom: 8px;
  width: 100%;
  display: flex;
  flex-direction: row;
  justify-content: space-between;

  ${(props: CommonThemeProps) => {
    const fontWeights = getFontWeights(props);
    return `
      font-weight: ${fontWeights?.semibold}
    `;
  }}
`;

export const EmptyRectangle = styled.div`
  width: 116px;
  height: 28px;
  border-radius: 4px;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);
    return `
      border: 1px dashed ${colors?.gray[300]};
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

export const TagWrapper = styled.div`
  width: 100%;
  margin-bottom: 8px;
`;

export const IconButtonWrapper = styled.div`
  cursor: pointer;
`;
