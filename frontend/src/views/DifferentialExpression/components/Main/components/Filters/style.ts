import styled from "@emotion/styled";
import {
  CommonThemeProps,
  ComplexFilter,
  ComplexFilterInputDropdown,
  ComplexFilterPopper,
  getColors,
  TagFilter,
  fontBodyXs,
  fontCapsXxs,
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

export const Wrapper = styled.div<CommonThemeProps>`
  display: flex;
  flex-direction: row;
  width: 100%;
  flex-wrap: wrap;
  column-gap: 16px;
  padding-top: 12px;
  padding-left: 12px;
  padding-right: 12px;
  margin-bottom: 8px;
  background-color: #f8f8f8;

  // (atarashansky): remove the blue chips from the filters
  div > div > .MuiChip-root {
    display: none;
  }

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);
    return `
      border: 1px solid ${colors?.gray[200]};
      border-radius: 6px;
    `;
  }}
`;

export const ClearButtonWrapper = styled.div`
  ${fontCapsXxs}
  font-weight: 600;
  width: 100%;
  justify-content: flex-end;
  display: flex;
  flex-direction: row;
  cursor: pointer;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[400]}
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
  min-width: 75px;

  button {
    padding: 0;
    > span > span {
      font-weight: 500 !important;
      font-size: 14px !important;
      line-height: 20px !important;
      color: #767676 !important;
    }
  }
` as typeof ComplexFilter;

export const StyledTagFilter = styled(TagFilter)`
  ${fontBodyXs}
  height: 28px;
  border-radius: 4px;
`;

export const TagWrapper = styled.div`
  width: 100%;
  margin-bottom: 8px;
  border-bottom: 2px solid #eaeaea;
  min-height: 127px;
`;
