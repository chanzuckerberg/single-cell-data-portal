import styled from "@emotion/styled";
import {
  Dropdown,
  fontBodyS,
  fontBodyXxxs,
  getColors,
  getFontWeights,
  getSpaces,
} from "@czi-sds/components";

export const Label = styled("div")`
  ${fontBodyS}

  ${(props) => {
    const fontWeights = getFontWeights(props);
    const spaces = getSpaces(props);

    return `
      font-weight: ${fontWeights?.semibold};
      margin-bottom: ${spaces?.xxs}px;
    `;
  }}
`;

export const Wrapper = styled("div")`
  display: flex;
  flex-direction: column;
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
` as typeof Dropdown;

export const FilterLabel = styled("label")`
  ${fontBodyXxxs}

  ${(props) => {
    const colors = getColors(props);
    const fontWeights = getFontWeights(props);
    const spaces = getSpaces(props);

    return `
      color: ${colors?.gray["500"]};
      font-weight: ${fontWeights?.medium};
      margin-bottom: ${spaces?.xxs}px;
    `;
  }}
`;
