import styled from "@emotion/styled";
import {
  Dropdown,
  fontBodyS,
  fontBodyXxxs,
  getColors,
  getFontWeights,
} from "czifui";

export const Label = styled("div")`
  ${fontBodyS}

  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.semibold};
      margin-bottom: 8px;
    `;
  }}
`;

export const FilterLabel = styled("label")`
  ${fontBodyXxxs}

  ${(props) => {
    const colors = getColors(props);
    const fontWeights = getFontWeights(props);

    return `
      color: ${colors?.gray["500"]};
      font-weight: ${fontWeights?.medium};
    `;
  }}
`;

export const FilterWrapper = styled("div")`
  display: flex;
  flex-direction: column;
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
` as typeof Dropdown;
