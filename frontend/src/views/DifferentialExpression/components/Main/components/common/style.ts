import styled from "@emotion/styled";
import {
  Dropdown,
  fontBodyS,
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
  max-width: 216px;
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
` as typeof Dropdown;
