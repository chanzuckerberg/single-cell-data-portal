import styled from "@emotion/styled";
import { Dropdown, fontBodyXxxs, getColors, getFontWeights } from "czifui";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const Label = styled("span")`
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

export const LabelWrapper = styled("div")`
  display: flex;
`;

export const FlexDiv = styled.div`
  padding-top: 4px;
  padding-left: 4px;
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
` as typeof Dropdown;
