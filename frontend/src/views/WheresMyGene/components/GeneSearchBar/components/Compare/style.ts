import styled from "@emotion/styled";
import {
  Dropdown,
  fontBodyS,
  fontBodyXxs,
  getColors,
  getFontWeights,
} from "czifui";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  margin: -30px 0;
`;

export const CompareLabel = styled.label`
  ${fontBodyS}

  margin-top: 40px;

  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.semibold};
    `;
  }}
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
  height: 30px;
  min-width: 135px;
`;
