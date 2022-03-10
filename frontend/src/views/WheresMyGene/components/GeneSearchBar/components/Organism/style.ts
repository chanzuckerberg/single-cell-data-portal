import styled from "@emotion/styled";
import { Dropdown, fontBodyXxxs, getColors } from "czifui";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const StyledDropdown = styled(Dropdown)`
  height: 30px;
  width: 135px;
`;

export const Label = styled.label`
  ${fontBodyXxxs}

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;
