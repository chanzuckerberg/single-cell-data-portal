import styled from "styled-components";
import { buttonStyle, disabledButtonStyle } from "../../common/style";

export const StyledButton = styled.button`
  ${buttonStyle}

  &:disabled {
    ${disabledButtonStyle}
  }
`;
