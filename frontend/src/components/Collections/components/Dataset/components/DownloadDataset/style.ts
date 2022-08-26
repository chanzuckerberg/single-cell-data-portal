import styled from "@emotion/styled";
import { buttonStyle, disabledButtonStyle } from "../../common/style";

export const StyledButton = styled.button`
  ${buttonStyle}

  &:disabled {
    ${disabledButtonStyle}
  }
`;
