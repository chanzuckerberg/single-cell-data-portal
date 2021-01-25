import { DARK_GRAY } from "src/components/common/theme";
import styled from "styled-components";
import { buttonStyle } from "../../common/style";

export const StyledButton = styled.button`
  ${buttonStyle}

  &:disabled {
    border: 1px solid ${DARK_GRAY.A};
    color: ${DARK_GRAY.A};

    &:hover {
      cursor: not-allowed;
      background-color: white;
    }
  }
`;
