import styled from "@emotion/styled";
import { fontBodyXxs, getColors, IconButton } from "czifui";

export const StyledButtonDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 0 10px;
`;

export const StyledLabel = styled.div`
  ${fontBodyXxs}

  white-space: nowrap;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;

export const StyledIconButton = styled(IconButton)`
  width: 30px;
  height: 30px;
`;
