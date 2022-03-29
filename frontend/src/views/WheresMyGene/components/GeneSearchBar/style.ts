import styled from "@emotion/styled";
import { fontBodyXxxs, getColors } from "czifui";

export const Container = styled.div`
  width: 80vw;
  margin-bottom: 20px;
`;

export const ActionWrapper = styled.div`
  display: flex;
  gap: 16px;
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
