import styled from "@emotion/styled";
import { fontBodyS, fontBodyXxs, getColors, getFontWeights } from "czifui";

export const Container = styled.div`
  width: 80vw;
  margin-bottom: 20px;
`;

export const ActionWrapper = styled.div`
  display: flex;
  gap: 16px;
`;

export const Label = styled.label`
  ${fontBodyXxs}

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;

export const OrganismLabel = styled.label`
  ${fontBodyS}

  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.semibold};
    `;
  }}
`;

export const LoadingIndicatorWrapper = styled.div`
  display: flex;
  align-items: center;
`;
