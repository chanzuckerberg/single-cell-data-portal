import styled from "@emotion/styled";
import { fontBodyXs, getFontWeights, fontBodyS, getColors } from "czifui";

export const GeneInfoWrapper = styled.div``;

export const GeneSummary = styled.div`
  ${fontBodyXs}

  padding: 16px 0;
  font-weight: 500;
  color: black;
`;

export const GeneSynonymsWrapper = styled.div`
  display: flex;
  padding-bottom: 16px;
  align-items: center;
`;

export const GeneSynonymsLabel = styled.div`
  ${fontBodyXs}
  ${(props) => {
    const colors = getColors(props);
    const fontWeights = getFontWeights(props);

    return `
      color: ${colors?.gray[500]};
      font-weight: ${fontWeights?.regular};
    `;
  }}
`;

export const GeneSynonyms = styled.div`
  padding: 4px;
  color: black;

  ${fontBodyXs}
  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.regular};
    `;
  }}
`;

export const GeneUrl = styled.div`
  ${fontBodyS}
  font-weight: 500;
  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.primary[400]};
      `;
  }}
`;

export const GeneName = styled.div`
  ${fontBodyS}
  font-weight: 500;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const InfoButtonWrapper = styled.div`
  padding-right: 8px;
  display: flex;
`;

export const WarningBanner = styled.div`
  padding: 8px;
  display: flex;
  align-items: center;
  margin-top: 8px;

  ${fontBodyXs}
  ${(props) => {
    const colors = getColors(props);
    return `
      background-color: ${colors?.warning[100]};
      svg {
        fill: ${colors?.warning[400]}
      }
    `;
  }}
`;
