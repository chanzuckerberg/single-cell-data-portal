import styled from "@emotion/styled";
import {
  fontBodyXs,
  getFontWeights,
  fontBodyS,
  getColors,
  Callout,
} from "czifui";

export const GeneInfoWrapper = styled.div``;

export const GeneSummary = styled.div`
  ${fontBodyXs}

  padding: 8px 0;
  font-weight: 400;
  color: black;
`;

export const GeneSynonymsWrapper = styled.div`
  display: flex;
  padding-bottom: 8px;
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
  padding: 4px 8px;
  color: black;

  ${fontBodyXs}
  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.regular};
    `;
  }}
`;

export const GeneUrl = styled.a`
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
  color: black;
  padding-top: 4px;
`;

export const StyledCallout = styled(Callout)`
  ${fontBodyXs}
  align-items: center;
  width: 100%;
  margin-top: 8px;
  margin-bottom: 0;
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
