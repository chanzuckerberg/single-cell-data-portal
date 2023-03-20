import styled from "@emotion/styled";
import {
  fontBodyXs,
  getFontWeights,
  fontBodyS,
  getColors,
  fontHeaderL,
} from "czifui";

export const GeneInfoWrapper = styled.div`
  position: fixed;
  border-radius: 3px 3px 0px 0px;
  padding: 0px 20px 20px 0px;
  background: white;
  box-shadow: 0px 0px 3px 2px rgba(153, 153, 153, 0.2);
  z-index: 2;
`;

export const GeneHeader = styled.p`
  font-weight: 500;
  ${fontBodyS}

  ${(props) => {
    const colors = getColors(props);

    return `
        color: ${colors?.gray[500]};
        `;
  }}
`;

export const GeneSymbol = styled.h1`
  color: black;
  ${fontHeaderL}
  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
        font-weight: ${fontWeights?.semibold};
        `;
  }}
`;

export const Content = styled.p`
  font-weight: 500;
  color: black;
  ${fontBodyXs}
`;

export const SynHeader = styled.span`
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

export const Synonyms = styled.span`
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

export const Link = styled.a`
  font-weight: 500;
  ${fontBodyS}
  ${(props) => {
    const colors = getColors(props);

    return `
        color: ${colors?.primary[400]};
        `;
  }}
`;

export const WarningBanner = styled.div`
  padding: 8px;
  display: flex;
  align-items: center;
  span {
    margin-left: 10px;
  }
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
