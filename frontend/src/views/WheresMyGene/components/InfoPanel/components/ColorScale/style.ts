import styled from "@emotion/styled";
import { fontBodyS, getFontWeights } from "czifui";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  padding-top: 15px;
`;

export const Label = styled("span")`
  ${fontBodyS};

  ${(props) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.semibold};
    `;
  }}
`;

export const LabelWrapper = styled("div")`
  display: flex;
`;

export const FlexDiv = styled.div`
  padding-top: 4px;
  padding-left: 4px;
`;
