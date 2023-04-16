import styled from "@emotion/styled";
import {
  fontHeaderXl,
  fontBodyS,
  getColors,
  CommonThemeProps,
  Button,
} from "czifui";

export const StepHeader = styled.div`
  ${fontHeaderXl}

  margin-bottom: 14px;
  margin-top: 6px;
`;

export const WordPop = styled.span<CommonThemeProps>`
  ${(props) => {
    const colors = getColors(props);
    return `
        color: ${colors?.primary[400]};
    `;
  }}
`;

export const StepSubHeader = styled.div`
  ${fontBodyS}

  margin-bottom: 32px;

  ${(props) => {
    const colors = getColors(props);
    return `
        color: ${colors?.gray[500]};
    `;
  }}
`;

export const FiltersWrapper = styled.div`
  min-width: 216px;
  margin-bottom: 42px;
`;

export const NextButton = styled(Button)`
  margin-bottom: 50px;
`;
