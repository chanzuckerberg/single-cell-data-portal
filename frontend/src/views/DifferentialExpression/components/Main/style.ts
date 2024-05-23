import styled from "@emotion/styled";
import { keyframes } from "@emotion/react";
import {
  fontHeaderXl,
  fontBodyS,
  getColors,
  CommonThemeProps,
  Button,
  fontBodyXxs,
} from "@czi-sds/components";

const LEFT_PANEL_WIDTH = "60vw";
const RIGHT_PANEL_WIDTH = "40vw";

export const TwoPanelLayout = styled.div`
  display: flex;
  flex-direction: row;

  .leftPanel {
    border-right: 1px solid #ccc;
    width: ${LEFT_PANEL_WIDTH};
    min-width: fit-content;
  }

  .rightPanel {
    width: ${RIGHT_PANEL_WIDTH};
    min-width: fit-content;
  }
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  margin-left: 120px;
  margin-top: 24px;
  margin-right: 24px;
`;

export const StepHeader = styled.div`
  ${fontHeaderXl}

  margin-bottom: 14px;
  margin-top: 6px;
`;

export const CellGroupTitle = styled.div`
  ${fontBodyS}
  font-weight: 600;
  margin-bottom: 8px;
  margin-top: 25px;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
`;

export const CellCountTitle = styled.div`
  ${fontBodyXxs}
  color: #767676;
  display: flex;
  flex-direction: row;
  align-items: center;
  gap: 4px;
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

export const RunButton = styled(Button)`
  ${fontBodyS}
  margin-bottom: 50px;
  margin-top: 59px;
  width: 256px;
`;

export const RunButtonWrapper = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: flex-end;
`;

export const FlexRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  column-gap: 16px;
`;

export const QuerySelectorWrapper = styled.div`
  width: 698px;
`;

const spin = keyframes`
  from {
    transform: rotate(0deg);
  }
  to {
    transform: rotate(360deg);
  }
`;

export const Spinner = styled.div`
  border: 2px solid rgba(0, 0, 0, 0.2);
  width: 12px;
  height: 12px;
  border-radius: 50%;
  border-left-color: #000;
  animation: ${spin} 0.5s linear infinite;
`;
