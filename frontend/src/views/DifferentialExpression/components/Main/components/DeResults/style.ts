import styled from "@emotion/styled";
import {
  Icon,
  fontBodyM,
  fontBodyS,
  fontHeaderL,
  fontHeaderXl,
} from "@czi-sds/components";
import { gray500, primary400 } from "src/common/theme";

const TABLE_WIDTH = "386px";

export const ButtonLabel = styled.span`
  ${fontBodyM}
  color: ${primary400};
  font-weight: 500;
`;

export const ResultsHeaderWrapper = styled.div`
  display: flex;
  column-gap: 8px;
  align-items: center;
  flex-direction: row;
  margin-bottom: 16px;
  justify-content: space-between;
  width: ${TABLE_WIDTH};
`;

interface ButtonsWrapperProps {
  disabled?: boolean;
}
export const ButtonsWrapper = styled.div<ButtonsWrapperProps>`
  display: flex;
  column-gap: 7px;
  align-items: center;
  flex-direction: row;
  ${({ disabled }) =>
    disabled
      ? `
    opacity: 0.5;
    cursor: default;
  `
      : `
    cursor: pointer;
  `}
`;

export const InstructionsWrapper = styled.div`
  word-wrap: break-word;
  white-space: pre-wrap;
  max-width: 368px;
  margin-left: 16px;
  margin-top: 16px;
`;
export const InstructionsHeader = styled.div`
  ${fontHeaderL}
`;
export const InstructionsBody = styled.div`
  ${fontBodyS}
  margin-top: 24px;
  color: ${gray500};
`;

export const ResultsWrapper = styled.div`
  margin-left: 40px;
  margin-top: 24px;
  padding-right: 20px;
`;

export const ResultsHeader = styled.div`
  ${fontHeaderXl}
`;

export const StyledIcon = styled(Icon)`
  width: 22px;
  height: 22px;
  flex-shrink: 0;
`;

export const FlexRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  column-gap: 16px;
`;
