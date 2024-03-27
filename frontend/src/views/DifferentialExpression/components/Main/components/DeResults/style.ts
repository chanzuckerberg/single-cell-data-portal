import { HTMLTable } from "@blueprintjs/core";
import styled from "@emotion/styled";
import {
  Button,
  fontBodyL,
  fontBodyS,
  fontBodyXs,
  fontHeaderL,
  fontHeaderXl,
  getColors,
} from "@czi-sds/components";
import { fontBodyXxs } from "@czi-sds/components";
import { gray500 } from "src/common/theme";

export const CopyGenesButton = styled(Button)`
  ${fontBodyXxs}
  font-weight: 500;
  margin-left: -5px;
`;

export const TableWrapper = styled.div`
  width: 386px;
`;

export const StyledHTMLTable = styled(HTMLTable)`
  & thead td {
    ${fontBodyS}
    font-weight: 600;
  }
  & td:nth-of-type(3) {
    text-align: end;
  }

  & tr:not(thead tr) {
    border-bottom: 1px solid #cccccc;
    height: 10px;
  }
`;

export const QueryGroupTitle = styled.div`
  ${fontBodyXs}

  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const QueryGroupSubTitle = styled.div`
  ${fontBodyL}
  font-weight: 400;
`;

export const ButtonsWrapper = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: 8px;
`;

export const RestartButton = styled(Button)`
  margin-bottom: 30px;
`;

export const BackButton = styled(Button)`
  margin-bottom: 30px;
`;

export const StyledLink = styled.a`
  color: black;

  &:hover {
    color: black;
    text-decoration: none;
  }
`;

export const ButtonWrapper = styled.div`
  display: flex;
  justify-content: flex-end;
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
`;
export const ResultsHeader = styled.div`
  ${fontHeaderXl}
  margin-bottom: 16px;
`;
