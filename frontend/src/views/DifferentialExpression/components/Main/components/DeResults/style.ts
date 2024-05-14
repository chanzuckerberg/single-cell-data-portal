import { HTMLTable } from "@blueprintjs/core";
import styled from "@emotion/styled";
import {
  Button,
  Tag,
  fontBodyL,
  fontBodyS,
  fontBodyXs,
  fontHeaderL,
  fontHeaderXl,
  getColors,
} from "@czi-sds/components";
import { fontBodyXxs } from "@czi-sds/components";
import { gray100, gray500 } from "src/common/theme";
import { TextField } from "@mui/material";

export const CopyGenesButton = styled(Button)`
  ${fontBodyXxs}
  font-weight: 500;
  margin-left: -5px;
`;

const TABLE_WIDTH = "386px";
export const TableWrapper = styled.div`
  width: ${TABLE_WIDTH};
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

export const StyledTextField = styled(TextField)`
  height: 32px;
  max-width: 140px;
  margin-top: 4px;
  & .MuiInputBase-root {
    padding: 0;
    height: 32px;
  }
  & .MuiInputLabel-root {
    margin-top: -8px;
    z-index: 0;
  }
`;

export const TableHeaderWrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const EffectSizeHeaderWrapper = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: 8px;
  align-items: center;
`;

export const CellGroupWrapper = styled.div`
  height: 92px;
  width: ${TABLE_WIDTH};
  margin-bottom: 8px;

  padding: 12px;
  background-color: ${gray100};
`;

export const CellGroupTitleWrapper = styled.span`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: center;
`;

export const CellGroupTitle = styled.span`
  ${fontBodyS}
  font-weight: 600;
`;

export const EffectSizeIndicator = styled.span`
  ${fontBodyXxs}
  color: #959595;
`;

export const FilterTagsWrapper = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: 8px;
  margin-top: 8px;
`;

export const StyledTag = styled(Tag)`
  font-weight: 400;
  padding: 2px 8px;
  background-color: rgba(0, 0, 0, 0.05);
  & .MuiChip-label {
    color: #000000;
    ${fontBodyXs}
  }
`;
