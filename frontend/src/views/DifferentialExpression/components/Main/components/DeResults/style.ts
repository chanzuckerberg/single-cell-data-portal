import { HTMLTable } from "@blueprintjs/core";
import styled from "@emotion/styled";
import {
  Button,
  CommonThemeProps,
  fontBodyL,
  fontBodyS,
  fontBodyXs,
  getColors,
} from "czifui";
import { fontBodyXxs } from "czifui";

export const CopyGenesButton = styled(Button)`
  ${fontBodyXxs}
  font-weight: 500;
  margin-left: -5px;
`;

export const TableWrapper = styled.div`
  width: 375px;
  max-height: 501px;
  min-height: 300px;
  overflow-y: scroll;
  margin-bottom: 40px;
`;

export const StyledHTMLTable = styled(HTMLTable)`
  & thead td {
    ${fontBodyS}
    // change this to czif font weight?
    font-weight: 600;
  }
  & td:nth-of-type(3) {
    text-align: end;
  }

  & tr:not(thead tr) {
    border-bottom: 1px solid #cccccc;
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

interface Props extends CommonThemeProps {
  isError?: boolean;
}

export const NoDeGenesContainer = styled("div")`
  margin-top: 16px;
  background: #f8f8f8;

  width: 100%;
  height: 120px;

  display: flex;
  flex-direction: column;

  justify-content: center;
  text-align: center;
`;

export const NoDeGenesHeader = styled("span")`
  ${fontBodyS}
  font-weight: 500;
`;

export const NoDeGenesDescription = styled("span")<Props>`
  ${fontBodyXxs}
  ${(props) => {
    const colors = getColors(props);
    return `
    color: ${props.isError ? colors?.error[400] : colors?.gray[500]};
    font-weight: ${props.isError ? 600 : "unset"};
    `;
  }}
`;
