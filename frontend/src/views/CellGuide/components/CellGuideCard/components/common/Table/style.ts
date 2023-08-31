import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyS,
  getColors,
  getSpaces,
} from "@czi-sds/components";
import { spacesL, spacesM, spacesS } from "src/common/theme";

export const StyledTable = styled.table`
  width: 100%;
  table-layout: auto;
`;

export const TableWrapper = styled.div`
  overflow-y: auto;
`;

export const StyledHead = styled.thead`
  height: 24px;
  cursor: default;
`;

export const StyledHeadCell = styled.th`
  ${fontBodyS}
  font-weight: 500;
  white-space: nowrap;
  ${(props) => {
    const colors = getColors(props);
    const spacings = getSpaces(props);

    return `
    color: ${colors?.gray[500]};
    padding: ${spacings?.s}px ${spacings?.m}px ${spacings?.s}px ${spacings?.m}px;`;
  }}
`;

interface StyledRowProps {
  highlight: boolean;
}
export const StyledRow = styled.tr<StyledRowProps>`
  background-color: ${(props) => (props.highlight ? "#F8F8F8" : "white")};
`;

export enum PaddingType {
  None = 0,
  Medium = 1,
  Large = 2,
  MediumLarge = 3,
}

interface StyledCellProps extends CommonThemeProps {
  minWidth?: number;
  maxWidth?: number;
  addPadding?: PaddingType;
}

export const StyledCell = styled.td<StyledCellProps>`
  ${fontBodyS}
  font-weight: 400;
  padding-top: ${spacesS}px;
  padding-bottom: ${spacesS}px;
  padding-right: ${(props) =>
    props.addPadding === PaddingType.None
      ? 0
      : props.addPadding === PaddingType.Medium
      ? spacesM(props)
      : props.addPadding === PaddingType.Large
      ? spacesL(props)
      : props.addPadding === PaddingType.MediumLarge
      ? (spacesM(props) ?? 0) + (spacesL(props) ?? 0)
      : spacesM(props)}px;
  min-width: ${(props) => props.minWidth ?? 120}px;
  max-width: ${(props) => props.maxWidth}px;
  word-break: break-word;
  overflow-wrap: break-word;
  vertical-align: top;

  a {
    display: inline-block;
    overflow-wrap: break-word;
  }
`;
