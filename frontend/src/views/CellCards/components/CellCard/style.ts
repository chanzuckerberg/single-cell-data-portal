import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyS,
  fontHeaderM,
  fontHeaderXxl,
  getColors,
} from "czifui";
import { Tag } from "czifui";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  padding: 32px 0px 32px 40px;
  width: 1079px;
`;

export const SearchBarWrapper = styled.div`
  width: 320px;
  margin-bottom: 20px;
`;

export const CellCardHeader = styled.div`
  display: flex;
  column-gap: 8px;
  flex-direction: row;
  align-items: center;
`;

export const CellCardName = styled.div`
  ${fontHeaderXxl}
  font-weight: 700;
`;

export const StyledTag = styled(Tag)`
  height: 24px;
  margin: 0;
`;

export const Divider = styled.div<CommonThemeProps>`
  width: 100%;
  height: 0.5px;
  margin-top: 8px;
  margin-bottom: 8px;

  ${(props) => {
    const colors = getColors(props);
    return `background-color: ${colors?.gray[300]}`;
  }}
`;

export const CellCardDescription = styled.div`
  ${fontBodyS}
  font-wieght: 400;
  white-space: pre-wrap;
`;

export const TableTitle = styled.div`
  ${fontHeaderM}
  font-weight: 600;
  margin-bottom: 8px;
`;

export const TableTitleWrapper = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin-top: 48px;
`;

export const WmgLink = styled.a`
  ${fontBodyS}
  font-weight: 500;
  ${(props) => {
    const colors = getColors(props);
    return `color: ${colors?.primary[400]}`;
  }}
`;
