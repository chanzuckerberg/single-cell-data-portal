import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS, fontHeaderXxl, getColors } from "czifui";
import { Tag } from "czifui";

export const TOP_PADDING_PX = 32;

export const CellCardsView = styled.div`
  display: flex;
  flex-direction: row;
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  padding: ${TOP_PADDING_PX}px 0px 32px 40px;
  width: 1080px;
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
