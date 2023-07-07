import { CommonThemeProps, fontBodyS, getColors } from "@czi-sds/components";
import styled from "@emotion/styled";

const DIVIDER_WIDTH = 2;

const gray200 = (props: CommonThemeProps) => getColors(props)?.gray[200];

export const TableTitleOuterWrapper = styled.div`
  display: flex;
  justify-content: space-between;
  align-items: center;
  width: 100%;
`;

export const PublicationLinkWrapper = styled.div`
  display: flex;
  column-gap: 2px;
`;

export const StyledHeadCellContent = styled.div`
  display: flex;
  flex-direction: row;
`;

export const TableSelectorRow = styled.div`
  display: flex;
  flex-direction: row;
  gap: 16px;
  border-bottom: ${DIVIDER_WIDTH}px solid ${gray200};
`;

interface TableSelectorButtonProps extends CommonThemeProps {
  isActive: boolean;
}
export const TableSelectorButton = styled.button<TableSelectorButtonProps>`
  ${fontBodyS}
  position: relative;
  background: none;
  border: none;
  cursor: pointer;
  padding: 0;
  margin: 0;
  font-weight: 600;
  color: ${(props) =>
    `${props.isActive ? "#000000" : getColors(props)?.gray[500]}`};

  &::after {
    content: "";
    position: absolute;
    bottom: -2px;
    left: 0;
    width: 100%;
    height: 2px;
    background-color: ${(props) =>
      `${props.isActive ? getColors(props)?.primary[400] : gray200}`};
    transition: background-color 0.3s ease;
  }
`;
