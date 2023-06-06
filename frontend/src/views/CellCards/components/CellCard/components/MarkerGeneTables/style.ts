import { CommonThemeProps, fontBodyS, getColors } from "@czi-sds/components";
import styled from "@emotion/styled";

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

export const TableSelectorRow = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: 16px;
`;

interface TableSelectorButtonProps extends CommonThemeProps {
  isActive: boolean;
}
export const TableSelectorButton = styled.button<TableSelectorButtonProps>`
  ${fontBodyS}
  background: none;
  border: none;
  cursor: pointer;
  padding: 0;
  margin: 0;
  font-weight: 600;

  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${props.isActive ? "#000000" : colors?.gray[500]};
    `;
  }}
`;
