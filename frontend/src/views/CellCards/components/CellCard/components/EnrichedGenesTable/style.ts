import styled from "@emotion/styled";

interface TableTitleInnerWrapperProps {
  columnGap?: number;
}

export const TableTitleInnerWrapper = styled.div<TableTitleInnerWrapperProps>`
  display: flex;
  align-items: center;
  ${(props) => {
    const columnGap = props.columnGap ?? 8;
    return `
      column-gap: ${columnGap}px;
    `;
  }}
`;

export const TableTitleOuterWrapper = styled.div`
  display: flex;
  justify-content: space-between;
  align-items: center;
  width: 100%;
`;
