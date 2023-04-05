import styled from "@emotion/styled";

interface Props {
  isSortable: boolean;
}

export const SortIcon = styled.span`
  align-items: center;
  color: inherit;
  display: flex;
  height: 16px;
  justify-content: center;
  width: 16px;

  .MuiSvgIcon-root {
    color: inherit;
  }
`;

export const LeftAlignCell = styled("span")`
  display: flex;
  justify-content: flex-start;
`;

export const RightAlignCell = styled("span")`
  display: flex;
  justify-content: flex-end;

  ${SortIcon} {
    order: -1;
  }
`;

export const HeaderCell = styled("span")<Props>`
  align-items: center;
  cursor: ${({ isSortable }) => (isSortable ? "pointer" : "default")};
  display: flex;
  gap: 4px;

  &:hover {
    color: ${({ isSortable }) => (isSortable ? `#000000` : undefined)};
  }

  span {
    min-width: 0; /* facilitates breaking of word on columns; flex default for min width is "auto" */
  }
`;
