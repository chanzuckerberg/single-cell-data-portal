import React, { ReactNode } from "react";
import {
  HeaderCell as Header,
  LeftAlignCell,
  RightAlignCell,
} from "src/components/common/Grid/components/HeaderCell/style";
import { TableSortByToggleProps } from "react-table";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";

interface Props {
  alignment?: ALIGNMENT;
  children: ReactNode;
  isSortable: boolean;
  onSort: TableSortByToggleProps["onClick"];
  sortIcon?: ReactNode;
  tag?: ReactNode;
}

export default function HeaderCell({
  alignment = ALIGNMENT.LEFT,
  children,
  isSortable,
  onSort,
  sortIcon,
  tag,
}: Props): JSX.Element {
  const Cell = alignment === ALIGNMENT.LEFT ? LeftAlignCell : RightAlignCell;
  return (
    <Cell>
      <Header isSortable={isSortable} onClick={onSort}>
        <span>{children}</span>
        {tag}
        {sortIcon}
      </Header>
    </Cell>
  );
}
