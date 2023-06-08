import React, { ReactNode } from "react";
import {
  Header,
  HeaderCell as Cell,
  SortIcon,
} from "src/components/common/Grid/components/HeaderCell/style";
import { TableSortByToggleProps } from "react-table";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";
import { Icon } from "@czi-sds/components";

interface Props {
  alignment?: ALIGNMENT;
  isSortable: boolean;
  isSorted: boolean;
  isSortedDesc: boolean | undefined;
  label: ReactNode;
  onSort: TableSortByToggleProps["onClick"];
  tag?: ReactNode;
}

export default function HeaderCell({
  alignment = ALIGNMENT.LEFT,
  isSortable,
  isSorted,
  isSortedDesc,
  label,
  onSort,
  tag,
}: Props): JSX.Element {
  return (
    <Cell alignment={alignment} isSortable={isSortable} isSorted={isSorted}>
      <Header onClick={onSort}>
        {label}
        {/* Sort icon is displayable when the column is sortable - and the corresponding icon displayed is chevron down. */}
        {isSortable ? (
          <SortIcon>
            <Icon
              sdsIcon={isSortedDesc === false ? "chevronUp" : "chevronDown"} // Show chevron down when isSortedDesc is true or undefined.
              sdsSize="s"
              sdsType="interactive"
            />
          </SortIcon>
        ) : null}
      </Header>
      {tag}
    </Cell>
  );
}
