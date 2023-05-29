import { TableInstance } from "react-table";
import {
  Categories,
  TableCountSummary,
  Testable,
} from "src/components/common/Filter/common/entities";
import { Grid as StyledGrid } from "./style";
import React from "react";
import HeaderCell from "src/components/common/Grid/components/HeaderCell";
import CountAndTotal from "src/components/common/Grid/components/HeaderCell/components/CountAndTotal";

interface Props<T extends Categories> {
  className?: string;
  tableCountSummary?: TableCountSummary;
  tableInstance: TableInstance<T>;
}

export default function Grid<T extends Categories & Testable>({
  className,
  tableCountSummary,
  tableInstance,
}: Props<T>): JSX.Element {
  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } =
    tableInstance;
  return (
    <StyledGrid {...getTableProps()} className={className}>
      <thead>
        {headerGroups.map((headerGroup) => {
          const { key, ...restHeaderGroupProps } =
            headerGroup.getHeaderGroupProps();
          return (
            <tr key={key} {...restHeaderGroupProps}>
              {headerGroup.headers.map((column) => {
                const { key, ...restColumnHeaderProps } =
                  column.getHeaderProps();
                const { onClick } = column.getSortByToggleProps();
                return (
                  <HeaderCell
                    key={key}
                    alignment={column.alignment}
                    isSortable={column.canSort}
                    isSorted={column.isSorted}
                    isSortedDesc={column.isSortedDesc}
                    label={column.render("Header")}
                    onSort={onClick}
                    tag={
                      column.showCountAndTotal ? (
                        <CountAndTotal tableCountSummary={tableCountSummary} />
                      ) : undefined
                    }
                    {...restColumnHeaderProps}
                  />
                );
              })}
            </tr>
          );
        })}
      </thead>
      <tbody {...getTableBodyProps}>
        {rows.map((row) => {
          prepareRow(row);
          const { key, ...restRowProps } = row.getRowProps();
          return (
            <tr key={key} data-testid={row.original.testId} {...restRowProps}>
              {row.cells.map((cell) => {
                const { key, ...restCellProps } = cell.getCellProps();
                return (
                  <td key={key} {...restCellProps}>
                    {cell.render("Cell")}
                  </td>
                );
              })}
            </tr>
          );
        })}
      </tbody>
    </StyledGrid>
  );
}
