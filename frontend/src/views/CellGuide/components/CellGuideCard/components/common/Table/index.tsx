import React, { ReactNode } from "react";
import {
  StyledTable,
  StyledHead,
  StyledHeadCell,
  StyledRow,
  StyledCell,
  TableWrapper,
  PaddingType,
} from "./style";

interface TableProps<T> {
  columns: Array<Extract<keyof T, string>>;
  rows: T[];
  columnIdToName?: Record<
    Extract<keyof T, string>,
    string | React.ReactElement
  >;
  testId?: string;
  minWidths?: (number | undefined)[];
  maxWidths?: (number | undefined)[];
  rightPaddings?: PaddingType[];
}

// This is a generic table component that can be used to render any type of data.
function Table<T extends object>({
  columns,
  rows,
  columnIdToName,
  testId,
  minWidths,
  maxWidths,
  rightPaddings,
}: TableProps<T>) {
  return (
    <TableWrapper data-testid={testId}>
      <StyledTable>
        <StyledHead>
          <tr>
            {columns.map((column, index) => {
              const columnName = columnIdToName
                ? columnIdToName[column]
                : column.charAt(0).toUpperCase() + column.slice(1);

              return (
                <StyledHeadCell key={index} id={`${columnName}-header-cell`}>
                  {columnName}
                </StyledHeadCell>
              );
            })}
          </tr>
        </StyledHead>
        <tbody>
          {rows.map((row, rowIndex) => (
            <StyledRow key={rowIndex} highlight={rowIndex % 2 === 1}>
              {columns.map((column, cellIndex) => (
                <StyledCell
                  key={cellIndex}
                  minWidth={minWidths?.[cellIndex]}
                  maxWidth={maxWidths?.[cellIndex]}
                  addPadding={rightPaddings?.[cellIndex]}
                >
                  {row[column] as ReactNode}
                </StyledCell>
              ))}
            </StyledRow>
          ))}
        </tbody>
      </StyledTable>
    </TableWrapper>
  );
}

export default Table;
