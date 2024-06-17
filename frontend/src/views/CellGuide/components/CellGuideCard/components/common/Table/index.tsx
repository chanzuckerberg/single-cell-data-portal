import React, { ReactNode } from "react";
import {
  StyledTable,
  StyledHead,
  StyledHeadCell,
  StyledRow,
  StyledCell,
  TableWrapper,
} from "./style";
import { Tooltip } from "@czi-sds/components";

interface TableProps<T> {
  columns: Array<Extract<keyof T, string>>;
  rows: T[];
  columnIdToName?: Record<
    Extract<keyof T, string>,
    string | React.ReactElement
  >;
  testId?: string;
  hoverable?: boolean;
  columnIdToNumCharactersTruncateThreshold?: Partial<
    Record<Extract<keyof T, string>, number>
  >;
}

// This is a generic table component that can be used to render any type of data.
function Table<T extends object>({
  columns,
  rows,
  columnIdToName,
  testId,
  hoverable = true,
  columnIdToNumCharactersTruncateThreshold,
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
            <StyledRow
              key={rowIndex}
              hoverable={hoverable}
              highlight={rowIndex % 2 === 1}
            >
              {columns.map((column, cellIndex) => {
                const numCharactersTruncateThreshold =
                  columnIdToNumCharactersTruncateThreshold?.[column] ??
                  Number.POSITIVE_INFINITY;

                const text = row[column] as string;
                const truncatedText =
                  text.length > (numCharactersTruncateThreshold ?? 0)
                    ? `${text.slice(0, numCharactersTruncateThreshold)}...`
                    : text;
                const showTooltip =
                  text.length > (numCharactersTruncateThreshold ?? 0);

                return (
                  <Tooltip
                    key={cellIndex}
                    sdsStyle="light"
                    placement="left"
                    leaveDelay={0}
                    title={text}
                    disableHoverListener={!showTooltip}
                  >
                    <StyledCell>{truncatedText as ReactNode}</StyledCell>
                  </Tooltip>
                );
              })}
            </StyledRow>
          ))}
        </tbody>
      </StyledTable>
    </TableWrapper>
  );
}

export default Table;
