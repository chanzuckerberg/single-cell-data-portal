import React, { ReactNode } from "react";
import {
  StyledTable,
  StyledHead,
  StyledHeadCell,
  StyledRow,
  StyledCell,
  TableWrapper,
  StyledHeadCellContent,
} from "./style";
import HelpTooltip from "../HelpTooltip";

interface TableProps<T> {
  columns: Array<Extract<keyof T, string>>;
  rows: T[];
  columnIdToName?: Record<Extract<keyof T, string>, string>;
}

// This is a generic table component that can be used to render any type of data.
function Table<T extends object>({
  columns,
  rows,
  columnIdToName,
}: TableProps<T>) {
  return (
    <TableWrapper>
      <StyledTable>
        <StyledHead>
          <tr>
            {columns.map((column, index) => {
              const columnName = columnIdToName
                ? columnIdToName[column]
                : column.charAt(0).toUpperCase() + column.slice(1);

              return (
                <StyledHeadCell key={index}>
                  <StyledHeadCellContent>
                    <div>{columnName}</div>
                    {columnName === "Expression Score" ? (
                      <HelpTooltip dark text="Expression Score" />
                    ) : (
                      columnName === "% of Cells" && (
                        <HelpTooltip dark text="% of Cells" />
                      )
                    )}
                  </StyledHeadCellContent>
                </StyledHeadCell>
              );
            })}
          </tr>
        </StyledHead>
        <tbody>
          {rows.map((row, rowIndex) => (
            <StyledRow key={rowIndex} highlight={rowIndex % 2 === 1}>
              {columns.map((column, cellIndex) => (
                <StyledCell key={cellIndex}>
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
