import React, { ReactElement, useState } from "react";
import {
  StyledTable,
  StyledHead,
  StyledHeadCell,
  StyledRow,
  StyledCell,
  TableWrapper,
} from "./style";

interface TableProps<T> {
  columns: Array<Extract<keyof T, string>>;
  rows: T[];
}

type SortOrder = "asc" | "desc" | "none";

// This is a generic table component that can be used to render any type of data.
function Table<T extends object>({ columns, rows }: TableProps<T>) {
  const [sortColumn, setSortColumn] = useState<Extract<keyof T, string> | null>(
    null
  );
  const [sortOrder, setSortOrder] = useState<SortOrder>("none");

  const handleSort = (column: Extract<keyof T, string>) => {
    if (sortColumn === column) {
      setSortOrder((prevSortOrder) => {
        if (prevSortOrder === "asc") return "desc";
        if (prevSortOrder === "desc") return "none";
        return "asc";
      });
    } else {
      setSortColumn(column);
      setSortOrder("asc");
    }
  };

  let sortedRows = [...rows];
  if (sortOrder !== "none" && sortColumn !== null) {
    sortedRows.sort((a, b) => {
      if (sortColumn === null) return 0;
      if (a[sortColumn] < b[sortColumn]) return sortOrder === "asc" ? -1 : 1;
      if (a[sortColumn] > b[sortColumn]) return sortOrder === "asc" ? 1 : -1;
      return 0;
    });
  }

  return (
    <TableWrapper>
      <StyledTable>
        <StyledHead>
          <tr>
            {columns.map((column, index) => (
              <StyledHeadCell key={index} onClick={() => handleSort(column)}>
                {column.charAt(0).toUpperCase() + column.slice(1)}
              </StyledHeadCell>
            ))}
          </tr>
        </StyledHead>
        <tbody>
          {sortedRows.map((row, rowIndex) => (
            <StyledRow key={rowIndex} highlight={rowIndex % 2 === 1}>
              {columns.map((column, cellIndex) => (
                <StyledCell key={cellIndex}>{row[column]}</StyledCell>
              ))}
            </StyledRow>
          ))}
        </tbody>
      </StyledTable>
    </TableWrapper>
  );
}

export default Table;
