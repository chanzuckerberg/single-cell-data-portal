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

function Table<T extends object>({ columns, rows }: TableProps<T>) {
  return (
    <TableWrapper>
      <StyledTable>
        <StyledHead>
          <tr>
            {columns.map((column, index) => (
              <StyledHeadCell key={index}>
                {column.charAt(0).toUpperCase() + column.slice(1)}
              </StyledHeadCell>
            ))}
          </tr>
        </StyledHead>
        <tbody>
          {rows.map((row, rowIndex) => (
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
