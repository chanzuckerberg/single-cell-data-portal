import { TableInstance } from "react-table";
import { Categories } from "src/components/common/Filter/common/entities";
import { Grid as StyledGrid } from "./style";

interface Props<T extends Categories> {
  className?: string;
  tableInstance: TableInstance<T>;
}

export default function Grid<T extends Categories>({
  className,
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
                return (
                  <th key={key} {...restColumnHeaderProps}>
                    {column.render("Header")}
                  </th>
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
            <tr key={key} {...restRowProps}>
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
