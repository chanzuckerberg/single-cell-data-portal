import { FC } from "react";
import { TableInstance } from "react-table";
import { CollectionRow } from "src/components/common/Filter/common/entities";

interface Props {
  tableInstance: TableInstance<CollectionRow>;
}

// TODO(cc) rename to CollectionsGrid and rename existing CollectionsGrid to MyCollectionsGrid.
const FilteredCollectionsGrid: FC<Props> = ({ tableInstance }) => {
  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } =
    tableInstance;
  return (
    <table {...getTableProps()} style={{ fontSize: "14px" }}>
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
    </table>
  );
};

export default FilteredCollectionsGrid;
