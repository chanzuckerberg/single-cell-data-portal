import { Cell as ICell, Column, useBlockLayout, useTable } from "react-table";
import { FixedSizeList } from "react-window";
import { CellTypeAndGenes } from "../../common/types";
import AsterChart from "./components/AsterChart";
import { Loader } from "./style";

interface Props {
  columns: Column<CellTypeAndGenes>[];
  data: CellTypeAndGenes[];
}

export default function TreeTable({ columns, data }: Props): JSX.Element {
  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } =
    useTable(
      {
        columns,
        data,
      },
      useBlockLayout
    );

  return (
    <div {...getTableProps()}>
      <div>
        {headerGroups.map((headerGroup) => {
          return (
            // eslint-disable-next-line react/jsx-key -- getHeaderGroupProps already has `key` prop
            <div {...headerGroup.getHeaderGroupProps()}>
              {headerGroup.headers.map((header) => {
                return (
                  // eslint-disable-next-line react/jsx-key -- getHeaderProps already has `key` prop
                  <div
                    {...header.getHeaderProps([
                      {
                        style: {
                          textOrientation: "sideways",
                          transform: "translateX(-5px) rotateZ(180deg)",
                          writingMode: "vertical-lr",
                        },
                      },
                    ])}
                  >
                    {header.render("Header")}
                  </div>
                );
              })}
            </div>
          );
        })}
      </div>
      <div {...getTableBodyProps()}>
        <FixedSizeList
          height={600}
          width={1350}
          itemCount={rows.length}
          itemSize={35}
        >
          {({ index, style }) => {
            const row = rows[index];
            prepareRow(row);

            return (
              // eslint-disable-next-line react/jsx-key -- getRowProps already has `key` prop
              <div
                {...row.getRowProps([{ style: { ...style, margin: "8px 0" } }])}
              >
                {row.cells.map((cell) => {
                  const { row, column } = cell;

                  return (
                    <Cell
                      key={`row-${row.id}-column-${column.id}`}
                      cell={cell}
                    />
                  );
                })}
              </div>
            );
          }}
        </FixedSizeList>
      </div>
    </div>
  );
}

interface CellProps {
  cell: ICell<CellTypeAndGenes>;
}

function Cell({ cell }: CellProps) {
  const { value, column } = cell;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any -- Figure out how to extend column type
  const { isLoading } = column as any;

  let Component = null;

  if (typeof value === "string") {
    Component = cell.render("Cell");
  } else if (value === undefined) {
    Component = isLoading ? <Loader /> : null;
  } else {
    const { me: meanExpression, pc: percentCells } = value;

    if (meanExpression !== undefined) {
      Component = (
        <AsterChart colorValue={meanExpression} degreeValue={percentCells} />
      );
    }
  }

  return (
    // eslint-disable-next-line react/jsx-key -- getCellProps already has key
    <div {...cell.getCellProps()}>{Component}</div>
  );
}
