import { useEffect, useRef, useState } from "react";
import { Cell as ICell, Column, useBlockLayout, useTable } from "react-table";
import { FixedSizeList } from "react-window";
import { CellTypeAndGenes } from "../../common/types";
import AsterChart from "./components/AsterChart";

interface Props {
  columns: Column<CellTypeAndGenes>[];
  data: CellTypeAndGenes[];
}

export default function TreeTable({ columns, data }: Props): JSX.Element {
  const tableContentRef = useRef<HTMLDivElement>(null);
  const [tableContentRect, setTableContentRect] =
    useState<DOMRectReadOnly | null>(null);

  useEffect(() => {
    const resizeObserver = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setTableContentRect(entry.contentRect);
      }
    });

    const current = tableContentRef.current;

    if (current) {
      resizeObserver.observe(current);
    }

    return () => {
      if (!current) return;
      resizeObserver.unobserve(current);
    };
  }, []);

  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } =
    useTable(
      {
        columns,
        data,
      },
      useBlockLayout
    );

  return (
    <div
      {...getTableProps({
        style: {
          display: "flex",
          flexDirection: "column",
          height: "100%",
          width: "100%",
        },
      })}
    >
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
                          maxHeight: "100px",
                          overflow: "auto",
                          textOrientation: "sideways",
                          transform: "translateX(-5px) rotateZ(180deg)",
                          whiteSpace: "nowrap",
                          writingMode: "vertical-lr",
                        },
                        // (thuang): HTML attributes not typed in react-table
                        ...{ title: header.Header },
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
      <div {...getTableBodyProps({ style: { flex: 1 } })} ref={tableContentRef}>
        <FixedSizeList
          height={(tableContentRect?.height || 200) * 0.97}
          width={tableContentRect?.width || 200}
          itemCount={rows.length}
          itemSize={35}
        >
          {({ index, style }) => {
            const row = rows[index];
            prepareRow(row);

            return (
              // eslint-disable-next-line react/jsx-key -- getRowProps already has `key` prop
              <div
                {...row.getRowProps([
                  {
                    style: {
                      ...style,
                      margin: "8px 0",
                    },
                  },
                ])}
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
  const { value } = cell;

  let Component = null;

  if (typeof value === "string") {
    Component = cell.render("Cell");
  } else if (value === undefined) {
    Component = null;
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
    <div
      {...cell.getCellProps({
        style: { overflow: "auto", whiteSpace: "nowrap" },
      })}
    >
      {Component}
    </div>
  );
}
