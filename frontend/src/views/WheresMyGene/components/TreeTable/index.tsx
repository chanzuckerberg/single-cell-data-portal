import { memo, useEffect, useRef, useState } from "react";
import { Cell as ICell, Column, useBlockLayout, useTable } from "react-table";
import { VariableSizeGrid } from "react-window";
import { CellTypeAndGenes } from "../../common/types";
import AsterChart from "./components/AsterChart";

interface Props {
  columns: Column<CellTypeAndGenes>[];
  data: CellTypeAndGenes[];
}

export default memo(function TreeTable({ columns, data }: Props): JSX.Element {
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
      <div {...getTableBodyProps({ style: { flex: 1 } })} ref={tableContentRef}>
        <VariableSizeGrid
          height={(tableContentRect?.height || 200) * 0.97}
          width={tableContentRect?.width || 200}
          rowCount={rows.length + 1}
          rowHeight={(index) => (index ? 35 : 100)}
          columnCount={columns.length}
          columnWidth={(index) => (index ? 20 : 200)}
        >
          {({ rowIndex, columnIndex, style }) => {
            if (rowIndex === 0) {
              const headerGroup = headerGroups[0];
              const header = headerGroup.headers[columnIndex];

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
            }

            const row = rows[rowIndex - 1];

            if (row) {
              prepareRow(row);
            }

            const cell = row.cells[columnIndex];
            const { column } = cell;

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
                <CellContent
                  key={`row-${row.id}-column-${column.id}`}
                  cell={cell}
                />
              </div>
            );
          }}
        </VariableSizeGrid>
      </div>
    </div>
  );
});

interface CellContentProps {
  cell: ICell<CellTypeAndGenes>;
}

const CellContent = memo(function Cell_({ cell }: CellContentProps) {
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
        <AsterChart
          colorValue={roundNumber(meanExpression)}
          degreeValue={roundNumber(percentCells)}
        />
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
});

function roundNumber(num: number): number {
  return Math.round(num * 100) / 100;
}
