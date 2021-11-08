import { interpolateViridis } from "d3-scale-chromatic";
import { useBlockLayout, useExpanded, useTable } from "react-table";
import AsterChart from "./components/AsterChart";
import { Square } from "./style";

interface Props {
  columns: any;
  data: any;
}

export default function TreeTable({ columns, data }: Props): JSX.Element {
  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } =
    useTable(
      {
        columns,
        data,
      },
      useExpanded,
      useBlockLayout
    );

  return (
    <div {...getTableProps()}>
      <div>
        {headerGroups.map((headerGroup) => {
          return (
            // eslint-disable-next-line react/jsx-key -- getHeaderGroupProps already has key
            <div {...headerGroup.getHeaderGroupProps()}>
              {headerGroup.headers.map((header) => {
                return (
                  // eslint-disable-next-line react/jsx-key -- getHeaderProps already has key
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
        {rows.map((row) => {
          prepareRow(row);

          return (
            // eslint-disable-next-line react/jsx-key -- getRowProps already has key
            <div {...row.getRowProps([{ style: { margin: "8px 0" } }])}>
              {row.cells.map((cell) => {
                return <Cell key={cell.key} cell={cell} />;
              })}
            </div>
          );
        })}
      </div>
    </div>
  );
}

function Cell({ cell }) {
  const { value } = cell;

  let Component = null;

  if (typeof value === "string") {
    Component = cell.render("Cell");
  } else {
    const { negEntropy, relativeExpression, proportionalExpression } = value;

    if (negEntropy !== undefined) {
      Component = <Square color={interpolateViridis(negEntropy)} />;
    }

    if (relativeExpression !== undefined) {
      Component = (
        <AsterChart
          relativeExpression={relativeExpression}
          proportionalExpression={proportionalExpression}
        />
      );
    }
  }

  return (
    // eslint-disable-next-line react/jsx-key -- getCellProps already has key
    <div {...cell.getCellProps()}>{Component}</div>
  );
}
