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
import { ROUTES } from "src/common/constants/routes";

interface TableProps<T> {
  columns: Array<Extract<keyof T, string>>;
  rows: T[];
  columnIdToName?: Record<Extract<keyof T, string>, string>;
  testId?: string;
}

export const EXPRESSION_SCORE_TOOLTIP_TEST_ID = "expression-score-tooltip";
export const PERCENT_OF_CELLS_TOOLTIP_TEST_ID = "percent-of-cells-tooltip";

// This is a generic table component that can be used to render any type of data.
function Table<T extends object>({
  columns,
  rows,
  columnIdToName,
  testId,
}: TableProps<T>) {
  const expressionScoreTooltip = (
    <HelpTooltip
      dark
      buttonDataTestId={EXPRESSION_SCORE_TOOLTIP_TEST_ID}
      text={
        <div>
          The expression score is the average{" "}
          <a
            href={ROUTES.WMG_DOCS_DATA_PROCESSING}
            target="_blank"
            rel="noreferrer noopener"
          >
            rankit-normalized gene expression
          </a>{" "}
          among cells in the cell type that have non-zero values.
        </div>
      }
    />
  );

  const percentOfCellsTooltip = (
    <HelpTooltip
      dark
      buttonDataTestId={PERCENT_OF_CELLS_TOOLTIP_TEST_ID}
      text={
        <div>
          Percentage of cells expressing a gene in the cell type. These numbers
          are calculated after cells with{" "}
          <a
            href={ROUTES.WMG_DOCS_DATA_PROCESSING}
            target="_blank"
            rel="noreferrer noopener"
          >
            low coverage and low expression values
          </a>{" "}
          have been filtered out.
        </div>
      }
    />
  );

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
                  <StyledHeadCellContent>
                    <div>{columnName}</div>
                    {columnName === "Expression Score"
                      ? expressionScoreTooltip
                      : columnName === "% of Cells" && percentOfCellsTooltip}
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
