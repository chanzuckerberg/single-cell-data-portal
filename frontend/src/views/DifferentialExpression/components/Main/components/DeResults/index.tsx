import { useState, useEffect, useContext } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import {
  InstructionsBody,
  InstructionsHeader,
  InstructionsWrapper,
  ResultsHeader,
  ResultsWrapper,
  TableWrapper,
} from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { clearSubmittedQueryGroups } from "src/views/DifferentialExpression/common/store/actions";
import { Pagination } from "@mui/material";
import Table from "src/views/CellGuide/components/CellGuideCard/components/common/Table";

interface DifferentialExpressionRow {
  name: string;
  pValue: string;
  effectSize: string;
}

interface Props {
  setIsLoading: (isLoading: boolean) => void;
}
export default function DeResults({ setIsLoading }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { data, isLoading } = useDifferentialExpression();
  const { differentialExpressionResults: rawDifferentialExpressionResults } =
    data;

  const [differentialExpressionResults, setDifferentialExpressionResults] =
    useState<DifferentialExpressionRow[]>([]);

  const {
    organismId,
    submittedQueryGroups,
    queryGroups,
    submittedQueryGroupsWithNames: queryGroupsWithNames,
  } = useContext(StateContext);

  useEffect(() => {
    if (!organismId || isLoading) return;

    // map ids to name
    const formattedDeResults = rawDifferentialExpressionResults.map(
      (diffExpResult) => {
        return {
          name: diffExpResult.gene_symbol,
          pValue: diffExpResult.p_value.toExponential(3),
          effectSize: diffExpResult.effect_size.toFixed(3),
        };
      }
    );

    setDifferentialExpressionResults(formattedDeResults);
  }, [rawDifferentialExpressionResults, isLoading, organismId]);

  const namesToShow: string[][] = [];
  const { queryGroup1, queryGroup2 } = queryGroupsWithNames ?? {};
  for (const [index, queryGroupWithNames] of [
    queryGroup1,
    queryGroup2,
  ].entries()) {
    namesToShow.push([]);
    for (const key in queryGroupWithNames) {
      for (const value of queryGroupWithNames[key as keyof QueryGroup]) {
        namesToShow[index].push(value);
      }
    }
  }

  useEffect(() => {
    if (!submittedQueryGroups || !dispatch) return;
    if (JSON.stringify(queryGroups) !== JSON.stringify(submittedQueryGroups))
      dispatch(clearSubmittedQueryGroups());
  }, [dispatch, queryGroups, submittedQueryGroups]);

  const showEmpty = !submittedQueryGroups || isLoading;

  useEffect(() => {
    setIsLoading(isLoading);
  }, [isLoading, setIsLoading]);

  console.log(differentialExpressionResults);
  return (
    <div>
      {!showEmpty ? (
        <DifferentialExpressionResultsTable
          results={differentialExpressionResults}
        />
      ) : (
        <InstructionsWrapper>
          <InstructionsHeader>Instructions</InstructionsHeader>
          <InstructionsBody>
            <ol>
              <li>
                Select a cell group of interest within the Cell Group 1 box by
                using the dropdown selectors.
                <br />
                <br />
                To copy the same selection over to Cell Group 2, click the copy
                button to the right of each dropdown in Cell Group 1.
              </li>
              <li>
                Within Cell Group 2, select a group that the cell group of
                interest will be compared to.
                <br />
                <br />
                To easily select the inverse of Cell Group 1, click the
                overlapping circle icon to the left of each dropdown in Cell
                Group 2.
              </li>
            </ol>
          </InstructionsBody>
        </InstructionsWrapper>
      )}
    </div>
  );
}

const ROWS_PER_PAGE = 15;

const columnIdToName: Record<keyof DifferentialExpressionRow, string> = {
  name: "Gene",
  pValue: "P-value",
  effectSize: "Effect Size",
};

interface DifferentialExpressionResultsTableProps {
  results: DifferentialExpressionRow[];
}
const DifferentialExpressionResultsTable = ({
  results,
}: DifferentialExpressionResultsTableProps) => {
  const [page, setPage] = useState(1);
  const pageCount = Math.ceil(results.length / ROWS_PER_PAGE);
  const handlePageChange = (
    _event: React.ChangeEvent<unknown>,
    page: number
  ) => {
    setPage(page);
  };

  return (
    <ResultsWrapper>
      <ResultsHeader>Results</ResultsHeader>
      <TableWrapper>
        <Table<DifferentialExpressionRow>
          columns={["name", "pValue", "effectSize"]}
          rows={results.slice((page - 1) * ROWS_PER_PAGE, page * ROWS_PER_PAGE)}
          columnIdToName={columnIdToName}
        />

        <Pagination count={pageCount} page={page} onChange={handlePageChange} />
      </TableWrapper>
    </ResultsWrapper>
  );
};
