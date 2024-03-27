import { useState, useEffect, useContext } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { StyledHTMLTable, TableWrapper } from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { clearSubmittedQueryGroups } from "src/views/DifferentialExpression/common/store/actions";
import { Pagination } from "@mui/material";

interface DifferentialExpressionRow {
  name: string;
  pValue: number;
  effectSize: number;
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
          pValue: diffExpResult.p_value,
          effectSize: diffExpResult.effect_size,
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

  const showEmpty = !submittedQueryGroups;

  useEffect(() => {
    setIsLoading(isLoading);
  }, [isLoading, setIsLoading]);

  if (showEmpty || isLoading) {
    return <div />;
  }

  return (
    <div>
      <DifferentialExpressionResultsTable
        results={differentialExpressionResults}
      />
    </div>
  );
}

const ROWS_PER_PAGE = 25;

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
    <TableWrapper>
      <StyledHTMLTable bordered={false}>
        <thead>
          <tr>
            <td>Gene </td>
            <td>P-value</td>
            <td>Effect size</td>
          </tr>
        </thead>
        <tbody>
          {results
            .slice((page - 1) * ROWS_PER_PAGE, page * ROWS_PER_PAGE)
            .map((result) => {
              const { name: symbol, pValue, effectSize } = result;
              return (
                <tr key={symbol}>
                  <td>{symbol}</td>
                  <td>{pValue.toPrecision(4)}</td>
                  <td>{effectSize.toPrecision(4)}</td>
                </tr>
              );
            })}
        </tbody>
      </StyledHTMLTable>
      <Pagination count={pageCount} page={page} onChange={handlePageChange} />
    </TableWrapper>
  );
};
