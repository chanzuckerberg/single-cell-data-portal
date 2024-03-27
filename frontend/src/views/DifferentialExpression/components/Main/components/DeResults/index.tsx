import { useState, useEffect, useContext } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import {
  StyledHTMLTable,
  TableWrapper,
  NoDeGenesContainer,
  NoDeGenesDescription,
  NoDeGenesHeader,
} from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { clearSubmittedQueryGroups } from "src/views/DifferentialExpression/common/store/actions";

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

interface DifferentialExpressionResultsTableProps {
  results: DifferentialExpressionRow[];
}
const DifferentialExpressionResultsTable = ({
  results,
}: DifferentialExpressionResultsTableProps) => {
  return (
    <TableWrapper>
      {results.length > 0 ? (
        <StyledHTMLTable bordered={false}>
          <thead>
            <tr>
              <td>Gene </td>
              <td>P-value</td>
              <td>Effect size</td>
            </tr>
          </thead>
          <tbody>
            {results.slice(0, 100).map((result) => {
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
      ) : (
        <NoDeGenesContainer>
          <NoDeGenesHeader>No Differentially Expressed Genes</NoDeGenesHeader>
          <NoDeGenesDescription>
            No differentially expressed genes for this query group.
          </NoDeGenesDescription>
        </NoDeGenesContainer>
      )}
    </TableWrapper>
  );
};
