import { Icon, Tag } from "@czi-sds/components";
import { useState, useEffect, useContext } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import {
  CopyGenesButton,
  QueryGroupSubTitle,
  QueryGroupTitle,
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

interface PathwayRow {
  geneSet: string;
  geneSymbols: string[];
  pValue: number;
  fdrQValue: number;
}

interface Props {
  setIsLoading: (isLoading: boolean) => void;
}
export default function DeResults({ setIsLoading }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { data: rawDifferentialExpressionResults, isLoading } =
    useDifferentialExpression();
  const {
    differentialExpressionResults1: rawDifferentialExpressionResults1,
    differentialExpressionResults2: rawDifferentialExpressionResults2,
    pathwayEnrichmentAnalysisResults1: rawPathwayEnrichmentAnalysisResults1,
    pathwayEnrichmentAnalysisResults2: rawPathwayEnrichmentAnalysisResults2,
  } = rawDifferentialExpressionResults;

  const [showEnrichedPathways, setShowEnrichedPathways] = useState(false);

  const [differentialExpressionResults1, setDifferentialExpressionResults1] =
    useState<DifferentialExpressionRow[]>([]);

  const [differentialExpressionResults2, setDifferentialExpressionResults2] =
    useState<DifferentialExpressionRow[]>([]);

  const [pathwayEnrichmentResults1, setPathwayEnrichmentResults1] = useState<
    PathwayRow[]
  >([]);

  const [pathwayEnrichmentResults2, setPathwayEnrichmentResults2] = useState<
    PathwayRow[]
  >([]);

  const toggleShowEnrichedPathways = () => {
    setShowEnrichedPathways(!showEnrichedPathways);
  };

  const {
    organismId,
    submittedQueryGroups,
    queryGroups,
    submittedQueryGroupsWithNames: queryGroupsWithNames,
  } = useContext(StateContext);

  useEffect(() => {
    if (!organismId || isLoading) return;

    // map ids to name
    const formattedDeResults1 = rawDifferentialExpressionResults1.map(
      (diffExpResult) => {
        return {
          name: diffExpResult.gene_symbol,
          pValue: diffExpResult.p_value,
          effectSize: diffExpResult.effect_size,
        };
      }
    );

    const formattedDeResults2 = rawDifferentialExpressionResults2.map(
      (diffExpResult) => {
        return {
          name: diffExpResult.gene_symbol,
          pValue: diffExpResult.p_value,
          effectSize: diffExpResult.effect_size,
        };
      }
    );
    setDifferentialExpressionResults1(formattedDeResults1);
    setDifferentialExpressionResults2(formattedDeResults2);

    const formattedPathwayResults1 = rawPathwayEnrichmentAnalysisResults1.map(
      (pathwayResult) => {
        return {
          geneSet: pathwayResult.gene_set,
          pValue: pathwayResult.p_value,
          fdrQValue: pathwayResult.fdr_q_value,
          geneSymbols: pathwayResult.gene_symbols.split(";"),
        };
      }
    );

    const formattedPathwayResults2 = rawPathwayEnrichmentAnalysisResults2.map(
      (pathwayResult) => {
        return {
          geneSet: pathwayResult.gene_set,
          pValue: pathwayResult.p_value,
          fdrQValue: pathwayResult.fdr_q_value,
          geneSymbols: pathwayResult.gene_symbols.split(";"),
        };
      }
    );
    setPathwayEnrichmentResults1(formattedPathwayResults1);
    setPathwayEnrichmentResults2(formattedPathwayResults2);
  }, [
    rawDifferentialExpressionResults1,
    rawDifferentialExpressionResults2,
    rawPathwayEnrichmentAnalysisResults1,
    rawPathwayEnrichmentAnalysisResults2,
    isLoading,
  ]);

  const handleCopyGenes1 = () => {
    const genes = differentialExpressionResults1.map((result) => result.name);
    navigator.clipboard.writeText(genes.join(", "));
  };
  const handleCopyGenes2 = () => {
    const genes = differentialExpressionResults2.map((result) => result.name);
    navigator.clipboard.writeText(genes.join(", "));
  };
  const copyGenesHandlers = [handleCopyGenes1, handleCopyGenes2];

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
  }, [queryGroups, submittedQueryGroups]);

  const showEmpty = !submittedQueryGroups;

  useEffect(() => {
    setIsLoading(isLoading);
  }, [isLoading]);

  if (showEmpty || isLoading) {
    return <div />;
  }

  return (
    <div>
      <div onClick={toggleShowEnrichedPathways}>Click me</div>
      {showEnrichedPathways
        ? [differentialExpressionResults1, differentialExpressionResults2].map(
            (results, index) => {
              return (
                <div>
                  <QueryGroupTitle>Query Group {index + 1}</QueryGroupTitle>
                  <QueryGroupSubTitle>
                    {namesToShow[index].join(", ")}
                  </QueryGroupSubTitle>
                  <DifferentialExpressionResultsTable
                    results={results}
                    copyGeneHandler={copyGenesHandlers[index]}
                  />
                </div>
              );
            }
          )
        : [pathwayEnrichmentResults1, pathwayEnrichmentResults2].map(
            (results, index) => {
              return (
                <div>
                  <QueryGroupTitle>Query Group {index + 1}</QueryGroupTitle>
                  <QueryGroupSubTitle>
                    {namesToShow[index].join(", ")}
                  </QueryGroupSubTitle>
                  <PathwayEnrichmentResultsTable results={results} />
                </div>
              );
            }
          )}
    </div>
  );
}

interface DifferentialExpressionResultsTableProps {
  results: DifferentialExpressionRow[];
  copyGeneHandler: () => void;
}
const DifferentialExpressionResultsTable = ({
  results,
  copyGeneHandler,
}: DifferentialExpressionResultsTableProps) => {
  return (
    <TableWrapper>
      {results.length > 0 ? (
        <StyledHTMLTable condensed bordered={false}>
          <thead>
            <tr>
              <td>
                <CopyGenesButton
                  onClick={copyGeneHandler}
                  sdsType="primary"
                  sdsStyle="minimal"
                  isAllCaps={false}
                  startIcon={
                    <Icon sdsIcon="copy" sdsSize="s" sdsType="button" />
                  }
                >
                  Copy Genes
                </CopyGenesButton>
              </td>
            </tr>
            <tr>
              <td>Gene </td>
              <td>P-value</td>
              <td>Effect size</td>
            </tr>
          </thead>
          <tbody>
            {results.map((result) => {
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

interface PathwayEnrichmentResultsTableProps {
  results: PathwayRow[];
}
const PathwayEnrichmentResultsTable = ({
  results,
}: PathwayEnrichmentResultsTableProps) => {
  return (
    <TableWrapper>
      {results.length > 0 ? (
        <StyledHTMLTable condensed bordered={false}>
          <thead>
            <tr>
              <td>Gene set </td>
              <td>P-value</td>
              <td>FDR Q-value</td>
              {/* <td>DE genes</td> */}
            </tr>
          </thead>
          <tbody>
            {results.map((result) => {
              const { geneSet, pValue, fdrQValue } = result;
              return (
                <tr key={geneSet}>
                  <td>{geneSet}</td>
                  <td>{pValue.toPrecision(4)}</td>
                  <td>{fdrQValue.toPrecision(4)}</td>
                  {/* <td>
                    {geneSymbols.map((geneSymbol) => (
                      <Tag label={geneSymbol} />
                    ))}
                  </td> */}
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
