import { Icon } from "czifui";
import { useState, useEffect, useContext } from "react";
import {
  OntologyTerm,
  useDifferentialExpression,
  usePrimaryFilterDimensions,
} from "src/common/queries/differentialExpression";
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

export default function DeResults(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { data: rawDifferentialExpressionResults, isLoading } =
    useDifferentialExpression();
  const {
    differentialExpressionResults1: rawDifferentialExpressionResults1,
    differentialExpressionResults2: rawDifferentialExpressionResults2,
  } = rawDifferentialExpressionResults;

  const [differentialExpressionResults1, setDifferentialExpressionResults1] =
    useState<DifferentialExpressionRow[]>([]);

  const [differentialExpressionResults2, setDifferentialExpressionResults2] =
    useState<DifferentialExpressionRow[]>([]);
  const { data, isLoading: isLoadingPrimaryFilters } =
    usePrimaryFilterDimensions();
  const { genes: rawGenes } = data || {};
  const {
    organismId,
    submittedQueryGroups,
    queryGroups,
    submittedQueryGroupsWithNames: queryGroupsWithNames,
  } = useContext(StateContext);

  useEffect(() => {
    if (!rawGenes || isLoadingPrimaryFilters || !organismId || isLoading)
      return;
    const genes = rawGenes[organismId];
    const genesById = genes.reduce((acc, gene) => {
      return acc.set(gene.id, gene);
    }, new Map<OntologyTerm["id"], OntologyTerm>());

    // map ids to name
    const formattedResults1 = rawDifferentialExpressionResults1.map(
      (diffExpResult) => {
        return {
          name: genesById.get(diffExpResult.gene_ontology_term_id)?.name ?? "", // nullish coalescing operator for type safety
          pValue: diffExpResult.p_value,
          effectSize: diffExpResult.effect_size,
        };
      }
    );

    const formattedResults2 = rawDifferentialExpressionResults2.map(
      (diffExpResult) => {
        return {
          name: genesById.get(diffExpResult.gene_ontology_term_id)?.name ?? "", // nullish coalescing operator for type safety
          pValue: diffExpResult.p_value,
          effectSize: diffExpResult.effect_size,
        };
      }
    );
    setDifferentialExpressionResults1(formattedResults1);
    setDifferentialExpressionResults2(formattedResults2);
  }, [
    rawDifferentialExpressionResults,
    isLoading,
    isLoadingPrimaryFilters,
    rawGenes,
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
  return (
    <div>
      {[differentialExpressionResults1, differentialExpressionResults2].map(
        (results, index) => {
          return (
            <div>
              <QueryGroupTitle>Query Group {index + 1}</QueryGroupTitle>
              <QueryGroupSubTitle>
                {namesToShow[index].join(", ")}
              </QueryGroupSubTitle>
              <TableWrapper>
                {results.length > 0 ? (
                  <StyledHTMLTable condensed bordered={false}>
                    <thead>
                      <tr>
                        <td>
                          <CopyGenesButton
                            onClick={copyGenesHandlers[index]}
                            sdsType="primary"
                            sdsStyle="minimal"
                            isAllCaps={false}
                            startIcon={
                              <Icon
                                sdsIcon="copy"
                                sdsSize="s"
                                sdsType="button"
                              />
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
                ) : isLoading ? (
                  <div />
                ) : (
                  <NoDeGenesContainer>
                    <NoDeGenesHeader>
                      No Differentially Expressed Genes
                    </NoDeGenesHeader>
                    <NoDeGenesDescription>
                      No differentially expressed genes for this query group.
                    </NoDeGenesDescription>
                  </NoDeGenesContainer>
                )}
              </TableWrapper>
            </div>
          );
        }
      )}
    </div>
  );
}
