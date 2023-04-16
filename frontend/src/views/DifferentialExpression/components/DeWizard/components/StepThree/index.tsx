import { Icon } from "czifui";
import { useState, useEffect, useContext, useCallback } from "react";
import {
  OntologyTerm,
  useDifferentialExpression,
  usePrimaryFilterDimensions,
} from "src/common/queries/differentialExpression";
import { StateContext } from "src/views/DifferentialExpression/common/store";
import {
  CopyGenesButton,
  MarkerStrengthContainer,
  MarkerStrengthLabel,
  StyledHTMLTable,
} from "./style";
interface Props {
  setStep: (step: number) => void;
}
interface DifferentialExpressionRow {
  name: string;
  pValue: number;
  effectSize: number;
}

export default function StepThree({ setStep }: Props): JSX.Element {
  const { data: rawDifferentialExpressionResults, isLoading } =
    useDifferentialExpression();
  const [differentialExpressionResults, setDifferentialExpressionResults] =
    useState<DifferentialExpressionRow[][]>([]);
  const { data, isLoading: isLoadingPrimaryFilters } =
    usePrimaryFilterDimensions();
  const { genes: rawGenes } = data || {};
  const { organismId, queryGroups } = useContext(StateContext);

  useEffect(() => {
    if (
      !rawGenes ||
      !isLoadingPrimaryFilters ||
      !isLoading ||
      !rawDifferentialExpressionResults.length
    )
      return;
    const genes = rawGenes[organismId || ""];
    const genesById = genes.reduce((acc, gene) => {
      return acc.set(gene.name.toLowerCase(), gene);
    }, new Map<OntologyTerm["id"], OntologyTerm>());

    // map ids to name
    const formattedResults = rawDifferentialExpressionResults.map(
      (diffExpResult) => {
        return diffExpResult.map((result) => {
          return {
            name: genesById.get(result.gene_ontology_term_id)?.name ?? "", // nullish coalescing operator for type safety
            pValue: result.p_value,
            effectSize: result.effect_size,
          };
        });
      }
    );
    setDifferentialExpressionResults(formattedResults);
  }, [
    rawDifferentialExpressionResults,
    isLoading,
    isLoadingPrimaryFilters,
    rawGenes,
  ]);

  const handleCopyGenes = useCallback(
    function handleCopyGenes_(index: number): () => void {
      return () => {
        const results = differentialExpressionResults[index];
        const genes = results.map((result) => result.name);
        navigator.clipboard.writeText(genes.join(", "));
      };
    },
    [differentialExpressionResults]
  );

  const handleGoBack = () => {
    setStep(2);
  };
  const handleStartOver = () => {
    setStep(1);
  };

  return (
    <div>
      {isLoading && <div>Loading...</div>}
      {differentialExpressionResults.map((results, index) => {
        return (
          <StyledHTMLTable condensed bordered={false}>
            <thead>
              <tr>
                <td>Gene </td>
                <td>P-value</td>
                <td>Effect size</td>
              </tr>
              <tr>
                <td>
                  <CopyGenesButton
                    onClick={handleCopyGenes(index)}
                    sdsType="primary"
                    sdsStyle="minimal"
                    isAllCaps={false}
                    startIcon={
                      <Icon sdsIcon="copy" sdsSize="s" sdsType="button" />
                    }
                  >
                    Copy
                  </CopyGenesButton>
                </td>
              </tr>
            </thead>
            <tbody>
              {results.map((result) => {
                const { name: symbol, pValue, effectSize } = result;
                return (
                  <tr key={symbol}>
                    <td>{symbol}</td>
                    <td>{pValue}</td>
                    <td>{effectSize.toPrecision(4)}</td>
                  </tr>
                );
              })}
            </tbody>
          </StyledHTMLTable>
        );
      })}
    </div>
  );
}
