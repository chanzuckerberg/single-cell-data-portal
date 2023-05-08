import { Intent } from "@blueprintjs/core";
import { LoadingIndicator } from "czifui";
import React, { useCallback, useContext, useMemo } from "react";
import { EVENTS } from "src/common/analytics/events";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGeneV2";
import Toast from "src/views/Collection/components/Toast";
import { DispatchContext, StateContext } from "../../common/store";
import { selectGenes } from "../../common/store/actions";
import { Gene } from "../../common/types";
import QuickSelect from "./components/QuickSelect";
import { ActionWrapper, Container, LoadingIndicatorWrapper } from "./style";

export default function GeneSearchBar({
  className,
}: {
  className?: string;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedGenes, selectedOrganismId } = useContext(StateContext);

  const { data, isLoading: isLoadingPrimaryFilters } =
    usePrimaryFilterDimensions();

  const { genes: rawGenes } = data || {};

  const genes: Gene[] = useMemo(() => {
    if (!rawGenes) return [];

    return rawGenes[selectedOrganismId || ""] || [];
  }, [rawGenes, selectedOrganismId]);

  /**
   * NOTE: key is gene name in lowercase
   */
  const genesByName = useMemo(() => {
    return genes.reduce((acc, gene) => {
      return acc.set(gene.name.toLowerCase(), gene);
    }, new Map<Gene["name"], Gene>());
  }, [genes]);

  const selectedGeneOptions: Gene[] = useMemo(() => {
    return selectedGenes.map((gene: string) => {
      return genesByName.get(gene.toLowerCase()) as Gene;
    });
  }, [selectedGenes, genesByName]);

  const handleGeneNotFound = useCallback((geneName: string): void => {
    Toast.show({
      intent: Intent.DANGER,
      message: `Gene not found: ${geneName}`,
    });
  }, []);

  return (
    <Container {...{ className }}>
      <ActionWrapper>
        <QuickSelect
          items={genes}
          itemsByName={genesByName}
          selected={selectedGeneOptions}
          multiple
          setSelected={handleSelectGenes}
          onItemNotFound={handleGeneNotFound}
          label=""
          text="Gene"
          dataTestId="add-gene-btn"
          placeholder="Search or paste comma separated gene names"
          isLoading={isLoadingPrimaryFilters}
          analyticsEvent={EVENTS.WMG_SELECT_GENE}
        />

        {isLoadingPrimaryFilters && (
          <LoadingIndicatorWrapper>
            <LoadingIndicator sdsStyle="tag" />
          </LoadingIndicatorWrapper>
        )}
      </ActionWrapper>
    </Container>
  );

  function handleSelectGenes(genes: Gene[]) {
    if (!dispatch) return;

    dispatch(selectGenes(genes.map((gene) => gene.name)));
  }
}
