import { Intent } from "@blueprintjs/core";
import { LoadingIndicator } from "czifui";
import React, { useCallback, useContext, useMemo } from "react";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  OntologyTerm,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import Toast from "src/views/Collection/components/Toast";
import { DispatchContext, StateContext } from "../../common/store";
import { selectGenes, selectTissues } from "../../common/store/actions";
import { Gene } from "../../common/types";
import Organism from "./components/Organism";
import QuickSelect from "./components/QuickSelect";
import { ActionWrapper, Container, LoadingIndicatorWrapper } from "./style";

interface Tissue {
  name: string;
}

export default function GeneSearchBar(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedGenes, selectedTissues, selectedOrganismId } =
    useContext(StateContext);

  const { data, isLoading } = usePrimaryFilterDimensions();

  const { genes: rawGenes, tissues } = data || {};

  const genes: Gene[] = useMemo(() => {
    if (!rawGenes) return [];

    return rawGenes[selectedOrganismId || ""] || [];
  }, [rawGenes, selectedOrganismId]);

  const genesByName = useMemo(() => {
    return genes.reduce((acc, gene) => {
      return acc.set(gene.name, gene);
    }, new Map<Gene["name"], Gene>());
  }, [genes]);

  const tissuesByName = useMemo(() => {
    const result = new Map<string, Tissue>();

    if (!tissues) return new Map<string, Tissue>();

    Object.values(tissues).forEach((tissueGroup) =>
      tissueGroup.reduce((acc, tissue) => {
        return acc.set(tissue.name, tissue);
      }, result)
    );

    return result;
  }, [tissues]);

  const flattenedTissues = useMemo((): Array<OntologyTerm> => {
    if (!tissues) return [];
    return Object.values(tissues).reduce((acc, tissueGroup) => {
      return acc.concat(
        // (thuang): Product requirement to exclude "cell culture" from the list
        // https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/2335
        tissueGroup.filter((tissue) => !tissue.name.includes("(cell culture)"))
      );
    }, new Array<OntologyTerm>());
  }, [tissues]);

  const selectedTissueOptions: Tissue[] = useMemo(() => {
    return selectedTissues.map((tissue: string) => {
      return tissuesByName.get(tissue) as Tissue;
    });
  }, [selectedTissues, tissuesByName]);

  const selectedGeneOptions: Gene[] = useMemo(() => {
    return selectedGenes.map((gene: string) => {
      return genesByName.get(gene) as Gene;
    });
  }, [selectedGenes, genesByName]);

  const handleGeneNotFound = useCallback((geneName: string): void => {
    Toast.show({
      intent: Intent.DANGER,
      message: `Gene not found: ${geneName}`,
    });
  }, []);

  return (
    <Container>
      <ActionWrapper>
        <Organism isLoading={isLoading} />

        <QuickSelect
          items={flattenedTissues || EMPTY_ARRAY}
          itemsByName={tissuesByName}
          multiple
          selected={selectedTissueOptions}
          setSelected={handleSelectTissues}
          label="Add Tissue"
          dataTestId="add-tissue"
          placeholder="Search"
          isLoading={isLoading}
          analyticsEvent={EVENTS.WMG_SELECT_TISSUE}
        />

        <QuickSelect
          items={genes}
          itemsByName={genesByName}
          selected={selectedGeneOptions}
          multiple
          setSelected={handleSelectGenes}
          onItemNotFound={handleGeneNotFound}
          label="Add Gene"
          dataTestId="add-gene"
          placeholder="Search or paste comma separated gene names"
          isLoading={isLoading}
          analyticsEvent={EVENTS.WMG_SELECT_GENE}
        />

        {isLoading && (
          <LoadingIndicatorWrapper>
            <LoadingIndicator sdsStyle="tag" />
          </LoadingIndicatorWrapper>
        )}
      </ActionWrapper>
    </Container>
  );

  function handleSelectTissues(tissues: Tissue[]) {
    if (!dispatch) return;

    dispatch(selectTissues(tissues.map((tissue) => tissue.name)));
  }

  function handleSelectGenes(genes: Gene[]) {
    if (!dispatch) return;

    dispatch(selectGenes(genes.map((gene) => gene.name)));
  }
}
