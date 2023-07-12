import { Intent } from "@blueprintjs/core";
import { LoadingIndicator } from "@czi-sds/components";
import React, {
  useCallback,
  useContext,
  useMemo,
  useState,
  useEffect,
} from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  usePrimaryFilterDimensions,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import Toast from "src/views/Collection/components/Toast";
import { DispatchContext, StateContext } from "../../common/store";
import {
  deleteAllGenes,
  selectGenes,
  selectTissues,
} from "../../common/store/actions";
import { Gene } from "../../common/types";
import QuickSelect from "./components/QuickSelect";
import {
  ActionWrapper,
  Container,
  LoadingIndicatorWrapper,
  StyledButtonWrapper,
  StyledClearButton,
} from "./style";

interface Tissue {
  name: string;
}

export default function GeneSearchBar({
  className,
}: {
  className?: string;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedGenes, selectedTissues, selectedOrganismId } =
    useContext(StateContext);

  const { data, isLoading: isLoadingPrimaryFilters } =
    usePrimaryFilterDimensions();
  const { data: filterData, isLoading: isLoadingFilters } =
    useFilterDimensions();
  const isLoading = isLoadingPrimaryFilters || isLoadingFilters;

  const { tissue_terms: filteredTissues } = filterData;
  const { genes: rawGenes, tissues: rawTissues } = data || {};

  const genes: Gene[] = useMemo(() => {
    if (!rawGenes) return [];

    return rawGenes[selectedOrganismId || ""] || [];
  }, [rawGenes, selectedOrganismId]);

  const [tissues, setTissues] = useState<Tissue[]>([]);

  useEffect(() => {
    if (rawTissues && filteredTissues.length) {
      const temp = rawTissues[selectedOrganismId || ""] || [];
      // (thuang): Product requirement to exclude "cell culture" from the list
      // https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/2335
      const newTissues = temp.filter((tissue) => {
        const notCellCulture = !tissue.name.includes("(cell culture)");
        const notFiltered =
          filteredTissues.map((val) => val.name).includes(tissue.name) ||
          selectedTissues?.includes(tissue.name);
        return notCellCulture && notFiltered;
      });
      setTissues(newTissues);
    }
  }, [
    rawTissues,
    filteredTissues,
    selectedOrganismId,
    setTissues,
    selectedTissues,
  ]);

  /**
   * NOTE: key is gene name in lowercase
   */
  const genesByName = useMemo(() => {
    return genes.reduce((acc, gene) => {
      return acc.set(gene.name.toLowerCase(), gene);
    }, new Map<Gene["name"], Gene>());
  }, [genes]);

  /**
   * NOTE: key is tissue name in lowercase
   */
  const tissuesByName = useMemo(() => {
    return tissues.reduce((acc, tissue) => {
      return acc.set(tissue.name.toLowerCase(), tissue);
    }, new Map<Tissue["name"], Tissue>());
  }, [tissues]);

  const selectedTissueOptions: Tissue[] = useMemo(() => {
    return (
      selectedTissues?.map((tissue: string) => {
        return tissuesByName.get(tissue.toLowerCase()) as Tissue;
      }) || []
    );
  }, [selectedTissues, tissuesByName]);

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
          items={tissues}
          itemsByName={tissuesByName}
          multiple
          selected={selectedTissueOptions}
          setSelected={handleSelectTissues}
          label=""
          text="Tissue"
          dataTestId="add-tissue-btn"
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
          label=""
          text="Gene"
          dataTestId="add-gene-btn"
          placeholder="Search or paste comma separated gene names"
          isLoading={isLoadingPrimaryFilters}
          analyticsEvent={EVENTS.WMG_SELECT_GENE}
        />

        {/* Clear Genes button */}
        {!selectedGenes.length || (
          <StyledButtonWrapper
            onClick={() => {
              if (dispatch) {
                track(EVENTS.WMG_CLEAR_GENES_CLICKED);

                dispatch(deleteAllGenes());
              }
            }}
          >
            <StyledClearButton
              data-testid="clear-genes-button"
              sdsType="primary"
              sdsStyle="minimal"
              isAllCaps={false}
            >
              Clear Genes
            </StyledClearButton>
          </StyledButtonWrapper>
        )}

        {isLoadingPrimaryFilters && (
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
