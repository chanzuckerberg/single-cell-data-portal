import { Intent } from "@blueprintjs/core";
import { LoadingIndicator } from "@czi-sds/components";
import React, { useCallback, useContext, useMemo } from "react";
import { EVENTS } from "src/common/analytics/events";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import Toast from "src/views/Collection/components/Toast";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  deleteAllGenes,
  selectGenes,
} from "src/views/WheresMyGene/common/store/actions";
import { Gene } from "src/views/WheresMyGene/common/types";
import QuickSelect from "./components/QuickSelect";
import GeneAutocomplete from "./components/GeneAutocomplete";
import { ActionWrapper, Container, LoadingIndicatorWrapper } from "./style";
import {
  StyledButtonWrapper,
  StyledClearButton,
} from "src/views/WheresMyGene/components/GeneSearchBar/style";
import { track } from "src/common/analytics";
import { StyledAutocomplete } from "src/views/WheresMyGene/components/HeatMap/style";
import { InputAdornment, TextField } from "@mui/material";
import SearchIcon from "@mui/icons-material/Search";
import { FixedSizeList, ListChildComponentProps } from "react-window";
import { Style } from "@blueprintjs/icons/lib/esm/generated-icons/16px/paths";

export default function GeneSearchBar({
  className,
}: {
  className?: string;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedGenes, selectedOrganismId } = useContext(StateContext);

  const { data, isLoading: isLoadingPrimaryFilters } =
    usePrimaryFilterDimensions(2); //temp version 2

  const { genes: rawGenes } = data || {};

  const genes: Gene[] = useMemo(() => {
    if (!rawGenes) return [];

    return rawGenes[selectedOrganismId || ""] || [];
  }, [rawGenes, selectedOrganismId]);

  const uniqueGeneNamesSet: Set<string> = new Set(
    genes.map((dictionary) => dictionary.name)
  );

  const uniqueGeneNamesArray: string[] = Array.from(uniqueGeneNamesSet);
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
        {/* <GeneAutocomplete
          items={genes}
          itemsByName={genesByName}
          selected={selectedGeneOptions}
          multiple
          setSelected={handleSelectGenes}
          onItemNotFound={handleGeneNotFound}
          label=""
          text="Caro"
          dataTestId="add-gene-btn"
          placeholder="Search or paste comma separated gene names"
          isLoading={isLoadingPrimaryFilters}
          analyticsEvent={EVENTS.WMG_SELECT_GENE}
        /> */}
        <StyledAutocomplete
          multiple
          renderInput={(params) => (
            <TextField
              {...params}
              InputProps={{
                ...params.InputProps,
                startAdornment: (
                  <InputAdornment position="start">
                    <SearchIcon />
                  </InputAdornment>
                ),
                endAdornment: undefined,
              }}
              placeholder="Add genes"
            ></TextField>
          )}
          options={uniqueGeneNamesArray}
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
              Clear All
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

  function handleSelectGenes(genes: Gene[]) {
    if (!dispatch) return;

    dispatch(selectGenes(genes.map((gene) => gene.name)));
  }
}
