import { Intent } from "@blueprintjs/core";
import { DefaultMenuSelectOption } from "czifui";
import Head from "next/head";
import { useCallback, useEffect, useMemo, useReducer, useState } from "react";
import { API } from "src/common/API";
import { EMPTY_OBJECT } from "src/common/constants/utils";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import SideBar from "src/components/common/SideBar";
import { Position } from "src/components/common/SideBar/style";
import { API_URL } from "src/configs/configs";
import Toast from "../Collection/components/Toast";
import { View } from "../globalStyle";
import {
  DispatchContext,
  INITIAL_STATE,
  reducer,
  StateContext,
} from "./common/store";
import {
  deleteSelectedGenesAndSelectedCellTypeIds,
  selectGenes,
  tissueCellTypesFetched,
} from "./common/store/actions";
import {
  CellTypeSummary,
  Filters as IFilters,
  Gene,
  GeneExpressionSummary,
  Tissue,
} from "./common/types";
import Beta from "./components/Beta";
import Filters from "./components/Filters";
import GeneFetcher from "./components/GeneFetcher";
import GeneSearchBar from "./components/GeneSearchBar";
import GetStarted from "./components/GetStarted";
import HeatMap from "./components/HeatMap";
import InfoPanel from "./components/InfoPanel";
import { SideBarPositioner, SideBarWrapper, Top, Wrapper } from "./style";

const INFO_PANEL_WIDTH_PX = 320;

const WheresMyGene = (): JSX.Element => {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);
  const [filters, setFilters] = useState<IFilters>(EMPTY_OBJECT);

  const { selectedGenes, selectedCellTypeIds } = state;

  /**
   * This holds ALL the geneData we have loaded from the API, including previously
   * and currently selected genes.
   * We use `selectedGeneData` to subset the data to only the genes that are
   * currently selected.
   */
  const [geneExpressionSummaries, setGeneExpressionSummaries] =
    useState<{ [geneName: string]: GeneExpressionSummary }>(EMPTY_OBJECT);

  /**
   * This holds ALL the cell type data we have loaded from the API
   */
  const [cellTypes, setCellTypes] =
    useState<{ [tissue: Tissue]: CellTypeSummary[] }>(EMPTY_OBJECT);

  /**
   * This holds only the CellTypeSummary objects that are currently selected in
   * `state.selectedCellTypeIds`.
   */
  const selectedCellTypes = useMemo(() => {
    const result: { [tissue: Tissue]: CellTypeSummary[] } = {};

    for (const [tissue, selectedIds] of Object.entries(selectedCellTypeIds)) {
      for (const selectedId of selectedIds) {
        const cellType = cellTypes[tissue].find(
          (cellType) => cellType.id === selectedId
        );

        if (cellType !== undefined) {
          const tissueCellTypes = result[tissue] || [];
          tissueCellTypes.push({ ...cellType, tissue });
          result[tissue] = tissueCellTypes;
        }
      }
    }

    return result;
  }, [selectedCellTypeIds, cellTypes]);

  /**
   * This indicates which tissues have less cell types than the API response,
   * indicating the user has deleted some cell types manually
   */
  const tissuesWithDeletedCellTypes = useMemo(() => {
    const result = [];

    for (const [tissue, tissueCellTypes] of Object.entries(cellTypes)) {
      if (selectedCellTypeIds[tissue]?.length < tissueCellTypes.length) {
        result.push(tissue);
      }
    }

    return result;
  }, [cellTypes, selectedCellTypeIds]);

  const selectedGeneData = useMemo(() => {
    return selectedGenes.map((geneName) => {
      return geneExpressionSummaries[geneName];
    });
  }, [selectedGenes, geneExpressionSummaries]);

  useEffect(() => {
    fetchCellTypes();

    async function fetchCellTypes(): Promise<void> {
      const response = await fetch(
        API_URL + API.WMG_CELL_TYPES,
        DEFAULT_FETCH_OPTIONS
      );

      // DEBUG
      // DEBUG
      // DEBUG
      // (thuang): Local test data
      // const response = await fetch(
      //   "https://wmg-prototype-data-dev-public.s3.amazonaws.com/lung-tissue-10x-human/lung_tissue_cell_types.json"
      // );

      const cellTypes = await response.json();

      // (thuang): Table y-axis defaults to descending order
      // .reverse() mutates the original array
      const cellTypeSummaries = cellTypes.reverse();

      setCellTypes({ lung: cellTypeSummaries });
      dispatch(tissueCellTypesFetched("lung", cellTypeSummaries));
    }
  }, [setCellTypes]);

  const handleGenesOnchange = useCallback(
    (genes: Gene[]) => {
      dispatch(selectGenes(genes.map((gene) => gene.name)));
    },
    [dispatch]
  );

  // Listen to delete keyboard press event
  useEffect(() => {
    document.addEventListener("keydown", handleKeyDown);

    return () => {
      document.removeEventListener("keydown", handleKeyDown);
    };

    function handleKeyDown(event: KeyboardEvent): void {
      if (event.code === "Backspace") {
        dispatch(deleteSelectedGenesAndSelectedCellTypeIds());
      }
    }
  }, [dispatch]);

  const handleFiltersChange = useCallback(
    (key: keyof IFilters, options: DefaultMenuSelectOption[] | null) => {
      setFilters((currentFilters) => ({
        ...currentFilters,
        [key]: options,
      }));
    },
    []
  );

  const hasSelectedGeneData = selectedGeneData.length > 0;
  const hasSelectedCellTypeIds = Object.keys(selectedCellTypeIds).length > 0;

  const shouldShowHeatMap = useMemo(() => {
    return hasSelectedGeneData && hasSelectedCellTypeIds;
  }, [hasSelectedGeneData, hasSelectedCellTypeIds]);

  return (
    <DispatchContext.Provider value={dispatch}>
      <StateContext.Provider value={state}>
        <Head>
          <title>cellxgene | Where&apos;s My Gene</title>
        </Head>

        <SideBar
          label="Filters"
          isOpen
          SideBarWrapperComponent={SideBarWrapper}
          SideBarPositionerComponent={SideBarPositioner}
        >
          <Filters filters={filters} onFiltersChange={handleFiltersChange} />
        </SideBar>

        <SideBar
          width={INFO_PANEL_WIDTH_PX}
          label="Info"
          isOpen
          position={Position.RIGHT}
          SideBarWrapperComponent={SideBarWrapper}
          SideBarPositionerComponent={SideBarPositioner}
        >
          <InfoPanel />
        </SideBar>

        <View hideOverflow>
          <Wrapper>
            <Top>
              <GeneSearchBar onGenesChange={handleGenesOnchange} />
              <Beta />
            </Top>

            {shouldShowHeatMap ? (
              <HeatMap
                cellTypes={selectedCellTypes}
                genes={selectedGenes}
                selectedGeneData={selectedGeneData}
                tissuesWithDeletedCellTypes={tissuesWithDeletedCellTypes}
                allTissueCellTypes={cellTypes}
              />
            ) : (
              <GetStarted />
            )}

            {selectedGenes.map((selectedGene) => {
              return (
                <GeneFetcher
                  name={selectedGene}
                  key={selectedGene}
                  onSuccess={handleGeneFetchSuccess}
                  onError={handleGeneFetchError}
                />
              );
            })}
          </Wrapper>
        </View>
      </StateContext.Provider>
    </DispatchContext.Provider>
  );

  function handleGeneFetchSuccess(geneData: GeneExpressionSummary) {
    setGeneExpressionSummaries((latestGeneData) => {
      return {
        ...latestGeneData,
        [geneData.name]: geneData,
      };
    });
  }

  function handleGeneFetchError(name: Gene["name"]) {
    Toast.show({
      intent: Intent.DANGER,
      message: `No data available for gene: ${name}`,
    });
  }
};

export default WheresMyGene;
