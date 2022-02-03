import { Intent } from "@blueprintjs/core";
import cloneDeep from "lodash/cloneDeep";
import debounce from "lodash/debounce";
import Head from "next/head";
import { useCallback, useEffect, useMemo, useReducer, useState } from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import SideBar from "src/components/common/SideBar";
import { Position } from "src/components/common/SideBar/style";
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
  CellTypeGeneExpressionSummaryData,
  CellTypeSummary,
  Gene,
  GeneExpressionSummary,
} from "./common/types";
import GeneFetcher from "./components/GeneFetcher";
import GeneSearchBar from "./components/GeneSearchBar";
import HeatMap from "./components/HeatMap";
import { Wrapper } from "./style";

const DEBOUNCE_MS = 2 * 1000;

const WheresMyGene = (): JSX.Element => {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);

  const { selectedGenes, selectedCellTypeIds } = state;

  /**
   * This holds ALL the geneData we have loaded from the API, including previously
   * and currently selected genes.
   * We use `selectedGeneData` to subset the data to only the genes that are
   * currently selected.
   */
  const [geneExpressionSummaries, setGeneExpressionSummaries] =
    useState<GeneExpressionSummary[]>(EMPTY_ARRAY);

  /**
   * This holds ALL the cell type data we have loaded from the API
   */
  const [cellTypes, setCellTypes] =
    useState<{ [tissue: string]: CellTypeSummary[] }>(EMPTY_OBJECT);

  /**
   * This holds only the CellTypeSummary objects that are currently selected in
   * `state.selectedCellTypeIds`.
   * NOTE: We also prepend corresponding tissues for their cell types for
   * rendering in the heat map.
   */
  const selectedCellTypes = useMemo(() => {
    const result = [];

    for (const [tissue, selectedIds] of Object.entries(selectedCellTypeIds)) {
      for (const selectedId of selectedIds) {
        const cellType = cellTypes[tissue].find(
          (cellType) => cellType.id === selectedId
        );

        if (cellType !== undefined) {
          result.push({ ...cellType, tissue });
        }
      }
      // (thuang): Tissue needs to be last in the list
      result.push({ id: tissue, name: tissue, tissue });
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

  /**
   * This is the formatted data that we use to render the heatmap.
   */
  const [cellTypeSummaries, setCellTypeSummaries] =
    useState<CellTypeSummary[]>(EMPTY_ARRAY);

  const selectedGeneData = useMemo(() => {
    return geneExpressionSummaries.filter((geneExpressionSummary) =>
      selectedGenes.some(
        (selectedGene) => selectedGene === geneExpressionSummary.name
      )
    );
  }, [selectedGenes, geneExpressionSummaries]);

  useEffect(() => {
    fetchCellTypes();

    async function fetchCellTypes(): Promise<void> {
      // const response = await fetch(
      //   API_URL + API.WMG_CELL_TYPES,
      //   DEFAULT_FETCH_OPTIONS
      // );

      // DEBUG
      // DEBUG
      // DEBUG
      // (thuang): Local test data
      const response = await fetch(
        "https://wmg-prototype-data-dev-public.s3.amazonaws.com/lung-tissue-10x-human/lung_tissue_cell_types.json"
      );

      const cellTypes = await response.json();

      // (thuang): Table y-axis defaults to descending order
      // .reverse() mutates the original array
      const cellTypeSummaries = cellTypes.reverse();
      setCellTypes({ lung: cellTypeSummaries });
      dispatch(tissueCellTypesFetched("lung", cellTypeSummaries));
    }
  }, [setCellTypes]);

  const debouncedIntegrateCellTypesAndGenes = useMemo(() => {
    return debounce(
      (cellTypes, geneData) => {
        setCellTypeSummaries(integrateCelTypesAndGenes(cellTypes, geneData));
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, []);

  // Cancel debounce when unmounting
  useEffect(() => {
    return () => debouncedIntegrateCellTypesAndGenes.cancel();
  }, [debouncedIntegrateCellTypesAndGenes]);

  /**
   * Performance optimization:
   * We only format and `setCellTypeSummaries()` after the watch list has stopped changing for
   * `DEBOUNCE_MS`
   */
  useEffect(() => {
    debouncedIntegrateCellTypesAndGenes(selectedCellTypes, selectedGeneData);
  }, [
    selectedGeneData,
    selectedCellTypes,
    debouncedIntegrateCellTypesAndGenes,
  ]);

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

  return (
    <DispatchContext.Provider value={dispatch}>
      <StateContext.Provider value={state}>
        <Head>
          <title>cellxgene | Where&apos;s My Gene</title>
        </Head>

        <SideBar label="Filters" isOpen>
          <span>
            Lorem ipsum dolor sit amet consectetur adipisicing elit. Doloribus
            autem deserunt assumenda repudiandae repellat quis sunt quae, aut
            vero rem itaque labore praesentium iure exercitationem minus iste
            laudantium sed aliquid.
          </span>
        </SideBar>

        <SideBar label="Filters" isOpen position={Position.RIGHT}>
          <span>
            Lorem ipsum dolor sit amet consectetur adipisicing elit. Doloribus
            autem deserunt assumenda repudiandae repellat quis sunt quae, aut
            vero rem itaque labore praesentium iure exercitationem minus iste
            laudantium sed aliquid.
          </span>
        </SideBar>

        <View hideOverflow>
          <Wrapper>
            <GeneSearchBar onGenesChange={handleGenesOnchange} />

            <HeatMap
              cellTypes={selectedCellTypes}
              data={cellTypeSummaries}
              genes={selectedGenes}
              tissuesWithDeletedCellTypes={tissuesWithDeletedCellTypes}
              allTissueCellTypes={cellTypes}
            />

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
    setGeneExpressionSummaries((latestGeneData) => [
      ...latestGeneData,
      geneData,
    ]);
  }

  function handleGeneFetchError(name: Gene["name"]) {
    Toast.show({
      intent: Intent.DANGER,
      message: `No data available for gene: ${name}`,
    });
  }
};

/**
 * Adds gene expressions to the cell types.
 */
function integrateCelTypesAndGenes(
  cellTypeSummaries: CellTypeSummary[],
  geneExpressionSummaries: GeneExpressionSummary[]
): CellTypeSummary[] {
  const geneMaps = geneExpressionSummaries.map((geneExpressionSummary) =>
    rawGeneDataToMap(geneExpressionSummary)
  );

  return cloneDeep(cellTypeSummaries).map((cellTypeSummary) => {
    const { id } = cellTypeSummary;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(id);

      if (columnData !== undefined) {
        cellTypeSummary.geneExpressions = {
          ...(cellTypeSummary.geneExpressions || {}),
          [name]: columnData,
        };
      }
    }

    return cellTypeSummary;
  });
}

function rawGeneDataToMap(
  gene: GeneExpressionSummary
): [string, Map<string, CellTypeGeneExpressionSummaryData>] {
  const { cellTypeGeneExpressionSummaries, name } = gene;

  return [
    name,
    new Map(cellTypeGeneExpressionSummaries?.map((row) => [row.id, row])),
  ];
}

export default WheresMyGene;
