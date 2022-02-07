import { Intent } from "@blueprintjs/core";
import cloneDeep from "lodash/cloneDeep";
import debounce from "lodash/debounce";
import Head from "next/head";
import { useCallback, useEffect, useMemo, useReducer, useState } from "react";
import { API } from "src/common/API";
import { EMPTY_ARRAY } from "src/common/constants/utils";
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
import { selectGenes } from "./common/store/actions";
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

  const { selectedGenes } = state;

  /**
   * This holds ALL the geneData we have loaded from the API, including previously
   * and currently selected genes.
   * We use `selectedGeneData` to subset the data to only the genes that are
   * currently selected.
   */
  const [geneExpressionSummaries, setGeneExpressionSummaries] =
    useState<GeneExpressionSummary[]>(EMPTY_ARRAY);

  const [cellTypes, setCellTypes] = useState<CellTypeSummary[]>(EMPTY_ARRAY);

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
      setCellTypes(cellTypes.reverse());
    }
  }, []);

  const debouncedIntegrateCellTypesAndGenes = useMemo(() => {
    return debounce(
      (cellTypes, geneData) => {
        setCellTypeSummaries(integrateCelTypesAndGenes(cellTypes, geneData));
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, []);

  /**
   * Performance optimization:
   * We only format and `setCellTypeSummaries()` after the watch list has stopped changing for
   * `DEBOUNCE_MS`
   */
  useEffect(() => {
    debouncedIntegrateCellTypesAndGenes(cellTypes, selectedGeneData);
  }, [selectedGeneData, cellTypes, debouncedIntegrateCellTypesAndGenes]);

  const handleGenesOnchange = useCallback(
    (genes: Gene[]) => {
      dispatch(selectGenes(genes.map((gene) => gene.name)));
    },
    [dispatch]
  );

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
              cellTypes={cellTypes}
              data={cellTypeSummaries}
              genes={selectedGenes}
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

function integrateCelTypesAndGenes(
  nodes: CellTypeSummary[],
  genes: GeneExpressionSummary[]
): CellTypeSummary[] {
  const geneMaps = genes.map((gene) => rawGeneDataToMap(gene));

  return cloneDeep(nodes).map((node) => {
    const { id } = node;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(id);

      if (columnData !== undefined) {
        node.geneExpressions = {
          ...(node.geneExpressions || {}),
          [name]: columnData,
        };
      }
    }

    return node;
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
