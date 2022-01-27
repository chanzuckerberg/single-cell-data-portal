import { Intent } from "@blueprintjs/core";
import cloneDeep from "lodash/cloneDeep";
import debounce from "lodash/debounce";
import Head from "next/head";
import { useEffect, useMemo, useState } from "react";
import { API } from "src/common/API";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { API_URL } from "src/configs/configs";
import Toast from "../Collection/components/Toast";
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
  /**
   * This is the genes that are currently selected.
   */
  const [genes, setGenes] = useState<Gene[]>(EMPTY_ARRAY);

  /**
   * This holds ALL the geneData we have loaded from the API, including previously
   * and currently selected genes.
   * We use `selectedGeneData` to subset the data to only the genes that are
   * currently selected.
   */
  const [geneData, setGeneData] =
    useState<GeneExpressionSummary[]>(EMPTY_ARRAY);
  const [cellTypes, setCellTypes] = useState<CellTypeSummary[]>(EMPTY_ARRAY);

  /**
   * This is the formatted data that we use to render the heatmap.
   */
  const [data, setData] = useState<CellTypeSummary[]>(EMPTY_ARRAY);

  const selectedGeneData = useMemo(() => {
    return geneData.filter((geneData) =>
      genes.some((gene) => gene.name === geneData.name)
    );
  }, [genes, geneData]);

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
      setCellTypes(cellTypes.reverse());
    }
  }, []);

  const debouncedIntegrateCellTypesAndGenes = useMemo(() => {
    return debounce(
      (cellTypes, geneData) => {
        setData(integrateCelTypesAndGenes(cellTypes, geneData));
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, []);

  /**
   * Performance optimization:
   * We only format and `setData()` after the watch list has stopped changing for
   * `DEBOUNCE_MS`
   */
  useEffect(() => {
    debouncedIntegrateCellTypesAndGenes(cellTypes, selectedGeneData);
  }, [selectedGeneData, cellTypes, debouncedIntegrateCellTypesAndGenes]);

  return (
    <>
      <Head>
        <title>cellxgene | Where&apos;s My Gene</title>
      </Head>

      <Wrapper>
        <GeneSearchBar onGenesChange={setGenes} />

        <HeatMap cellTypes={cellTypes} data={data} genes={genes} />

        {genes.map((gene) => {
          const { name } = gene;

          return (
            <GeneFetcher
              fetchedGenes={genes}
              name={name}
              key={name}
              onSuccess={handleGeneFetchSuccess}
              onError={handleGeneFetchError}
            />
          );
        })}
      </Wrapper>
    </>
  );

  function handleGeneFetchSuccess(geneData: GeneExpressionSummary) {
    setGeneData((latestGeneData) => [...latestGeneData, geneData]);
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
