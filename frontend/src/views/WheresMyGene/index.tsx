import { Intent } from "@blueprintjs/core";
import cloneDeep from "lodash/cloneDeep";
import Head from "next/head";
import { useEffect, useMemo, useState } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import Toast from "../Collection/components/Toast";
import {
  CellTypeAndGenes,
  Gene,
  GeneExpression,
  RawGeneExpression,
} from "./common/types";
import GeneFetcher from "./components/GeneFetcher";
import GeneSearchBar from "./components/GeneSearchBar";
import HeatMap from "./components/HeatMap";
import { Wrapper } from "./style";

const WheresMyGene = (): JSX.Element => {
  const [genes, setGenes] = useState<Gene[]>(EMPTY_ARRAY);
  const [geneData, setGeneData] = useState<RawGeneExpression[]>(EMPTY_ARRAY);
  const [cellTypes, setCellTypes] = useState<CellTypeAndGenes[]>(EMPTY_ARRAY);

  useEffect(() => {
    fetchCellTypes();

    async function fetchCellTypes(): Promise<void> {
      // const response = await fetch(
      //   API_URL + API.WMG_CELL_TYPES,
      //   DEFAULT_FETCH_OPTIONS
      // );

      const response = await fetch(
        "https://wmg-prototype-data-dev-public.s3.amazonaws.com/lung-tissue-10x-human/lung_tissue_cell_types.json"
      );

      const cellTypes = await response.json();

      setCellTypes(cellTypes);
    }
  }, []);

  const data = useMemo(
    () => integrateCelTypesAndGenes(cellTypes, geneData),
    [cellTypes, geneData]
  );

  return (
    <>
      <Head>
        <title>cellxgene | Where&apos;s My Gene</title>
      </Head>

      <Wrapper>
        {/* <TreeTable columns={columns} data={data} /> */}

        <HeatMap cellTypes={cellTypes} data={data} genes={genes} />

        <GeneSearchBar onGenesChange={setGenes} />

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

  function handleGeneFetchSuccess(geneData: RawGeneExpression) {
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
  nodes: CellTypeAndGenes[],
  genes: RawGeneExpression[]
): CellTypeAndGenes[] {
  const geneMaps = genes.map((gene) => rawGeneDataToMap(gene));

  return cloneDeep(nodes).map((node) => {
    const { id } = node;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(id);

      if (columnData !== undefined) {
        node.expressions = {
          ...(node.expressions || {}),
          [name]: columnData,
        };
      }
    }

    return node;
  });
}

function rawGeneDataToMap(
  gene: RawGeneExpression
): [string, Map<string, GeneExpression>] {
  const { cell_types, gene_name } = gene;

  return [gene_name, new Map(cell_types?.map((row) => [row.id, row]))];
}

export default WheresMyGene;
