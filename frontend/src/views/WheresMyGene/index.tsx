import { Intent } from "@blueprintjs/core";
import cloneDeep from "lodash/cloneDeep";
import Head from "next/head";
import { useMemo, useState } from "react";
import { Cell as ICell, Column } from "react-table";
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
import TreeTable from "./components/TreeTable";
import CELL_TYPES from "./mocks/lung_tissue_cell_types.json";
import { Wrapper } from "./style";

const WheresMyGene = (): JSX.Element => {
  const [genes, setGenes] = useState<Gene[]>(EMPTY_ARRAY);
  const [geneData, setGeneData] = useState<RawGeneExpression[]>(EMPTY_ARRAY);

  const data = useMemo(
    () => integrateCelTypesAndGenes(CELL_TYPES as CellTypeAndGenes[], geneData),
    [geneData]
  );

  const columns: Column<CellTypeAndGenes>[] = useMemo(() => {
    return [
      {
        Cell({ row }: ICell) {
          const {
            values: { name },
          } = row;

          return <span>{name}</span>;
        },
        Header: "",
        accessor: "name",
        minWidth: 200,
      },
      ...genes.map(({ name }) => {
        return {
          Header: name,
          accessor: `expressions.${name}` as unknown as "expressions",
          isLoading: !geneData.find((gene) => {
            return gene.gene_name === name;
          }),
          maxWidth: 20,
        };
      }),
    ];
  }, [genes, geneData]);

  return (
    <>
      <Head>
        <title>cellxgene | Where&apos;s My Gene</title>
      </Head>

      <Wrapper>
        <GeneSearchBar onGenesChange={setGenes} />

        <br />
        <br />

        <TreeTable columns={columns} data={data} />

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
