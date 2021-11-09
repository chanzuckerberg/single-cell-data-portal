import cloneDeep from "lodash/cloneDeep";
import Head from "next/head";
import { useMemo } from "react";
import GeneSearchBar from "./components/GeneSearchBar";
import TreeTable from "./components/TreeTable";
import { GENES, TISSUE } from "./mocks/brain";

const data = integrateTissuesAndGenes(TISSUE);

const WheresMyGene = (): JSX.Element => {
  const columns = useMemo(
    () => [
      {
        Cell({ row }) {
          const {
            values: { name },
          } = row;

          return <span>{name}</span>;
        },
        Header: "",
        accessor: "name",
        minWidth: 200,
      },
      ...GENES.map(({ name }) => {
        return {
          Header: name,
          accessor: name,
          maxWidth: 20,
        };
      }),
    ],
    []
  );

  return (
    <>
      <Head>
        <title>cellxgene | Where&apos;s My Gene</title>
      </Head>

      <GeneSearchBar />

      <TreeTable columns={columns} data={data} />
    </>
  );
};

interface Node {
  id: string;
  parentId?: string;
  name: string;
}

function integrateTissuesAndGenes(nodes: Node[]) {
  const geneMaps = GENES.map((gene) => rawGeneDataToMap(gene));

  return cloneDeep(nodes).map((node) => {
    const { id } = node;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(id);

      if (columnData !== undefined) {
        node[name] = columnData;
      }
    }

    return node;
  });
}

function rawGeneDataToMap(gene) {
  const { data, name } = gene;

  return [name, new Map(data.map((row) => [row.id, row]))];
}

export default WheresMyGene;
