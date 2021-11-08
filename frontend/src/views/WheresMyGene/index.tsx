import cloneDeep from "lodash/cloneDeep";
import Head from "next/head";
import { useMemo } from "react";
import TreeTable from "./components/TreeTable";
import { GENES, TISSUE } from "./mocks/brain";

const data = buildTissueTree(TISSUE).subRows;

const WheresMyGene = (): JSX.Element => {
  const columns = useMemo(
    () => [
      {
        Cell({ row }) {
          // DEBUG
          // DEBUG
          // DEBUG
          console.log("---row", row);

          const {
            depth,
            canExpand,
            isExpanded,
            getToggleRowExpandedProps,
            values: { name },
          } = row;

          return (
            <span
              {...getToggleRowExpandedProps({
                style: {
                  // We can even use the row.depth property
                  // and paddingLeft to indicate the depth
                  // of the row
                  paddingLeft: `${depth * 2}rem`,
                },
              })}
            >
              {canExpand ? (isExpanded ? "👇 " + name : "👉 " + name) : name}
            </span>
          );
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

      <TreeTable columns={columns} data={data} />
    </>
  );

  // function renderRelativeExpressionCell(row) {
  //   const {
  //     data: { id },
  //   } = row;
  //   // DEBUG
  //   // DEBUG
  //   // DEBUG
  //   console.log("row.data", row.data);
  //   console.log("-----gM5741", gM5741);

  //   const rowData = gM5741.get(id);

  //   const { negEntropy, relativeExpression } = rowData;

  //   return (
  //     <div
  //       style={{
  //         height: "100%",
  //         display: "flex",
  //         alignItems: "center",
  //         justifyContent: "center",
  //       }}
  //     >
  //       <Square color={interpolateViridis(negEntropy ?? relativeExpression)} />
  //     </div>
  //   );
  // }
};

interface Node {
  id: string;
  parentId?: string;
  name: string;
}

function buildTissueTree(nodes: Node[]) {
  const idToNodes = new Map();
  let rootNode = null;

  const geneMaps = GENES.map((gene) => rawGeneDataToMap(gene));

  for (const node of cloneDeep(nodes)) {
    const { id, parentId } = node;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(id);

      if (columnData !== undefined) {
        node[name] = columnData;
      }
    }

    if (!rootNode) {
      rootNode = node;
    }

    idToNodes.set(id, node);

    if (parentId) {
      const parentNode = idToNodes.get(parentId);

      if (!parentNode) {
        throw Error(`parentNode id does not exist: ${parentId}`);
      }

      parentNode.subRows = parentNode.subRows || [];
      parentNode.subRows.push(node);
    }
  }

  // DEBUG
  // DEBUG
  // DEBUG
  // Don't mutate the original data!!!!
  console.log("----><><><><><><>rootNode", rootNode);

  return rootNode;
}

function rawGeneDataToMap(gene) {
  const { data, name } = gene;

  return [name, new Map(data.map((row) => [row.id, row]))];
}

export default WheresMyGene;
