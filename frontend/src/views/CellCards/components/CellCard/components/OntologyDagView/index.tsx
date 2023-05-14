import React, { useMemo, MouseEventHandler } from "react";
import { Group } from "@visx/group";
import { Tree, hierarchy } from "@visx/hierarchy";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { LinkHorizontal } from "@visx/shape";
import { TableTitleWrapper, TableTitle } from "../common/style";
import { ONTOLOGY_SECTION_ID } from "../CellCardSidebar";
import { Zoom } from "@visx/zoom";
import { RectClipPath } from "@visx/clip-path";
import {
  CellOntologyTree as TreeNode,
  useCellOntologyTree,
} from "src/common/queries/cellCards";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";

const plum = "#71248e";
const lightpurple = "#374469";
const white = "#ffffff";
const black = "#000000";
export const background = white;

type HierarchyNode = HierarchyPointNode<TreeNode>;

interface NodeProps {
  node: HierarchyNode;
  handleClick: MouseEventHandler<SVGGElement>;
}

function RootNode({ node, handleClick }: NodeProps) {
  return (
    <Group
      top={node.x}
      left={node.y}
      onClick={handleClick}
      style={{ cursor: "pointer" }}
    >
      <circle r={12} fill={plum} />
      <text
        dy=".33em"
        fontSize={9}
        fontFamily="Arial"
        textAnchor="middle"
        style={{ pointerEvents: "none" }}
        fill={black}
      >
        {node.data.name}
      </text>
    </Group>
  );
}

function ParentNode({ node, handleClick }: NodeProps) {
  const width = 40;
  const height = 20;
  const centerX = -width / 2;
  const centerY = -height / 2;

  return (
    <Group top={node.x} left={node.y}>
      <rect
        height={height}
        width={width}
        y={centerY}
        x={centerX}
        fill={background}
        stroke={black}
        strokeWidth={1}
        onClick={handleClick}
        style={{ cursor: "pointer" }}
      />
      <text
        dy=".33em"
        fontSize={9}
        fontFamily="Arial"
        textAnchor="middle"
        style={{ pointerEvents: "none" }}
        fill={black}
      >
        {node.data.name}
      </text>
    </Group>
  );
}

/** Handles rendering Root, Parent, and other Nodes. */
function Node({ node, handleClick }: NodeProps) {
  const width = 40;
  const height = 20;
  const centerX = -width / 2;
  const centerY = -height / 2;
  const isRoot = node.depth === 0;
  const isParent = !!node.children;

  if (isRoot) return <RootNode node={node} handleClick={handleClick} />;
  if (isParent) return <ParentNode node={node} handleClick={handleClick} />;

  return (
    <Group top={node.x} left={node.y}>
      <rect
        height={height}
        width={width}
        y={centerY}
        x={centerX}
        fill={background}
        stroke={black}
        strokeWidth={1}
        strokeDasharray="2,2"
        strokeOpacity={0.6}
        rx={10}
        onClick={handleClick}
        style={{ cursor: "pointer" }}
      />
      <text
        dy=".33em"
        fontSize={9}
        fontFamily="Arial"
        textAnchor="middle"
        fill={black}
        style={{ pointerEvents: "none" }}
      >
        {node.data.name}
      </text>
    </Group>
  );
}

const defaultMargin = { top: 10, left: 80, right: 80, bottom: 10 };

export type TreeProps = {
  cellTypeId: string;
  width: number;
  height: number;
  margin?: { top: number; right: number; bottom: number; left: number };
};

export default function OntologyDagView({
  cellTypeId,
  width,
  height,
  margin = defaultMargin,
}: TreeProps) {
  const { data: rawTree } = useCellOntologyTree(cellTypeId);

  const router = useRouter();
  const data = useMemo(() => {
    if (!rawTree) return null;
    return hierarchy(rawTree);
  }, [rawTree]);

  const yMax = height - margin.top - margin.bottom;
  const xMax = width - margin.left - margin.right;

  // Customize nodeSize and separation
  const nodeSize = [15, 200]; // Increase width and height for more spacing

  const initialTransform = {
    scaleX: 1,
    scaleY: 1,
    translateX: 80,
    translateY: height / 2,
    skewX: 0,
    skewY: 0,
  };

  return (
    <>
      <TableTitleWrapper id={ONTOLOGY_SECTION_ID}>
        <TableTitle>Cell Ontology</TableTitle>
      </TableTitleWrapper>
      {data ? (
        <Zoom<SVGSVGElement>
          width={width}
          height={height}
          scaleXMin={0.1}
          scaleXMax={4}
          scaleYMin={0.1}
          scaleYMax={4}
          initialTransformMatrix={initialTransform}
        >
          {(zoom) => (
            <div className="relative">
              <svg width={width} height={height} ref={zoom.containerRef}>
                <RectClipPath id="zoom-clip" width={width} height={height} />
                <rect width={width} height={height} rx={14} fill={background} />
                <g transform={zoom.toString()}>
                  <Tree<TreeNode>
                    root={data}
                    size={[yMax, xMax]}
                    nodeSize={nodeSize as [number, number]}
                  >
                    {(tree) => (
                      <Group top={margin.top} left={margin.left}>
                        {tree.links().map((link, i) => (
                          <LinkHorizontal
                            key={`link-${i}`}
                            data={link}
                            stroke={lightpurple}
                            strokeWidth="1"
                            fill="none"
                          />
                        ))}
                        {tree.descendants().map((node, i) => (
                          <Node
                            key={`node-${i}`}
                            node={node}
                            handleClick={() => {
                              router.push(
                                `${ROUTES.CELL_CARDS}/${node.data.id.replace(
                                  ":",
                                  "_"
                                )}`
                              );
                            }}
                          />
                        ))}
                      </Group>
                    )}
                  </Tree>
                </g>
              </svg>
            </div>
          )}
        </Zoom>
      ) : null}
    </>
  );
}
