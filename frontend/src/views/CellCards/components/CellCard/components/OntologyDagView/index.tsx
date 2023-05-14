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

const primaryColor = "#0073FF";
const secondaryColor = "#999999";
const tertiaryColor = "#CCCCCC";
const highlightColor = "#3CB371";
const white = "#ffffff";
const black = "#000000";

type HierarchyNode = HierarchyPointNode<TreeNode>;

interface NodeProps {
  node: HierarchyNode;
  handleClick: MouseEventHandler<SVGGElement>;
  isTargetNode: boolean;
}

function RootNode({ node, handleClick, isTargetNode }: NodeProps) {
  return (
    <Group
      top={node.x}
      left={node.y}
      onClick={handleClick}
      style={{ cursor: "pointer" }}
    >
      <RectOrCircle
        node={node.data}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
      />
      <Text name={node.data.name} />
    </Group>
  );
}

function ParentNode({ node, handleClick, isTargetNode }: NodeProps) {
  return (
    <Group top={node.x} left={node.y}>
      <RectOrCircle
        node={node.data}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
      />
      <Text name={node.data.name} />
    </Group>
  );
}
interface RectOrCircleProps {
  handleClick: MouseEventHandler<SVGGElement>;
  isTargetNode: boolean;
  node: TreeNode;
}

function RectOrCircle({ node, handleClick, isTargetNode }: RectOrCircleProps) {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
  if (isTargetNode) {
    color = highlightColor;
  }
  const size = node.n_cells === 0 ? 4 : 8;
  return node.hasChildren ? (
    <circle
      r={size}
      fill={color}
      onClick={handleClick}
      style={{ cursor: "pointer" }}
    />
  ) : (
    <rect
      height={size * 2}
      width={size * 2}
      y={-size}
      x={-size}
      fill={color}
      onClick={handleClick}
      style={{ cursor: "pointer" }}
    />
  );
}

interface TextProps {
  name: string;
}
function Text({ name }: TextProps) {
  return (
    <text
      dy=".33em"
      fontSize={9}
      fontFamily="Arial"
      textAnchor="left"
      fill={black}
      dx={15}
      style={{ pointerEvents: "none" }}
    >
      {name}
    </text>
  );
}
/** Handles rendering Root, Parent, and other Nodes. */
function Node({ node, handleClick, isTargetNode }: NodeProps) {
  const isRoot = node.depth === 0;
  const isParent = !!node.children;

  if (isRoot)
    return (
      <RootNode
        node={node}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
      />
    );
  if (isParent)
    return (
      <ParentNode
        node={node}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
      />
    );

  return (
    <Group top={node.x} left={node.y}>
      <RectOrCircle
        node={node.data}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
      />
      <Text name={node.data.name} />
    </Group>
  );
}

const defaultMargin = { top: 10, left: 0, right: 0, bottom: 10 };

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
  const nodeSize = [15, 1000 / (data?.height ?? 1)]; // Increase width and height for more spacing

  const initialTransform = {
    scaleX: 1,
    scaleY: 1,
    translateX: 10,
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
                <rect width={width} height={height} rx={14} fill={white} />
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
                            stroke={secondaryColor}
                            strokeWidth="1"
                            fill="none"
                          />
                        ))}
                        {tree.descendants().map((node, i) => (
                          <Node
                            key={`node-${i}`}
                            node={node}
                            isTargetNode={
                              node.data.id.replace("_", ":") === cellTypeId
                            }
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
