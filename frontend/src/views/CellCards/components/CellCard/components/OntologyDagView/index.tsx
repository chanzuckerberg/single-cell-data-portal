import React, {
  useMemo,
  MouseEventHandler,
  useEffect,
  useState,
  useLayoutEffect,
} from "react";
import { Group } from "@visx/group";
import { localPoint } from "@visx/event";
import { useTooltip, useTooltipInPortal } from "@visx/tooltip";
import { Tree, hierarchy } from "@visx/hierarchy";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { LinkHorizontal } from "@visx/shape";
import NodeGroup from "react-move/NodeGroup";

import { TableTitleWrapper, TableTitle } from "../common/style";
import { ONTOLOGY_SECTION_ID } from "../CellCardSidebar";
import { Zoom } from "@visx/zoom";
import { RectClipPath } from "@visx/clip-path";
import {
  CellOntologyTree as TreeNode,
  useCellOntologyTree,
  useCellOntologyTreeState,
} from "src/common/queries/cellCards";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";
import {
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import { StyledLegendText } from "./style";

const primaryColor = "#0073FF";
const secondaryColor = "#999999";
const tertiaryColor = "#CCCCCC";
const highlightColor = "#3CB371";
const white = "#F8F8F8";
const black = "#000000";

interface TreeNodeWithState extends TreeNode {
  isExpanded?: boolean;
  showAllChildren?: boolean;
  x0?: number;
  y0?: number;
  x?: number;
  y?: number;
}

type HierarchyNode = HierarchyPointNode<TreeNodeWithState>;

interface NodeProps {
  node: HierarchyNode;
  handleClick?: MouseEventHandler<SVGGElement>;
  isTargetNode: boolean;
  handleMouseOver: (
    event: React.MouseEvent<SVGElement>,
    datum: TreeNodeWithState
  ) => void;
  handleMouseOut: () => void;
  left: number;
  top: number;
  opacity: number;
  key: string;
}

function RootNode({
  node,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  left,
  top,
  key,
  opacity,
}: NodeProps) {
  const router = useRouter();
  return (
    <Group
      top={top}
      left={left}
      style={{ cursor: "pointer" }}
      key={key}
      opacity={opacity}
    >
      <RectOrCircle
        key={key}
        node={node.data}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
      <g
        onClick={() => {
          router.push(`${ROUTES.CELL_CARDS}/${node.data.id.replace(":", "_")}`);
        }}
      >
        <a>
          <Text name={node.data.name} />
        </a>
      </g>
    </Group>
  );
}

function ParentNode({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  left,
  top,
  key,
  opacity,
}: NodeProps) {
  const router = useRouter();
  return (
    <Group top={top} left={left} key={key} opacity={opacity}>
      <RectOrCircle
        node={node.data}
        key={key}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
      <g
        onClick={() => {
          router.push(`${ROUTES.CELL_CARDS}/${node.data.id.replace(":", "_")}`);
        }}
      >
        <a>
          <Text name={node.data.name} />
        </a>
      </g>
    </Group>
  );
}
interface RectOrCircleProps {
  handleClick?: MouseEventHandler<SVGGElement>;
  isTargetNode: boolean;
  node: TreeNodeWithState;
  handleMouseOver: (
    event: React.MouseEvent<SVGElement>,
    datum: TreeNodeWithState
  ) => void;
  handleMouseOut: () => void;
}

const smallSize = 4;
const largeSize = 8;

function RectOrCircle({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
}: RectOrCircleProps) {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
  if (isTargetNode) {
    color = highlightColor;
  }
  const size = node.n_cells === 0 ? smallSize : largeSize;
  let stroke = "none";
  let strokeDasharray;
  if (node.name === "") {
    color = white;
    stroke = black;
    //strokeDasharray = "1.5";
  }

  const cursor = "pointer";

  const onMouseOver =
    node.name === ""
      ? undefined
      : (event: React.MouseEvent<SVGElement>) => {
          handleMouseOver(event, node);
        };
  const onMouseOut = node.name === "" ? undefined : handleMouseOut;
  return node?.children?.length ? (
    <circle
      r={size}
      fill={color}
      onClick={handleClick}
      style={{ cursor: cursor }}
      stroke={stroke}
      strokeWidth={0.5}
      strokeDasharray={strokeDasharray}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
    />
  ) : (
    <rect
      height={size * 2}
      width={size * 2}
      y={-size}
      x={-size}
      fill={color}
      onClick={handleClick}
      style={{ cursor: cursor }}
      stroke={stroke}
      strokeWidth={0.5}
      strokeDasharray={strokeDasharray}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
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
    >
      {name}
    </text>
  );
}
/** Handles rendering Root, Parent, and other Nodes. */
function Node({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  left,
  top,
  key,
  opacity,
}: NodeProps) {
  const isRoot = node.depth === 0;
  const isParent = !!node.children;

  if (isRoot)
    return (
      <RootNode
        node={node}
        left={left}
        top={top}
        key={key}
        opacity={opacity}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
    );
  if (isParent)
    return (
      <ParentNode
        node={node}
        left={left}
        top={top}
        opacity={opacity}
        key={key}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
    );

  const router = useRouter();

  return (
    <Group top={top} left={left} key={key} opacity={opacity}>
      <g
        onClick={() => {
          router.push(`${ROUTES.CELL_CARDS}/${node.data.id.replace(":", "_")}`);
        }}
      >
        <a>
          <Text name={node.data.name} />
        </a>
      </g>
      <RectOrCircle
        key={key}
        node={node.data}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
    </Group>
  );
}

const defaultMargin = { top: 50, left: 0, right: 0, bottom: 0 };

interface LegendProps {
  xPos: number;
  yPos: number;
}
const InCorpusLegend = ({ xPos, yPos }: LegendProps) => {
  return (
    <g>
      <StyledLegendText x={xPos + 2.5} y={yPos + 30}>
        Yes
      </StyledLegendText>
      <circle
        cx={xPos + 12.5}
        cy={yPos + 15}
        fill={primaryColor}
        r={largeSize}
      />
      <StyledLegendText x={xPos} y={yPos}>
        In Corpus
      </StyledLegendText>
      <StyledLegendText x={xPos + 33} y={yPos + 30}>
        No
      </StyledLegendText>
      <circle
        cx={xPos + 40}
        cy={yPos + 15}
        fill={tertiaryColor}
        r={largeSize}
      />
    </g>
  );
};

const DescendantsLegend = ({ xPos, yPos }: LegendProps) => {
  return (
    <g>
      <defs>
        <pattern
          id="crosshatch"
          width="5"
          height="5"
          patternUnits="userSpaceOnUse"
        >
          <path d="M 5 0 L 0 5" stroke="#999999" strokeWidth="0.5" />
        </pattern>
      </defs>
      <StyledLegendText x={xPos + 2.5} y={yPos + 30}>
        Yes
      </StyledLegendText>
      <circle
        cx={xPos + 12.5}
        cy={yPos + 15}
        fill="url(#crosshatch)"
        stroke="#999999"
        strokeWidth={1}
        r={largeSize}
      />
      <StyledLegendText x={xPos - 5} y={yPos}>
        Descendants
      </StyledLegendText>
      <StyledLegendText x={xPos + 33.5} y={yPos + 30}>
        No
      </StyledLegendText>
      <rect
        x={xPos + 40 - largeSize + 2}
        y={yPos + 15 - largeSize + 1}
        fill="url(#crosshatch)"
        stroke="#999999"
        strokeWidth={1}
        width={largeSize * 2 - 2}
        height={largeSize * 2 - 2}
      />
    </g>
  );
};

const CollapsedNodesLegend = ({ xPos, yPos }: LegendProps) => {
  return (
    <g>
      <rect
        x={xPos + 30 - largeSize + 2}
        y={yPos + 15 - largeSize + 1}
        fill={white}
        stroke="#999999"
        strokeWidth={1}
        width={largeSize * 2 - 2}
        height={largeSize * 2 - 2}
      />
      <StyledLegendText x={xPos - 5} y={yPos}>
        Hidden terms
      </StyledLegendText>
    </g>
  );
};

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
  const initialTransformMatrixDefault = {
    scaleX: 1,
    scaleY: 1,
    translateX: 10,
    translateY: height / 2,
    skewX: 0,
    skewY: 0,
  };
  const [triggerRender, setTriggerRender] = useState(false);
  const { data: rawTree } = useCellOntologyTree();
  const { data: initialTreeState } = useCellOntologyTreeState(cellTypeId);
  const expandedNodes = initialTreeState?.isExpandedNodes ?? [];
  const notShownWhenExpandedNodes =
    initialTreeState?.notShownWhenExpandedNodes ?? {};

  const treeData: TreeNodeWithState | null = useMemo(() => {
    const newTree = rawTree ? JSON.parse(JSON.stringify(rawTree)) : null;
    if (newTree) {
      setNodesStateInRawTree(newTree, expandedNodes);
    }
    return newTree;
  }, [rawTree, expandedNodes, notShownWhenExpandedNodes]);

  const [initialTransformMatrix, setInitialTransformMatrix] = useState<
    typeof initialTransformMatrixDefault
  >(initialTransformMatrixDefault);
  const [centerNodeCoords, setCenterNodeCoords] = useState<
    [number, number] | null
  >(null);

  const data = useMemo(() => {
    if (!treeData) return null;
    return hierarchy(treeData, (d) => {
      const isExpanded = d.isExpanded;

      const newChildren: TreeNodeWithState[] = [];
      if (isExpanded) {
        const appendDummy = d.id in notShownWhenExpandedNodes;
        for (const child of d?.children ?? []) {
          if (
            d.showAllChildren ||
            !appendDummy ||
            (appendDummy && !notShownWhenExpandedNodes[d.id].includes(child.id))
          ) {
            newChildren.push(child);
          }
        }
        if (appendDummy && !d.showAllChildren) {
          newChildren.push({
            id: `${d.id}-dummy-child`,
            name: "",
            n_cells: 0,
            n_cells_rollup: 0,
            isExpanded: false,
          });
        }
      }
      return newChildren.length ? newChildren : null;
    });
  }, [treeData, triggerRender]);

  const yMax = height - margin.top - margin.bottom;
  const xMax = width - margin.left - margin.right;

  // Customize nodeSize and separation
  const nodeSize = [15, 200]; // Increase width and height for more spacing

  useEffect(() => {
    setCenterNodeCoords(null);
    return () => {
      setCenterNodeCoords(null);
      setInitialTransformMatrix(initialTransformMatrixDefault);
      hideTooltip();
    };
  }, [cellTypeId]);

  // This effect is used to center the node corresponding to the cell type id
  // useLayoutEffect is used to ensure that node coordinates have been populated by the renderer prior to painting
  const useIsomorphicLayoutEffect =
    typeof window !== "undefined" ? useLayoutEffect : useEffect;
  useIsomorphicLayoutEffect(() => {
    if (data) {
      let targetNode = data
        .descendants()
        .find((node) => node.data.id === cellTypeId && node.data.children) as
        | TreeNodeWithState
        | undefined;
      if (!targetNode) {
        targetNode = data
          .descendants()
          .find((node) => node.data.id === cellTypeId) as
          | TreeNodeWithState
          | undefined;
      }
      if (targetNode) {
        if (targetNode.x && targetNode.y) {
          setCenterNodeCoords([targetNode.x, targetNode.y]);
          setInitialTransformMatrix({
            scaleX: 1,
            scaleY: 1,
            translateX: width / 2 - targetNode.y,
            translateY: height / 2 - targetNode.x,
            skewX: 0,
            skewY: 0,
          });
        }
      }
    }
  }, [data, cellTypeId, width, height]);

  const {
    tooltipData,
    tooltipLeft,
    tooltipTop,
    tooltipOpen,
    showTooltip,
    hideTooltip,
  } = useTooltip<{ n_cells: number; n_cells_rollup: number }>();

  const { TooltipInPortal, containerRef } = useTooltipInPortal({
    detectBounds: true,
    scroll: true,
  });
  const handleMouseOver = (
    event: React.MouseEvent<SVGElement>,
    datum: TreeNodeWithState
  ) => {
    if (event.target instanceof SVGElement) {
      if (event.target.ownerSVGElement !== null) {
        const coords = localPoint(event.target.ownerSVGElement, event);
        if (coords) {
          showTooltip({
            tooltipLeft: coords.x,
            tooltipTop: coords.y,
            tooltipData: {
              n_cells: datum.n_cells,
              n_cells_rollup: datum.n_cells_rollup,
            },
          });
        }
      }
    }
  };

  return (
    <>
      <TableTitleWrapper id={ONTOLOGY_SECTION_ID}>
        <TableTitle>Cell Ontology</TableTitle>
      </TableTitleWrapper>
      {data ? (
        <Zoom<SVGSVGElement>
          key={centerNodeCoords ? "centered" : "initial"}
          width={width}
          height={height}
          scaleXMin={0.25}
          scaleXMax={4}
          scaleYMin={0.25}
          scaleYMax={4}
          initialTransformMatrix={initialTransformMatrix}
          wheelDelta={(event: WheelEvent | React.WheelEvent<Element>) => {
            const newScale = event.deltaY > 0 ? 0.95 : 1.05;
            return { scaleX: newScale, scaleY: newScale };
          }}
        >
          {(zoom) => (
            <div
              ref={containerRef}
              onMouseDown={() => {
                hideTooltip();
              }}
            >
              {tooltipOpen && (
                <TooltipInPortal
                  // set this to random so it correctly updates with parent bounds
                  key={Math.random()}
                  top={tooltipTop}
                  left={tooltipLeft}
                >
                  <div>
                    <b>{tooltipData?.n_cells}</b> cells
                    <br />
                    <b>{tooltipData?.n_cells_rollup}</b> descendant cells
                  </div>
                </TooltipInPortal>
              )}
              <svg
                width={width}
                height={height}
                ref={zoom.containerRef}
                style={{ cursor: zoom.isDragging ? "grabbing" : "grab" }}
              >
                <RectClipPath id="zoom-clip" width={width} height={height} />
                <rect width={width} height={height} rx={14} fill={white} />
                <g transform={zoom.toString()}>
                  <Tree<TreeNodeWithState>
                    root={data}
                    size={[yMax, xMax]}
                    nodeSize={nodeSize as [number, number]}
                  >
                    {(tree) => (
                      <Group top={margin.top} left={margin.left}>
                        <NodeGroup
                          data={tree.links()}
                          keyAccessor={(d) =>
                            `${d.source.data.id}_${d.target.data.id}`
                          }
                          start={({ source }) => {
                            return {
                              source: {
                                x: source.data.x0 ?? source.x,
                                y: source.data.y0 ?? source.y,
                              },
                              target: {
                                x: source.data.x0 ?? source.x,
                                y: source.data.y0 ?? source.y,
                              },
                            };
                          }}
                          enter={({ source, target }) => {
                            return {
                              source: {
                                x: [source.x],
                                y: [source.y],
                              },
                              target: {
                                x: [target.x],
                                y: [target.y],
                              },
                            };
                          }}
                          update={({ source, target }) => {
                            return {
                              source: {
                                x: [source.x],
                                y: [source.y],
                              },
                              target: {
                                x: [target.x],
                                y: [target.y],
                              },
                            };
                          }}
                          leave={({ source, target }) => {
                            const collapsedParent = findCollapsedParent(source);
                            return collapsedParent
                              ? {
                                  source: {
                                    x: [collapsedParent.data.x0],
                                    y: [collapsedParent.data.y0],
                                  },
                                  target: {
                                    x: [collapsedParent.data.x0],
                                    y: [collapsedParent.data.y0],
                                  },
                                }
                              : {
                                  source: {
                                    x: [source.x],
                                    y: [source.y],
                                  },
                                  target: {
                                    x: [target.x],
                                    y: [target.y],
                                  },
                                };
                          }}
                        >
                          {(nodes) => {
                            return (
                              <Group>
                                {nodes.map(({ key, state }) => {
                                  return (
                                    <LinkHorizontal
                                      key={key}
                                      data={state}
                                      stroke={secondaryColor}
                                      strokeWidth="1"
                                      fill="none"
                                    />
                                  );
                                })}
                              </Group>
                            );
                          }}
                        </NodeGroup>
                        <NodeGroup
                          data={tree.descendants()}
                          keyAccessor={(
                            d: HierarchyPointNode<TreeNodeWithState>
                          ) => d.data.id}
                          start={(
                            node: HierarchyPointNode<TreeNodeWithState>
                          ) => {
                            return node.parent
                              ? {
                                  top: node.parent.x,
                                  left: node.parent.y,
                                  opacity: 0,
                                }
                              : {
                                  top: 0,
                                  left: 0,
                                  opacity: 0,
                                };
                          }}
                          enter={(
                            node: HierarchyPointNode<TreeNodeWithState>
                          ) => {
                            return {
                              top: [node.x],
                              left: [node.y],
                              opacity: [1],
                            };
                          }}
                          update={(
                            node: HierarchyPointNode<TreeNodeWithState>
                          ) => {
                            return {
                              top: [node.x],
                              left: [node.y],
                              opacity: [1],
                            };
                          }}
                          leave={(
                            node: HierarchyPointNode<TreeNodeWithState>
                          ) => {
                            const collapsedParent = node.parent
                              ? findCollapsedParent(node.parent)
                              : null;
                            const collapsedParentPrevPos = collapsedParent
                              ? {
                                  x: collapsedParent.data.x0,
                                  y: collapsedParent.data.y0,
                                }
                              : {
                                  x: 0,
                                  y: 0,
                                };
                            return {
                              top: [collapsedParentPrevPos.x],
                              left: [collapsedParentPrevPos.y],
                              opacity: [0],
                            };
                          }}
                        >
                          {(nodes) => {
                            return (
                              <Group>
                                {nodes.map(({ key, data: node, state }) => (
                                  <Node
                                    key={key}
                                    node={node}
                                    isTargetNode={
                                      node.data.id.split("__").at(0) ===
                                      cellTypeId
                                    }
                                    handleMouseOver={handleMouseOver}
                                    handleMouseOut={hideTooltip}
                                    left={state.left}
                                    top={state.top}
                                    opacity={state.opacity}
                                    handleClick={() => {
                                      if (!node.data.isExpanded) {
                                        node.data.x0 = node.x;
                                        node.data.y0 = node.y;
                                      }
                                      if (node.data.name === "") {
                                        if (node.parent) {
                                          node.parent.data.showAllChildren =
                                            true;
                                        }
                                      }
                                      node.data.isExpanded =
                                        !node.data.isExpanded;
                                      if (!node.data.isExpanded) {
                                        node.data.showAllChildren = false;
                                      }
                                      setTriggerRender(!triggerRender);
                                    }}
                                  />
                                ))}
                              </Group>
                            );
                          }}
                        </NodeGroup>
                      </Group>
                    )}
                  </Tree>
                </g>
                <g>
                  <rect
                    x={width - 260}
                    y={-10}
                    width={240}
                    height={60}
                    fill={white}
                  />
                  <path
                    d={`M ${width - 260} ${-10} L ${width - 260} ${45} Q ${
                      width - 260
                    } ${50} ${width - 255} ${50} L ${width} ${50}`}
                    stroke="black"
                    fill="transparent"
                  />

                  <InCorpusLegend xPos={width - 240} yPos={10} />
                  <DescendantsLegend xPos={width - 160} yPos={10} />
                  <CollapsedNodesLegend xPos={width - 80} yPos={10} />
                </g>
              </svg>
            </div>
          )}
        </Zoom>
      ) : (
        <TableUnavailableContainer>
          <TableUnavailableHeader>
            Cell ontology visualization unavailable
          </TableUnavailableHeader>
          <TableUnavailableDescription>
            The cell ontology visualization for this cell type is currently not
            available. Please check again later.
          </TableUnavailableDescription>
        </TableUnavailableContainer>
      )}
    </>
  );
}

function setNodesStateInRawTree(
  graph: TreeNodeWithState,
  isExpandedNodes: string[]
) {
  graph.isExpanded = isExpandedNodes.includes(graph.id);
  if (graph.children) {
    for (const child of graph.children) {
      setNodesStateInRawTree(child, isExpandedNodes);
    }
  }
}

export function findCollapsedParent(
  node: HierarchyPointNode<TreeNodeWithState>
): HierarchyPointNode<TreeNodeWithState> | null {
  if (!node.data.isExpanded) {
    return node;
  } else if (node.parent) {
    return findCollapsedParent(node.parent);
  } else {
    return null;
  }
}
