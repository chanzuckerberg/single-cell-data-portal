import React, {
  useMemo,
  MouseEventHandler,
  useEffect,
  useState,
  useRef,
} from "react";
import { Group } from "@visx/group";
import { Global } from "@emotion/react";
import { localPoint } from "@visx/event";
import { useTooltip, useTooltipInPortal } from "@visx/tooltip";
import { Tree, hierarchy } from "@visx/hierarchy";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { LinkHorizontal } from "@visx/shape";
import NodeGroup from "react-move/NodeGroup";
import FullscreenIcon from "@material-ui/icons/Fullscreen";
import FullscreenExitIcon from "@material-ui/icons/FullscreenExit";
import { TableTitleWrapper, TableTitle } from "../common/style";
import { ONTOLOGY_SECTION_ID } from "../CellCardSidebar";
import { Zoom } from "@visx/zoom";
import { RectClipPath } from "@visx/clip-path";
import {
  CellOntologyTree as TreeNode,
  useCellOntologyTree,
  useCellOntologyTreeState,
  useCellTypesById,
} from "src/common/queries/cellCards";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";
import {
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../common/style";
import {
  FullscreenButton,
  HoverContainer,
  StyledLegendText,
  TooltipInPortalStyle,
} from "./style";
import { useFullScreen } from "../FullScreenProvider";

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
  animationKey: string;
  maxWidth: number;
  isInCorpus: boolean;
}

function RootNode({
  node,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  handleClick,
  left,
  top,
  animationKey,
  opacity,
  maxWidth,
  isInCorpus,
}: NodeProps) {
  const router = useRouter();
  return (
    <Group
      top={top}
      left={left}
      style={{ cursor: "pointer" }}
      key={animationKey}
      opacity={opacity}
    >
      <RectOrCircle
        animationKey={animationKey}
        node={node.data}
        isTargetNode={isTargetNode}
        handleClick={handleClick}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
      {isInCorpus ? (
        <g
          onClick={() => {
            router.push(
              `${ROUTES.CELL_CARDS}/${node.data.id
                .replace(":", "_")
                .split("__")
                .at(0)}`
            );
          }}
        >
          <a>
            <Text
              isInCorpus={isInCorpus}
              name={node.data.name}
              maxWidth={maxWidth}
            />
          </a>
        </g>
      ) : (
        <Text
          isInCorpus={isInCorpus}
          name={node.data.name}
          maxWidth={maxWidth}
        />
      )}
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
  animationKey,
  opacity,
  maxWidth,
  isInCorpus,
}: NodeProps) {
  const router = useRouter();
  return (
    <Group top={top} left={left} key={animationKey} opacity={opacity}>
      <RectOrCircle
        node={node.data}
        animationKey={animationKey}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
      {isInCorpus ? (
        <g
          onClick={() => {
            router.push(
              `${ROUTES.CELL_CARDS}/${node.data.id
                .replace(":", "_")
                .split("__")
                .at(0)}`
            );
          }}
        >
          <a>
            <Text
              isInCorpus={isInCorpus}
              name={node.data.name}
              maxWidth={maxWidth}
            />
          </a>
        </g>
      ) : (
        <Text
          isInCorpus={isInCorpus}
          name={node.data.name}
          maxWidth={maxWidth}
        />
      )}
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
  animationKey: string;
}

const smallSize = 4;
const largeSize = 8;

function RectOrCircle({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  animationKey,
}: RectOrCircleProps) {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
  if (isTargetNode) {
    color = highlightColor;
  }
  const size = node.n_cells === 0 ? smallSize : largeSize;
  if (node.id.startsWith("dummy-child")) {
    color = primaryColor;
  }

  const onMouseOver = node.id.startsWith("dummy-child")
    ? undefined
    : (event: React.MouseEvent<SVGElement>) => {
        handleMouseOver(event, node);
      };
  const onMouseOut = node.name.startsWith("dummy-child")
    ? undefined
    : handleMouseOut;
  return node?.children?.length || node.id.startsWith("dummy-child") ? (
    <circle
      r={size}
      fill={color}
      key={animationKey}
      onClick={handleClick}
      style={{ cursor: "pointer" }}
      strokeWidth={0.5}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
    />
  ) : (
    <rect
      height={size * 2}
      width={size * 2}
      y={-size}
      x={-size}
      key={animationKey}
      fill={color}
      onClick={handleClick}
      style={{ cursor: "pointer" }}
      strokeWidth={0.5}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
    />
  );
}

interface TextProps {
  name: string;
  maxWidth: number;
  isInCorpus: boolean;
}
function Text({ isInCorpus, name, maxWidth }: TextProps) {
  const textRef = useRef<SVGTextElement>(null);

  useEffect(() => {
    const textElement = textRef.current;
    let textContent = name;
    if (textElement) {
      while (textElement.getComputedTextLength() > maxWidth) {
        textContent = textContent.slice(0, -1);
        textElement.textContent = textContent + "...";
      }
    }
  }, [name]);

  return (
    <text
      ref={textRef}
      dy=".33em"
      fontSize={9}
      fontFamily="Arial"
      textAnchor="left"
      fontStyle={
        isInCorpus || name.endsWith("cell types") ? "normal" : "italic"
      }
      fontWeight={isInCorpus ? "bold" : "normal"}
      fill={black}
      dx={10}
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
  animationKey,
  opacity,
  maxWidth,
  isInCorpus,
}: NodeProps) {
  const isRoot = node.depth === 0;
  const isParent = !!node.children;

  if (isRoot)
    return (
      <RootNode
        node={node}
        left={left}
        top={top}
        animationKey={animationKey}
        opacity={opacity}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
        handleClick={handleClick}
        maxWidth={maxWidth}
        isInCorpus={isInCorpus}
      />
    );
  if (isParent)
    return (
      <ParentNode
        node={node}
        left={left}
        top={top}
        opacity={opacity}
        animationKey={animationKey}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
        maxWidth={maxWidth}
        isInCorpus={isInCorpus}
      />
    );

  const router = useRouter();

  return (
    <Group top={top} left={left} key={animationKey} opacity={opacity}>
      {isInCorpus ? (
        <g
          onClick={() => {
            router.push(
              `${ROUTES.CELL_CARDS}/${node.data.id
                .replace(":", "_")
                .split("__")
                .at(0)}`
            );
          }}
        >
          <a>
            <Text
              isInCorpus={isInCorpus}
              name={node.data.name}
              maxWidth={maxWidth}
            />
          </a>
        </g>
      ) : (
        <Text
          isInCorpus={isInCorpus}
          name={node.data.name}
          maxWidth={maxWidth}
        />
      )}
      <RectOrCircle
        animationKey={`${animationKey}-rect-or-circle`}
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
  skinnyMode: boolean;
};

export default function OntologyDagView({ skinnyMode, cellTypeId }: TreeProps) {
  const margin = defaultMargin;
  const [width, setWidth] = useState(1000);
  const [height, setHeight] = useState(500);

  const [resizeWidth, setResizeWidth] = useState(1000);

  const {
    isFullScreen,
    screenDimensions,
    enableFullScreen,
    disableFullScreen,
  } = useFullScreen();

  useEffect(() => {
    const skinnyAdjustment = skinnyMode ? 0 : 120 + 240;
    const width = Math.min(
      1000,
      window.innerWidth - 40 - 40 - skinnyAdjustment
    );
    setResizeWidth(width);
    if (!isFullScreen) setWidth(width);

    const handleResize = () => {
      const skinnyAdjustment = skinnyMode ? 0 : 120 + 240;
      const width = Math.min(
        1000,
        window.innerWidth - 40 - 40 - skinnyAdjustment
      );
      setResizeWidth(width);
      if (!isFullScreen) setWidth(width);
    };
    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  }, [isFullScreen, skinnyMode]);

  useEffect(() => {
    let newWidth = resizeWidth;
    let newHeight = 500;
    if (screenDimensions.width > 0 && isFullScreen) {
      newWidth = screenDimensions.width;
    }
    if (screenDimensions.height > 0 && isFullScreen) {
      newHeight = screenDimensions.height;
    }
    setWidth(newWidth);
    setHeight(newHeight);
  }, [screenDimensions, isFullScreen, resizeWidth]);

  const initialTransformMatrixDefault = {
    scaleX: 1,
    scaleY: 1,
    translateX: 10,
    translateY: height / 2,
    skewX: 0,
    skewY: 0,
  };
  const [triggerRender, setTriggerRender] = useState(false);
  const [duration, setDuration] = useState(0);
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
          const numHiddenChildren =
            (d.children?.length ?? 0) - newChildren.length;
          newChildren.push({
            id: `dummy-child-${d.id}`,
            name: `${numHiddenChildren} cell types`,
            n_cells: 0,
            n_cells_rollup: 0,
            isExpanded: false,
          });
        } else if (d.showAllChildren) {
          if (d.children) {
            const firstId = d.children[0].id;
            newChildren.sort((a, b) => {
              if (a.id === firstId) return -1;
              if (b.id === firstId) return 1;
              return a.id.localeCompare(b.id);
            });
          }
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

  useEffect(() => {
    if (data) {
      data.each((node) => {
        const pointNode = node as HierarchyPointNode<TreeNodeWithState>;
        if (pointNode.x && pointNode.y) {
          node.data.x0 = pointNode.x;
          node.data.y0 = pointNode.y;
        }
      });
      let targetNode = data
        .descendants()
        .find(
          (node) =>
            node.data.id.split("__").at(0) === cellTypeId && node.data.children
        ) as HierarchyPointNode<TreeNodeWithState> | undefined;
      if (!targetNode) {
        targetNode = data
          .descendants()
          .find((node) => node.data.id.split("__").at(0) === cellTypeId) as
          | HierarchyPointNode<TreeNodeWithState>
          | undefined;
      }
      if (targetNode) {
        if (targetNode.x !== undefined && targetNode.y !== undefined) {
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

  useEffect(() => {
    setTriggerRender(!triggerRender);
  }, [isFullScreen]);

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
  const cellTypesById = useCellTypesById() || {};

  return (
    <>
      <Global styles={TooltipInPortalStyle} />
      <TableTitleWrapper id={ONTOLOGY_SECTION_ID}>
        <TableTitle>Cell Ontology</TableTitle>
      </TableTitleWrapper>
      {data && initialTreeState ? (
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
            <HoverContainer
              height={height}
              width={width}
              ref={containerRef}
              onMouseDown={() => {
                hideTooltip();
              }}
              isFullScreen={isFullScreen}
            >
              <FullscreenButton
                onClick={isFullScreen ? disableFullScreen : enableFullScreen}
              >
                {isFullScreen ? <FullscreenExitIcon /> : <FullscreenIcon />}
              </FullscreenButton>
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
                style={{
                  cursor: zoom.isDragging ? "grabbing" : "grab",
                  position: "absolute",
                  top: 0,
                  left: 0,
                }}
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
                              opacity: 0,
                              timing: { duration },
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
                              opacity: [1],
                              timing: { duration },
                            };
                          }}
                          update={({ source, target }) => {
                            const collapsedParent = findCollapsedParent(source);
                            return collapsedParent
                              ? {
                                  source: {
                                    x: [collapsedParent.x],
                                    y: [collapsedParent.y],
                                  },
                                  target: {
                                    x: [collapsedParent.x],
                                    y: [collapsedParent.y],
                                  },
                                  opacity: [1],
                                  timing: { duration },
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
                                  opacity: [1],
                                  timing: { duration },
                                };
                          }}
                          leave={() => {
                            return { opacity: [0], timing: { duration: 0 } };
                          }}
                        >
                          {(nodes) => {
                            return (
                              <Group>
                                {nodes.map(({ key, state }) => {
                                  return (
                                    <Group key={key} opacity={state.opacity}>
                                      <LinkHorizontal
                                        key={key}
                                        data={state}
                                        stroke={secondaryColor}
                                        strokeWidth="1"
                                        fill="none"
                                      />
                                    </Group>
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
                                  top: node.parent.data.x0 ?? node.parent.x,
                                  left: node.parent.data.y0 ?? node.parent.y,
                                  opacity: 0,
                                  timing: { duration },
                                }
                              : {
                                  top: node.x,
                                  left: node.y,
                                  opacity: 0,
                                  timing: { duration },
                                };
                          }}
                          enter={(
                            node: HierarchyPointNode<TreeNodeWithState>
                          ) => {
                            return {
                              top: [node.x],
                              left: [node.y],
                              opacity: [1],
                              timing: { duration },
                            };
                          }}
                          update={(
                            node: HierarchyPointNode<TreeNodeWithState>
                          ) => {
                            return {
                              top: [node.x],
                              left: [node.y],
                              opacity: [1],
                              timing: { duration },
                            };
                          }}
                          leave={() => {
                            return { opacity: [0], timing: { duration: 0 } };
                          }}
                        >
                          {(nodes) => {
                            return (
                              <Group>
                                {nodes.map(({ key, data: node, state }) => {
                                  const isInCorpus =
                                    (node.data.id
                                      .replace(":", "_")
                                      .split("__")
                                      .at(0)
                                      ?.replace("_", ":") ?? "") in
                                    cellTypesById;
                                  return (
                                    <Node
                                      key={`${key}-node`}
                                      isInCorpus={isInCorpus}
                                      animationKey={key}
                                      node={node}
                                      isTargetNode={
                                        node.data.id.split("__").at(0) ===
                                        cellTypeId
                                      }
                                      handleMouseOver={handleMouseOver}
                                      handleMouseOut={hideTooltip}
                                      maxWidth={180}
                                      left={state.left}
                                      top={state.top}
                                      opacity={state.opacity}
                                      handleClick={() => {
                                        if (duration === 0) {
                                          setDuration(250);
                                        }
                                        if (
                                          node.data.id.startsWith("dummy-child")
                                        ) {
                                          if (node.parent) {
                                            node.parent.data.showAllChildren =
                                              true;
                                          }
                                        }
                                        node.data.isExpanded =
                                          !node.data.isExpanded;
                                        if (!node.data.isExpanded) {
                                          node.data.showAllChildren = false;
                                          collapseAllDescendants(node);
                                        }
                                        setTriggerRender(!triggerRender);
                                      }}
                                    />
                                  );
                                })}
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
                    width={260}
                    height={60}
                    fill={white}
                    rx={4}
                  />
                  {/* <path
                    d={`M ${width - 260} ${-10} L ${width - 260} ${45} Q ${
                      width - 260
                    } ${50} ${width - 255} ${50} L ${width} ${50}`}
                    stroke="black"
                    strokeWidth={0.25}
                    fill="transparent"
                  /> */}

                  <InCorpusLegend xPos={width - 240} yPos={10} />
                  <DescendantsLegend xPos={width - 160} yPos={10} />
                  <CollapsedNodesLegend xPos={width - 80} yPos={10} />
                </g>
              </svg>
            </HoverContainer>
          )}
        </Zoom>
      ) : (
        <TableUnavailableContainer>
          <TableUnavailableHeader>
            Cell ontology visualization unavailable
          </TableUnavailableHeader>
          <TableUnavailableDescription>
            The cell ontology visualization for this cell type is currently not
            available as it is not a descendant of "animal cell".
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

function collapseAllDescendants(node: HierarchyPointNode<TreeNodeWithState>) {
  if (node.children) {
    for (const child of node.children) {
      child.data.isExpanded = false;
      collapseAllDescendants(child);
    }
  }
}
