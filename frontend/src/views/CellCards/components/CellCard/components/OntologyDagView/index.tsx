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
import { Divider } from "../../style";
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
const white = "#ffffff";
const black = "#000000";

type HierarchyNode = HierarchyPointNode<TreeNode>;

interface NodeProps {
  node: HierarchyNode;
  handleClick: MouseEventHandler<SVGGElement>;
  isTargetNode: boolean;
  handleMouseOver: (
    event: React.MouseEvent<SVGElement>,
    datum: TreeNode
  ) => void;
  handleMouseOut: () => void;
}

interface TreeNodeWithCoords extends TreeNode {
  x: number;
  y: number;
}

function RootNode({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
}: NodeProps) {
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
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
      <Text name={node.data.name} />
    </Group>
  );
}

function ParentNode({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
}: NodeProps) {
  return (
    <Group top={node.x} left={node.y}>
      <RectOrCircle
        node={node.data}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
      <Text name={node.data.name} />
    </Group>
  );
}
interface RectOrCircleProps {
  handleClick: MouseEventHandler<SVGGElement>;
  isTargetNode: boolean;
  node: TreeNode;
  handleMouseOver: (
    event: React.MouseEvent<SVGElement>,
    datum: TreeNode
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
  if (node.id === "") {
    color = white;
    stroke = black;
    //strokeDasharray = "1.5";
  }

  const cursor = node.id !== "" ? "pointer" : "default";
  const clickHandler = node.id !== "" ? handleClick : undefined;

  const onMouseOver =
    node.id === ""
      ? undefined
      : (event: React.MouseEvent<SVGElement>) => {
          handleMouseOver(event, node);
        };
  const onMouseOut = node.id === "" ? undefined : handleMouseOut;
  return node.hasChildren ? (
    <circle
      r={size}
      fill={color}
      onClick={clickHandler}
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
      onClick={clickHandler}
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
      style={{ pointerEvents: "none" }}
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
}: NodeProps) {
  const isRoot = node.depth === 0;
  const isParent = !!node.children;

  if (isRoot)
    return (
      <RootNode
        node={node}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
    );
  if (isParent)
    return (
      <ParentNode
        node={node}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
    );

  return (
    <Group top={node.x} left={node.y}>
      <Text name={node.data.name} />
      <RectOrCircle
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

  const { data: rawTree } = useCellOntologyTree(cellTypeId);
  const [initialTransformMatrix, setInitialTransformMatrix] = useState<
    typeof initialTransformMatrixDefault
  >(initialTransformMatrixDefault);
  const [centerNodeCoords, setCenterNodeCoords] = useState<
    [number, number] | null
  >(null);

  const router = useRouter();
  const data = useMemo(() => {
    if (!rawTree) return null;
    return hierarchy(rawTree);
  }, [rawTree]);

  const yMax = height - margin.top - margin.bottom;
  const xMax = width - margin.left - margin.right;

  // Customize nodeSize and separation
  const nodeSize = [15, 1000 / (data?.height ?? 1)]; // Increase width and height for more spacing

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
        | TreeNodeWithCoords
        | undefined;
      if (!targetNode) {
        targetNode = data
          .descendants()
          .find((node) => node.data.id === cellTypeId) as
          | TreeNodeWithCoords
          | undefined;
      }
      if (targetNode) {
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
    datum: TreeNode
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
      <Divider />
      {data ? (
        <Zoom<SVGSVGElement>
          key={centerNodeCoords ? "centered" : "initial"}
          width={width}
          height={height}
          scaleXMin={0.1}
          scaleXMax={4}
          scaleYMin={0.1}
          scaleYMax={4}
          initialTransformMatrix={initialTransformMatrix}
        >
          {(zoom) => (
            <div
              ref={containerRef}
              className="relative"
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
                            isTargetNode={node.data.id === cellTypeId}
                            handleMouseOver={handleMouseOver}
                            handleMouseOut={hideTooltip}
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
                <g>
                  <rect
                    x={width - 250}
                    y={0}
                    width={240}
                    height={50}
                    fill="white"
                    stroke="black"
                    strokeWidth={1}
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
