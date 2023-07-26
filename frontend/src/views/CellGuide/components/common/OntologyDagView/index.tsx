import React, { useMemo, useEffect, useState } from "react";
import { Group } from "@visx/group";
import { Global } from "@emotion/react";
import { useTooltip, useTooltipInPortal } from "@visx/tooltip";
import { Tree, hierarchy } from "@visx/hierarchy";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import FullscreenIcon from "@mui/icons-material/Fullscreen";
import FullscreenExitIcon from "@mui/icons-material/FullscreenExit";
import {
  TableTitleWrapper,
  TableTitle,
} from "../../CellGuideCard/components/common/style";
import { Zoom } from "@visx/zoom";
import { RectClipPath } from "@visx/clip-path";
import {
  TissueCountsPerCellType,
  useCellOntologyTree,
  useCellOntologyTreeState,
  useCellOntologyTreeStateTissue,
} from "src/common/queries/cellGuide";
import {
  TableUnavailableContainer,
  TableUnavailableHeader,
  TableUnavailableDescription,
} from "../../CellGuideCard/components/common/style";
import {
  FullscreenButton,
  HoverContainer,
  TooltipInPortalStyle,
  StyledSVG,
} from "./style";
import { useFullScreen } from "../FullScreenProvider";
import {
  defaultMargin,
  backgroundColor,
  DEFAULT_ONTOLOGY_HEIGHT,
  DEFAULT_ONTOLOGY_WIDTH,
  NODE_SPACINGS,
} from "./common/constants";
import { TreeNodeWithState } from "./common/types";
import Legend from "./components/Legend";
import AnimatedNodes from "./components/AnimatedNodes";
import AnimatedLinks from "./components/AnimatedLinks";
import {
  LEFT_RIGHT_PADDING_PX,
  LEFT_RIGHT_PADDING_PX_SKINNY_MODE,
  SIDEBAR_COLUMN_GAP_PX,
} from "../../CellGuideCard/style";

export const CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP =
  "cell-guide-card-ontology-dag-view-tooltip";
export const CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW =
  "cell-guide-card-ontology-dag-view";
export const CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON =
  "cell-guide-card-ontology-dag-view-fullscreen-button";
export const CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER =
  "cell-guide-card-ontology-dag-view-hover-container";

interface BaseTreeProps {
  skinnyMode?: boolean;
  initialWidth?: number;
  initialHeight?: number;
}

interface CellTypeIdProps extends BaseTreeProps {
  cellTypeId: string;
  tissueId?: never;
  tissueName?: never;
}

interface TissueIdProps extends BaseTreeProps {
  cellTypeId?: never;
  tissueId: string;
  tissueName: string;
}

type TreeProps = CellTypeIdProps | TissueIdProps;

// This determines the initial Zoom position and scale
const initialTransformMatrixDefault = {
  scaleX: 1,
  scaleY: 1,
  translateX: 10,
  translateY: 500 / 2,
  skewX: 0,
  skewY: 0,
};

export default function OntologyDagView({
  cellTypeId,
  tissueId,
  tissueName,
  skinnyMode,
  initialHeight,
  initialWidth,
}: TreeProps) {
  skinnyMode = cellTypeId ? skinnyMode : true;

  const defaultHeight = initialHeight ?? DEFAULT_ONTOLOGY_HEIGHT;
  const defaultWidth = initialWidth ?? DEFAULT_ONTOLOGY_WIDTH;

  const [width, setWidth] = useState(defaultWidth);
  const [height, setHeight] = useState(defaultHeight);

  const [initialTransformMatrix, setInitialTransformMatrix] = useState<
    typeof initialTransformMatrixDefault
  >(initialTransformMatrixDefault);

  // This toggler is used for centering the Zoom component on the target cell type.
  // It triggers a re-render of the Zoom component so that updates to the initialTransformMatrix
  // take effect.
  const [centeredNodeCoords, setCenteredNodeCoords] = useState<boolean>(false);

  // This is used to store the desired resized width of the ontology view
  // while full screen mode is active.
  const [resizeWidth, setResizeWidth] = useState(defaultWidth);

  const {
    isFullScreen,
    screenDimensions,
    enableFullScreen,
    disableFullScreen,
  } = useFullScreen();

  // Handle the resizing of the ontology view when the screen is resized
  useEffect(() => {
    const skinnyAdjustment = skinnyMode ? 0 : SIDEBAR_COLUMN_GAP_PX + 240;

    const leftRightPadding = skinnyMode
      ? LEFT_RIGHT_PADDING_PX_SKINNY_MODE
      : LEFT_RIGHT_PADDING_PX;

    const width = Math.min(
      defaultWidth,
      window.innerWidth - leftRightPadding * 2 - skinnyAdjustment
    );
    setResizeWidth(width);
    if (!isFullScreen) setWidth(width);

    const handleResize = () => {
      const skinnyAdjustment = skinnyMode ? 0 : SIDEBAR_COLUMN_GAP_PX + 240;

      // Account for the padding on the left and right of the CellGuideCard component
      const width = Math.min(
        defaultWidth,
        window.innerWidth - leftRightPadding * 2 - skinnyAdjustment
      );
      // Always set the resize width, but only set the width if not in full screen mode
      setResizeWidth(width);
      if (!isFullScreen) setWidth(width);
    };
    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  }, [isFullScreen, skinnyMode, defaultWidth]);

  // Handle the resizing of the ontology view when full screen mode is toggled
  useEffect(() => {
    let newWidth = resizeWidth;
    let newHeight = defaultHeight;
    if (screenDimensions.width > 0 && isFullScreen) {
      newWidth = screenDimensions.width;
    }
    if (screenDimensions.height > 0 && isFullScreen) {
      newHeight = screenDimensions.height;
    }
    setWidth(newWidth);
    setHeight(newHeight);
  }, [screenDimensions, isFullScreen, resizeWidth, defaultHeight]);

  // This is used to trigger a re-render of the ontology view
  const [triggerRender, setTriggerRender] = useState(false);
  const toggleTriggerRender = () => {
    setTriggerRender(!triggerRender);
  };

  // Animation duration - initially zero so the animation doesn't play on load
  const [duration, setDuration] = useState(0);

  // The raw cell ontology tree data. This is called "rawTree" because it does not contain
  // the "isExpanded" property that is used to track the expanded state of each node, along with
  // other properties like the positions of the nodes.
  const { data: rawTree } = useCellOntologyTree();

  const entityId = cellTypeId ?? tissueId ?? "";
  const useTreeState = cellTypeId
    ? useCellOntologyTreeState
    : useCellOntologyTreeStateTissue;
  // Contains the nodes that are initially expanded (expandedNodes) and the nodes that are collapsed
  // by default when their parents are expanded (notShownWhenExpanded).
  // Also optionally contains cell counts for each node in the tissue ontology view (tissueCounts).
  const { data: initialTreeState } = useTreeState(entityId);

  // Populate the tree data structure nodes with "isExpanded".
  const treeData: TreeNodeWithState | null = useMemo(() => {
    const newTree = rawTree ? JSON.parse(JSON.stringify(rawTree)) : null;
    if (newTree) {
      const expandedNodes = initialTreeState?.isExpandedNodes ?? [];
      setIsExpandedInRawTree(newTree, expandedNodes);

      const tissueCounts = initialTreeState?.tissueCounts ?? {};
      if (Object.keys(tissueCounts).length > 0) {
        setCellCountsInRawTree(newTree, tissueCounts);
      }
    }
    return newTree;
  }, [rawTree, initialTreeState]);

  // Create the tree data structure
  const data = useMemo(() => {
    if (!treeData) return null;
    return hierarchy(treeData, (d) => {
      const newChildren: TreeNodeWithState[] = [];
      if (d.isExpanded) {
        const notShownWhenExpandedNodes =
          initialTreeState?.notShownWhenExpandedNodes ?? {};
        // If this node is a key in `notShownWhenExpandedNodes`, then it is a parent to nodes
        // that should be collapsed. Therefore, set `appendDummy` to `true`.
        // The term "dummy" here is used to describe the collapsed node containing hidden terms.
        // The text label of a "dummy" is the number of hidden terms, e.g. 52 cell types.
        const appendDummy = d.id in notShownWhenExpandedNodes;

        // `showAllChildren` is a flag that is set to `true` when the user clicks on a collapsed node.
        // It indicates that all of the children of the node should be shown.
        for (const child of d?.children ?? []) {
          if (
            d.showAllChildren ||
            !appendDummy ||
            // This is a child that should be shown because it is in a path to the target node.
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
          // Sorting such that the first child in the `newChildren` array is the
          // same as the first child in the original array (`d.children`) seems to
          // ensure that the position of the visible node is preserved when its dummy
          // sibling node is expanded. It makes for smoother animations.
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
  }, [treeData, triggerRender, initialTreeState]);

  // Cleanup when the cell type id changes
  useEffect(() => {
    setCenteredNodeCoords(false);
    return () => {
      setCenteredNodeCoords(false);
      setInitialTransformMatrix(initialTransformMatrixDefault);
      disableFullScreen();
      hideTooltip();
    };
  }, [cellTypeId, tissueId]);

  // This useEffect is used to set the initial transform matrix when the tree data changes.
  // The initial transform matrix is used to center the tree view on the target node.
  // This hook is also used to set the initial positions of the nodes in the tree to `x0` and `y0`.
  useEffect(() => {
    if (data) {
      // This block sets the initial positions of the nodes in the tree to `x0` and `y0`.
      // This is used when animating transitions in node positions.
      data.each((node) => {
        const pointNode = node as HierarchyPointNode<TreeNodeWithState>;
        if (pointNode.x && pointNode.y) {
          node.data.x0 = pointNode.x;
          node.data.y0 = pointNode.y;
        }
      });
      // Now, we find the target node that has children visible.
      // By construction, only one copy of the target node in the tree will have children visible.
      // The target node is the node corresponding to the cell type id of the CellGuideCard.
      let targetNode = data
        .descendants()
        .find(
          (node) =>
            node.data.id.split("__").at(0) === cellTypeId && node.data.children
        ) as HierarchyPointNode<TreeNodeWithState> | undefined;
      // If no target nodes have children, just pick any target node.
      if (!targetNode) {
        targetNode = data
          .descendants()
          .find((node) => node.data.id.split("__").at(0) === cellTypeId) as
          | HierarchyPointNode<TreeNodeWithState>
          | undefined;
      }
      // If the target node is found and its position is known, set the initial transform matrix.
      // This will always be false when in tissue mode.
      if (targetNode) {
        if (targetNode.x !== undefined && targetNode.y !== undefined) {
          setInitialTransformMatrix({
            scaleX: 1,
            scaleY: 1,
            translateX: width / 2 - targetNode.y,
            translateY: height / 2 - targetNode.x,
            skewX: 0,
            skewY: 0,
          });
          setCenteredNodeCoords(true);
        }
      }
    }
  }, [data, cellTypeId, width, height]);

  // Trigger re-render when full screen mode is toggled
  useEffect(() => {
    toggleTriggerRender();
  }, [isFullScreen]);

  // Hover over node tooltip
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

  const yMax = height - defaultMargin.top - defaultMargin.bottom;
  const xMax = width - defaultMargin.left - defaultMargin.right;

  return (
    <div data-testid={CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW}>
      <Global styles={TooltipInPortalStyle} />
      {!!cellTypeId && (
        <TableTitleWrapper>
          <TableTitle>Cell Ontology</TableTitle>
        </TableTitleWrapper>
      )}
      {data && initialTreeState ? (
        <Zoom<SVGSVGElement>
          key={centeredNodeCoords ? "centered" : "initial"}
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
              data-testid={CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER}
            >
              <FullscreenButton
                data-testid={
                  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON
                }
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
                  <div data-testid={CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP}>
                    <b>{tooltipData?.n_cells}</b>
                    {" cells"}
                    {tissueName ? ` in ${tissueName}` : ""}
                    {tooltipData?.n_cells !== tooltipData?.n_cells_rollup && (
                      <>
                        <br />
                        <b>{tooltipData?.n_cells_rollup}</b>
                        {" descendant cells"}
                        {tissueName ? ` in ${tissueName}` : ""}
                      </>
                    )}
                  </div>
                </TooltipInPortal>
              )}
              <StyledSVG
                width={width}
                height={height}
                ref={zoom.containerRef}
                isDragging={zoom.isDragging}
              >
                <RectClipPath id="zoom-clip" width={width} height={height} />
                <rect
                  width={width}
                  height={height}
                  rx={14}
                  fill={backgroundColor}
                />
                <g transform={zoom.toString()}>
                  <Tree<TreeNodeWithState>
                    root={data}
                    size={[yMax, xMax]}
                    nodeSize={NODE_SPACINGS as [number, number]}
                  >
                    {(tree) => (
                      <Group top={defaultMargin.top} left={defaultMargin.left}>
                        <AnimatedLinks tree={tree} duration={duration} />
                        <AnimatedNodes
                          tree={tree}
                          cellTypeId={cellTypeId}
                          duration={duration}
                          setDuration={setDuration}
                          toggleTriggerRender={toggleTriggerRender}
                          showTooltip={showTooltip}
                          hideTooltip={hideTooltip}
                        />
                      </Group>
                    )}
                  </Tree>
                </g>
                <Legend width={width} />
              </StyledSVG>
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
            available as it is not a descendant of &quot;animal cell&quot;.
          </TableUnavailableDescription>
        </TableUnavailableContainer>
      )}
    </div>
  );
}

function setIsExpandedInRawTree(
  graph: TreeNodeWithState,
  isExpandedNodes: string[]
) {
  graph.isExpanded = isExpandedNodes.includes(graph.id);
  if (graph.children) {
    for (const child of graph.children) {
      setIsExpandedInRawTree(child, isExpandedNodes);
    }
  }
}

function setCellCountsInRawTree(
  graph: TreeNodeWithState,
  tissueCounts: TissueCountsPerCellType
) {
  const cellTypeId = graph.id.split("__").at(0) ?? graph.id;

  graph.n_cells = tissueCounts[cellTypeId]?.n_cells ?? graph.n_cells;
  graph.n_cells_rollup =
    tissueCounts[cellTypeId]?.n_cells_rollup ?? graph.n_cells_rollup;

  if (graph.children) {
    for (const child of graph.children) {
      setCellCountsInRawTree(child, tissueCounts);
    }
  }
}
