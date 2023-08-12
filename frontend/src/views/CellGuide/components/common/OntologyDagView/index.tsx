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
  InitialCellOntologyTreeStateResponse,
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

// Both cellTypeId and tissueId (with tissueName) are mandatory.
interface BothCellAndTissueIdProps extends BaseTreeProps {
  cellTypeId: string;
  tissueId: string;
  tissueName: string;
}

type TreeProps = BothCellAndTissueIdProps | CellTypeIdProps | TissueIdProps;

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

  const { data: initialTreeStateCell } = useCellOntologyTreeState(
    cellTypeId ?? ""
  );
  const { data: initialTreeStateTissue } = useCellOntologyTreeStateTissue(
    tissueId ?? ""
  );
  const initialTreeState: InitialCellOntologyTreeStateResponse | undefined =
    useMemo(() => {
      let initialTreeState;
      if (initialTreeStateCell && initialTreeStateTissue) {
        // Intersection of isExpandedNodes
        const isExpandedIntersection =
          initialTreeStateCell.isExpandedNodes.filter((node) =>
            initialTreeStateTissue.isExpandedNodes.includes(node)
          );

        // Union of notShownWhenExpandedNodes
        const notShownWhenExpandedUnion: {
          [key: string]: string[];
        } = { ...initialTreeStateCell.notShownWhenExpandedNodes };

        for (const key of Object.keys(
          initialTreeStateTissue.notShownWhenExpandedNodes
        )) {
          if (notShownWhenExpandedUnion[key]) {
            // If the key exists in both, merge the arrays and de-duplicate
            notShownWhenExpandedUnion[key] = Array.from(
              new Set([
                ...notShownWhenExpandedUnion[key],
                ...initialTreeStateTissue.notShownWhenExpandedNodes[key],
              ])
            );
          } else {
            // If the key exists only in the tissue data, add it to the union
            notShownWhenExpandedUnion[key] =
              initialTreeStateTissue.notShownWhenExpandedNodes[key];
          }
        }

        initialTreeState = {
          isExpandedNodes: isExpandedIntersection,
          notShownWhenExpandedNodes: notShownWhenExpandedUnion,
          tissueCounts: initialTreeStateTissue.tissueCounts, // Only specified for tissueId
        };
      } else if (initialTreeStateCell) {
        initialTreeState = initialTreeStateCell;
      } else if (initialTreeStateTissue) {
        initialTreeState = initialTreeStateTissue;
      }
      return initialTreeState;
    }, [initialTreeStateCell, initialTreeStateTissue]);

  // Populate the tree data structure nodes with "isExpanded".
  const treeData: TreeNodeWithState | null = useMemo(() => {
    const newTree = rawTree ? JSON.parse(JSON.stringify(rawTree)) : null;
    if (newTree && initialTreeState) {
      if (initialTreeState.isExpandedNodes) {
        setIsExpandedInRawTree(newTree, initialTreeState.isExpandedNodes);
      }
      if (initialTreeState.tissueCounts) {
        setCellCountsInRawTree(newTree, initialTreeState.tissueCounts);
        deleteNodesWithNoDescendants(newTree);
      }
    }
    return newTree;
  }, [rawTree, initialTreeState]);

  // Create the tree data structure
  const data = useMemo(() => {
    if (!treeData) return null;
    return hierarchy(treeData, (d) => {
      if (d.isExpanded && d.children && initialTreeState) {
        const newChildren: TreeNodeWithState[] = [];
        const notShownWhenExpandedNodes =
          initialTreeState.notShownWhenExpandedNodes;
        /**
         * If a node is a key in `notShownWhenExpandedNodes`, then it is a parent to nodes
         * that should be collapsed. The text label of this "dummy" node is the number of hidden terms,
         * e.g. 52 cell types.
         * `showAllChildren` is a flag that is set to `true` when the user clicks on a collapsed node.
         * It indicates that all of the children of the node should be shown.
         */

        for (const child of d.children) {
          if (
            d.showAllChildren ||
            !notShownWhenExpandedNodes[d.id]?.includes(child.id)
          ) {
            newChildren.push(child);
          }
        }

        const numHiddenChildren = d.children.length - newChildren.length;
        if (numHiddenChildren > 0) {
          newChildren.push({
            id: `dummy-child-${d.id}`,
            name: `${numHiddenChildren} cell types`,
            n_cells: 0,
            n_cells_rollup: 0,
            isExpanded: false,
          });
        }
        return newChildren;
      }
    });
  }, [treeData, triggerRender, initialTreeState]);

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
      if (
        targetNode &&
        targetNode.x !== undefined &&
        targetNode.y !== undefined
      ) {
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

  graph.n_cells = tissueCounts[cellTypeId]?.n_cells ?? 0;
  graph.n_cells_rollup = tissueCounts[cellTypeId]?.n_cells_rollup ?? 0;

  if (graph.children) {
    for (const child of graph.children) {
      setCellCountsInRawTree(child, tissueCounts);
    }
  }
}

function deleteNodesWithNoDescendants(graph: TreeNodeWithState): void {
  if (graph.children) {
    graph.children = graph.children.filter((child) => {
      deleteNodesWithNoDescendants(child);
      return child.n_cells_rollup !== 0;
    });
  }
}
