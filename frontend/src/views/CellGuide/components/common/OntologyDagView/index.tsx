import React, { useMemo, useEffect, useState } from "react";
import { Group } from "@visx/group";
import { Global } from "@emotion/react";
import { useTooltip, useTooltipInPortal } from "@visx/tooltip";
import { Tree, hierarchy } from "@visx/hierarchy";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import FullscreenIcon from "@mui/icons-material/Fullscreen";
import FullscreenExitIcon from "@mui/icons-material/FullscreenExit";
import { CellOntologyTreeResponse as TreeNode } from "src/common/queries/cellGuide";
import {
  TableTitleWrapper,
  TableTitle,
} from "../../CellGuideCard/components/common/style";
import { Zoom } from "@visx/zoom";
import { RectClipPath } from "@visx/clip-path";
import {
  CellOntologyTreeStateResponse,
  TissueCountsPerCellType,
  useCellOntologyTree,
  useCellOntologyTreeStateCellType,
  useCellOntologyTreeStateTissue,
  useMarkerGenePresenceQuery,
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
  StyledButtonIcon,
  RightAligned,
} from "./style";
import { useFullScreen } from "../FullScreenProvider";
import {
  defaultMargin,
  backgroundColor,
  NODE_SPACINGS,
} from "./common/constants";
import { TreeNodeWithState } from "./common/types";
import Legend from "./components/Legend";
import AnimatedNodes from "./components/AnimatedNodes";
import AnimatedLinks from "./components/AnimatedLinks";
import {
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_DEACTIVATE_MARKER_GENE_MODE,
  MINIMUM_NUMBER_OF_HIDDEN_CHILDREN_FOR_DUMMY_NODE,
} from "src/views/CellGuide/components/common/OntologyDagView/constants";
import {
  ALL_TISSUES,
  TISSUE_AGNOSTIC,
} from "../../CellGuideCard/components/MarkerGeneTables/constants";

interface TreeProps {
  skinnyMode?: boolean;
  inputWidth: number;
  inputHeight: number;
  selectedGene?: string;
  selectedOrganism?: string;
  cellTypeId?: string;
  tissueId: string;
  tissueName: string;
  selectGene?: (gene: string) => void;
}

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
  inputWidth,
  inputHeight,
  selectedGene,
  selectGene,
  selectedOrganism,
}: TreeProps) {
  const [width, setWidth] = useState(inputWidth);
  const [height, setHeight] = useState(inputHeight);

  const [initialTransformMatrix, setInitialTransformMatrix] = useState<
    typeof initialTransformMatrixDefault
  >(initialTransformMatrixDefault);

  // This toggler is used for centering the Zoom component on the target cell type.
  // It triggers a re-render of the Zoom component so that updates to the initialTransformMatrix
  // take effect.
  const [centeredNodeCoords, setCenteredNodeCoords] = useState<boolean>(false);

  const {
    isFullScreen,
    screenDimensions,
    enableFullScreen,
    disableFullScreen,
  } = useFullScreen();

  // Handle the resizing of the ontology view when full screen mode is toggled
  useEffect(() => {
    let newWidth = inputWidth;
    let newHeight = inputHeight;
    if (screenDimensions.width > 0 && isFullScreen) {
      newWidth = screenDimensions.width;
    }
    if (screenDimensions.height > 0 && isFullScreen) {
      newHeight = screenDimensions.height;
    }
    setWidth(newWidth);
    setHeight(newHeight);
  }, [screenDimensions, isFullScreen, inputWidth, inputHeight]);

  const { data: markerGenePresence, isLoading: isLoadingMarkerGenePresence } =
    useMarkerGenePresenceQuery();

  const selectedTissue =
    tissueName === TISSUE_AGNOSTIC ? ALL_TISSUES : tissueName;

  const cellTypesWithMarkerGeneStats: {
    [cellTypeId: string]: {
      me: number;
      pc: number;
      marker_score: number;
    };
  } | null = useMemo(() => {
    if (
      isLoadingMarkerGenePresence ||
      !markerGenePresence ||
      !selectedGene ||
      !selectedOrganism
    )
      return null;
    return (
      markerGenePresence?.[selectedGene]?.[selectedOrganism]?.[
        selectedTissue
      ]?.reduce((acc, markerGeneStats) => {
        const { cell_type_id, ...rest } = markerGeneStats;
        return {
          ...acc,
          [`${cell_type_id}__0`]: rest,
        };
      }, {}) ?? {}
    );
  }, [
    markerGenePresence,
    isLoadingMarkerGenePresence,
    selectedGene,
    selectedOrganism,
    selectedTissue,
  ]);

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

  const { data: initialTreeStateCell } = useCellOntologyTreeStateCellType(
    cellTypeId ?? ""
  );
  const { data: initialTreeStateTissue } = useCellOntologyTreeStateTissue(
    tissueId ?? ""
  );

  const parentMap = useMemo(() => {
    if (!rawTree) return null;
    const queue: TreeNodeWithState[] = [rawTree];
    const map = new Map<string, string>();
    while (queue.length > 0) {
      const node = queue.pop() as TreeNodeWithState;
      if (node.children) {
        for (const child of node.children) {
          map.set(child.id, node.id);
          queue.push(child);
        }
      }
    }
    return map;
  }, [rawTree]);

  const initialTreeState: CellOntologyTreeStateResponse | undefined =
    useMemo(() => {
      let initialTreeState;
      if (cellTypesWithMarkerGeneStats && parentMap) {
        const cellTypesWithMarkerGene = Object.keys(
          cellTypesWithMarkerGeneStats
        );
        initialTreeState = getInitialStateForSelectedGene(
          rawTree as TreeNode,
          parentMap,
          cellTypesWithMarkerGene
        );
      } else if (initialTreeStateCell && initialTreeStateTissue) {
        // When both cell and tissue tree states are available, inject the tissue tree counts
        // into the cell tree state.
        initialTreeState = {
          ...initialTreeStateCell,
          tissueCounts: initialTreeStateTissue.tissueCounts,
        };
      } else if (initialTreeStateCell) {
        initialTreeState = initialTreeStateCell;
      } else if (initialTreeStateTissue) {
        initialTreeState = initialTreeStateTissue;
      }
      return initialTreeState;
    }, [
      initialTreeStateCell,
      initialTreeStateTissue,
      rawTree,
      parentMap,
      cellTypesWithMarkerGeneStats,
    ]);

  // Populate the tree data structure nodes with "isExpanded".
  const treeData: TreeNodeWithState | null = useMemo(() => {
    const newTree = rawTree ? JSON.parse(JSON.stringify(rawTree)) : null;
    if (newTree && initialTreeState) {
      if (initialTreeState.isExpandedNodes) {
        setIsExpandedInRawTree(newTree, initialTreeState.isExpandedNodes);
      }
      if (initialTreeState.tissueCounts) {
        setCellCountsInRawTreeForTissueCardOntologyView(
          newTree,
          initialTreeState.tissueCounts
        );
        deleteNodesWithNoDescendantsForTissueCardOntologyView(newTree);
      }
    }
    return newTree;
  }, [rawTree, initialTreeState]);

  // Create the tree data structure
  const data = useMemo(() => {
    if (!treeData) return null;
    return hierarchy(treeData, (d) => {
      if (d.isExpanded && d.children && initialTreeState) {
        const notShownWhenExpandedNodes =
          initialTreeState.notShownWhenExpandedNodes;
        /**
         * If a node is a key in `notShownWhenExpandedNodes`, then it is a parent to nodes
         * that should be collapsed. The text label of this "dummy" node is the number of hidden terms,
         * e.g. 52 cell types.
         * `showAllChildren` is a flag that is set to `true` when the user clicks on a collapsed node.
         * It indicates that all of the children of the node should be shown.
         */

        const hiddenChildren = d.children.filter(
          (child) =>
            !d.showAllChildren &&
            notShownWhenExpandedNodes[d.id]?.includes(child.id)
        );

        const newChildren = d.children.filter(
          (child) =>
            d.showAllChildren ||
            !notShownWhenExpandedNodes[d.id]?.includes(child.id)
        );

        if (
          hiddenChildren.length <=
          MINIMUM_NUMBER_OF_HIDDEN_CHILDREN_FOR_DUMMY_NODE
        ) {
          newChildren.push(...hiddenChildren);
        } else if (hiddenChildren.length > 0) {
          newChildren.push({
            id: `dummy-child-${d.id}`,
            name: `${hiddenChildren.length} cell types`,
            n_cells: 0,
            n_cells_rollup: 0,
            isExpanded: false,
          });
        }
        return newChildren;
      }
    });
    /**
     * TODO: Even though `triggerRender` is not used in the memo, it is used here to force a re-render
     * when the node components update the tree data in-place when users collapse and expand nodes.
     * This is a known anti-pattern and will be addressed in later work.
     * See this ticket: https://github.com/chanzuckerberg/single-cell-data-portal/issues/5478
     */
    // eslint-disable-next-line react-hooks/exhaustive-deps
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

  // Hover over node tooltip
  const {
    tooltipData,
    tooltipLeft,
    tooltipTop,
    tooltipOpen,
    showTooltip,
    hideTooltip,
  } = useTooltip<{
    n_cells: number;
    n_cells_rollup: number;
    marker_score?: number;
    me?: number;
    pc?: number;
  }>();

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
      {data && initialTreeState && <Legend selectedGene={selectedGene} />}
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
              <RightAligned>
                <FullscreenButton
                  data-testid={
                    CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON
                  }
                  onClick={isFullScreen ? disableFullScreen : enableFullScreen}
                >
                  {isFullScreen ? <FullscreenExitIcon /> : <FullscreenIcon />}
                </FullscreenButton>
                {selectedGene && (
                  <StyledButtonIcon
                    aria-label={`deactivate ${selectedGene} marker gene ontology view mode`}
                    data-testid={
                      CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_DEACTIVATE_MARKER_GENE_MODE
                    }
                    sdsIcon="eyeClosed"
                    sdsSize="small"
                    sdsType="secondary"
                    onClick={() => selectGene && selectGene(selectedGene)}
                  />
                )}
              </RightAligned>

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
                    {tissueName ? ` in ${selectedTissue.toLowerCase()}` : ""}
                    {tooltipData?.n_cells !== tooltipData?.n_cells_rollup && (
                      <>
                        <br />
                        <b>{tooltipData?.n_cells_rollup}</b>
                        {" descendant cells"}
                        {tissueName
                          ? ` in ${selectedTissue.toLowerCase()}`
                          : ""}
                      </>
                    )}
                    {tooltipData &&
                      tooltipData.marker_score &&
                      tooltipData.me &&
                      tooltipData.pc && (
                        <>
                          <br />
                          <br />
                          <b>{selectedGene} stats</b>
                          <br />
                          {"Marker score: "}
                          <b>{tooltipData.marker_score.toFixed(2)}</b>

                          <br />
                          {"Expression score: "}
                          <b>{tooltipData.me.toFixed(2)}</b>
                          <br />
                          {"% of cells: "}
                          <b>{(tooltipData.pc * 100).toFixed(2)}</b>
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
                          cellTypesWithMarkerGeneStats={
                            cellTypesWithMarkerGeneStats
                          }
                        />
                      </Group>
                    )}
                  </Tree>
                </g>
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

function setCellCountsInRawTreeForTissueCardOntologyView(
  graph: TreeNodeWithState,
  tissueCounts: TissueCountsPerCellType
) {
  const cellTypeId = graph.id.split("__").at(0) ?? graph.id;
  // (alec) this makes it so that cell types that aren't present in the tissue will be removed
  // by deleteNodesWithNoDescendantsForTissueCardOntologyView
  graph.n_cells = tissueCounts[cellTypeId]?.n_cells ?? 0;
  graph.n_cells_rollup = tissueCounts[cellTypeId]?.n_cells_rollup ?? 0;

  if (graph.children) {
    for (const child of graph.children) {
      setCellCountsInRawTreeForTissueCardOntologyView(child, tissueCounts);
    }
  }
}

function deleteNodesWithNoDescendantsForTissueCardOntologyView(
  graph: TreeNodeWithState
): void {
  if (graph.children) {
    graph.children = graph.children.filter((child) => {
      deleteNodesWithNoDescendantsForTissueCardOntologyView(child);
      return child.n_cells_rollup !== 0;
    });
  }
}

function getInitialStateForSelectedGene(
  rawTree: TreeNodeWithState,
  parentMap: Map<string, string>,
  validCellTypes: string[]
): CellOntologyTreeStateResponse {
  // first get all the ancestors of each node in validCellTypes
  const ancestors = new Set<string>();
  for (const cellType of validCellTypes) {
    let current: string | undefined = cellType;
    while (current) {
      ancestors.add(current);
      current = parentMap.get(current);
    }
  }
  const validNodes = Array.from(ancestors);
  const tree: OntologyTree = JSON.parse(JSON.stringify(rawTree));

  truncateGraph(tree, validNodes);

  const isExpandedNodes = Array.from(new Set(getExpandedData(tree)));
  const notShownWhenExpandedNodes = getShownData(tree);

  let notShownWhenExpanded: ShownData = {};
  for (const i of notShownWhenExpandedNodes) {
    notShownWhenExpanded = { ...notShownWhenExpanded, ...i };
  }
  return {
    isExpandedNodes,
    notShownWhenExpandedNodes: notShownWhenExpanded,
  };
}

interface OntologyTree extends TreeNode {
  invalid_children_ids?: string[];
  parent?: string;
}

function filterChildren(
  graph: OntologyTree,
  validNodes: string[]
): OntologyTree[] {
  const newChildren: OntologyTree[] = [];
  const invalidChildrenIds: string[] = [];

  for (const child of graph.children ?? []) {
    if (validNodes.includes(child.id)) {
      newChildren.push(child);
    } else {
      invalidChildrenIds.push(child.id);
    }
  }

  if (invalidChildrenIds.length > 0) {
    newChildren.push({
      id: "",
      name: "",
      n_cells: 0,
      n_cells_rollup: 0,
      invalid_children_ids: invalidChildrenIds,
      parent: graph.id,
    });
  }

  return newChildren;
}

function truncateGraph(graph: OntologyTree, validNodes: string[]): void {
  const children = filterChildren(graph, validNodes);

  if (children.length === 0) {
    graph.children = undefined;
  } else {
    graph.children = children;
  }

  for (const child of graph.children ?? []) {
    if (child.id !== "") {
      truncateGraph(child, validNodes);
    }
  }
}

function getExpandedData(
  ontologyGraph: OntologyTree,
  isExpandedNodes: string[] = []
): string[] {
  if (ontologyGraph.children) {
    isExpandedNodes.push(ontologyGraph.id);
    for (const child of ontologyGraph.children) {
      getExpandedData(child, isExpandedNodes);
    }
  }
  return isExpandedNodes;
}

interface ShownData {
  [key: string]: string[];
}

function getShownData(
  graph: OntologyTree,
  notShownWhenExpandedNodes: ShownData[] = []
): ShownData[] {
  if (graph.children) {
    for (const child of graph.children) {
      if (
        child.id === "" &&
        child.invalid_children_ids &&
        child.invalid_children_ids.length > 0
      ) {
        notShownWhenExpandedNodes.push({
          [child.parent as string]: Array.from(
            new Set(child.invalid_children_ids)
          ),
        });
      } else {
        getShownData(child, notShownWhenExpandedNodes);
      }
    }
  }
  return notShownWhenExpandedNodes;
}
