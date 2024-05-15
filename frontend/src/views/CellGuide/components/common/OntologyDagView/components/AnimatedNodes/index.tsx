import NodeGroup from "react-move/NodeGroup";
import { Group } from "@visx/group";
import { localPoint } from "@visx/event";
import Node from "../Node";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import {
  CellType,
  MarkerGeneStatsByCellType,
  TreeNodeWithState,
} from "../../common/types";
import { useCellTypeMetadata } from "src/common/queries/cellGuide";
import { NODE_SPACINGS, TREE_ANIMATION_DURATION } from "../../common/constants";
import { EVENTS } from "src/common/analytics/events";
import { track } from "src/common/analytics";
import { MouseEventHandler, useState } from "react";
import { useRouter } from "next/router";
import { getCellTypeLink } from "src/views/CellGuide/common/utils";
import { NO_ORGAN_ID } from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";

interface AnimatedNodesProps {
  tree: HierarchyPointNode<TreeNodeWithState>;
  tissueId?: string;
  cellTypeId?: string;
  duration: number;
  setDuration: (duration: number) => void;
  toggleTriggerRender: () => void;
  showTooltip: (args: {
    tooltipLeft: number;
    tooltipTop: number;
    tooltipData: {
      n_cells: number;
      n_cells_rollup: number;
      marker_score?: number;
      me?: number;
      pc?: number;
    };
  }) => void;
  hideTooltip: () => void;
  cellTypesWithMarkerGeneStats: MarkerGeneStatsByCellType | null;
  setCellInfoCellType?: (props: CellType | null) => void;
  setLastNodeClicked: (nodeId: string) => void;
}

interface AnimationState {
  top: number;
  left: number;
  opacity: number;
  timing: { duration: number };
}

interface AnimationNode {
  key: string;
  data: HierarchyPointNode<TreeNodeWithState>;
  state: AnimationState;
}

export default function AnimatedNodes({
  tree,
  cellTypeId,
  tissueId = NO_ORGAN_ID,
  duration,
  setDuration,
  toggleTriggerRender,
  showTooltip,
  hideTooltip,
  cellTypesWithMarkerGeneStats,
  setCellInfoCellType,
  setLastNodeClicked,
}: AnimatedNodesProps) {
  const [timerId, setTimerId] = useState<number | null>(null); // For hover event
  const router = useRouter();
  const handleAnimationEnd = (node: HierarchyPointNode<TreeNodeWithState>) => {
    // Update the starting position of the node to be its current position
    node.data.x0 = node.x;
    node.data.y0 = node.y;
  };
  const { data: cellTypeMetadata } = useCellTypeMetadata() || {};
  const handleNodeLabelClick = ({ cellTypeId, cellTypeName }: CellType) => {
    if (setCellInfoCellType) {
      setCellInfoCellType({ cellTypeId, cellTypeName });
    } else {
      const url = getCellTypeLink({
        tissueId,
        cellTypeId,
      });

      router.push(url);
    }
  };

  const handleMouseOver = (
    event: React.MouseEvent<HTMLDivElement>,
    datum: TreeNodeWithState
  ) => {
    if (!timerId) {
      /**
       * (thuang): Use window.setTimeout instead of setTimeout to avoid
       * NodeJS setTimeout type being used
       */
      const id = window.setTimeout(() => {
        track(EVENTS.CG_TREE_NODE_HOVER, {
          cell_type: datum.name,
        });
      }, 2 * 1000);
      setTimerId(id);
    }
    const ownerSVGElement = findSVGParent(event.target as Node);
    if (event.target instanceof HTMLDivElement && ownerSVGElement) {
      const coords = localPoint(ownerSVGElement, event);
      if (coords) {
        const marker_score =
          cellTypesWithMarkerGeneStats?.[datum.id]?.marker_score;
        const me = cellTypesWithMarkerGeneStats?.[datum.id]?.me;
        const pc = cellTypesWithMarkerGeneStats?.[datum.id]?.pc;
        showTooltip({
          tooltipLeft: coords.x,
          tooltipTop: coords.y,
          tooltipData: {
            n_cells: datum.n_cells,
            n_cells_rollup: datum.n_cells_rollup,
            marker_score,
            me,
            pc,
          },
        });
      }
    }
  };

  return (
    <NodeGroup
      data={tree.descendants()}
      keyAccessor={(d: HierarchyPointNode<TreeNodeWithState>) => d.data.id}
      start={(node: HierarchyPointNode<TreeNodeWithState>) => {
        // when the animation ends, update the start position of the
        // node to be its current position
        setTimeout(() => handleAnimationEnd(node), duration);
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
      enter={(node: HierarchyPointNode<TreeNodeWithState>) => {
        return {
          top: [node.x],
          left: [node.y],
          opacity: [1],
          timing: { duration },
        };
      }}
      update={(node: HierarchyPointNode<TreeNodeWithState>) => {
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
      {(nodes) =>
        renderNodes(nodes, cellTypesWithMarkerGeneStats, handleNodeLabelClick)
      }
    </NodeGroup>
  );

  function renderNode(
    animatedNode: {
      key: string;
      data: HierarchyPointNode<TreeNodeWithState>;
      state: AnimationState;
    },
    cellTypesWithMarkerGeneStats: MarkerGeneStatsByCellType | null,
    handleNodeLabelClick: (props: CellType) => void
  ) {
    const { key, data: node, state } = animatedNode;

    const isInCorpus =
      (node.data.id.split("__")[0]?.replace("_", ":") ?? "") in
      (cellTypeMetadata ?? {});

    return (
      <Node
        key={`${key}-node`}
        isInCorpus={isInCorpus}
        animationKey={key}
        node={node}
        cellTypesWithMarkerGeneStats={cellTypesWithMarkerGeneStats}
        isTargetNode={cellTypeId === node.data.id.split("__")[0]}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
        maxWidth={NODE_SPACINGS[1] - 20}
        left={state.left}
        top={state.top}
        opacity={state.opacity}
        handleClick={
          handleNodeClick as unknown as MouseEventHandler<HTMLDivElement>
        }
        handleNodeLabelClick={handleNodeLabelClick}
      />
    );

    function handleNodeClick() {
      if (duration === 0) {
        setDuration(TREE_ANIMATION_DURATION);
      }

      if (node.data.id.startsWith("dummy-child")) {
        track(EVENTS.CG_TREE_NODE_CLICKED, {
          cell_type: "multiple cell types",
        });
        if (node.parent) {
          node.parent.data.showAllChildren = true;
        }
      } else {
        track(EVENTS.CG_TREE_NODE_CLICKED, { cell_type: node.data.name });
      }

      node.data.isExpanded = !node.data.isExpanded;
      if (!node.data.isExpanded) {
        node.data.showAllChildren = false;
        collapseAllDescendants(node);
      }
      toggleTriggerRender();
      setLastNodeClicked(node.data.id);
    }
  }

  function handleMouseOut() {
    hideTooltip();
    if (timerId) {
      /**
       * (thuang): Use window.clearTimeout instead of clearTimeout to avoid
       * NodeJS clearTimeout type being used
       */
      window.clearTimeout(timerId);
      setTimerId(null);
    }
  }

  function renderNodes(
    nodes: AnimationNode[],
    cellTypesWithMarkerGeneStats: AnimatedNodesProps["cellTypesWithMarkerGeneStats"],
    handleNodeLabelClick: (props: CellType) => void
  ) {
    return (
      <Group>
        {nodes.map((node) =>
          renderNode(node, cellTypesWithMarkerGeneStats, handleNodeLabelClick)
        )}
      </Group>
    );
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

function findSVGParent(node: Node | null): SVGElement | undefined {
  while (node) {
    if (node instanceof SVGElement) {
      return node;
    }
    node = node.parentNode;
  }
  return undefined;
}
