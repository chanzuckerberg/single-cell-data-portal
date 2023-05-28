import NodeGroup from "react-move/NodeGroup";
import { Group } from "@visx/group";
import { localPoint } from "@visx/event";
import Node from "../Node";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { TreeNodeWithState } from "../../common/types";
import { useCellTypesById } from "src/common/queries/cellCards";
import { NODE_SPACINGS, TREE_ANIMATION_DURATION } from "../../common/constants";

interface AnimatedNodesProps {
  tree: HierarchyPointNode<TreeNodeWithState>;
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
    };
  }) => void;
  hideTooltip: () => void;
}
export default function AnimatedNodes({
  tree,
  cellTypeId,
  duration,
  setDuration,
  toggleTriggerRender,
  showTooltip,
  hideTooltip,
}: AnimatedNodesProps) {
  const cellTypesById = useCellTypesById() || {};
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
    <NodeGroup
      data={tree.descendants()}
      keyAccessor={(d: HierarchyPointNode<TreeNodeWithState>) => d.data.id}
      start={(node: HierarchyPointNode<TreeNodeWithState>) => {
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
      {(nodes) => {
        return (
          <Group>
            {nodes.map(({ key, data: node, state }) => {
              const isInCorpus =
                (node.data.id.split("__").at(0)?.replace("_", ":") ?? "") in
                cellTypesById;
                
              return (
                <Node
                  key={`${key}-node`}
                  isInCorpus={isInCorpus}
                  animationKey={key}
                  node={node}
                  isTargetNode={cellTypeId === node.data.id.split("__").at(0)}
                  handleMouseOver={handleMouseOver}
                  handleMouseOut={hideTooltip}
                  maxWidth={NODE_SPACINGS[1] - 20}
                  left={state.left}
                  top={state.top}
                  opacity={state.opacity}
                  handleClick={() => {
                    if (duration === 0) {
                      setDuration(TREE_ANIMATION_DURATION);
                    }
                    if (node.data.id.startsWith("dummy-child")) {
                      if (node.parent) {
                        node.parent.data.showAllChildren = true;
                      }
                    }
                    node.data.isExpanded = !node.data.isExpanded;
                    if (!node.data.isExpanded) {
                      node.data.showAllChildren = false;
                      collapseAllDescendants(node);
                    }
                    toggleTriggerRender();
                  }}
                />
              );
            })}
          </Group>
        );
      }}
    </NodeGroup>
  );
}

function collapseAllDescendants(node: HierarchyPointNode<TreeNodeWithState>) {
  if (node.children) {
    for (const child of node.children) {
      child.data.isExpanded = false;
      collapseAllDescendants(child);
    }
  }
}
