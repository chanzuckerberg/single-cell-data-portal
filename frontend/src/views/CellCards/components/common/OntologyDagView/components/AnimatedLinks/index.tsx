import NodeGroup from "react-move/NodeGroup";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { LinkHorizontal } from "@visx/shape";
import { Group } from "@visx/group";
import { TreeNodeWithState } from "../../common/types";
import { secondaryColor } from "../../common/constants";

interface AnimatedLinksProps {
  tree: HierarchyPointNode<TreeNodeWithState>;
  duration: number;
}

export default function AnimatedLinks({ tree, duration }: AnimatedLinksProps) {
  return (
    <NodeGroup
      data={tree.links()}
      keyAccessor={(d) => `${d.source.data.id}_${d.target.data.id}`}
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
  );
}

function findCollapsedParent(
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
