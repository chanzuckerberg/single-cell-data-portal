import { MouseEventHandler } from "react";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";
import { TreeNodeWithState } from "../../common/types";
import Text from "./components/Text";
import RectOrCircle from "./components/RectOrCircle";
import { StyledGroup } from "./style";

export const CELL_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL =
  "cell-card-ontology-dag-view-clickable-text-label";

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

export default function Node({
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
    <StyledGroup
      top={top}
      left={left}
      key={animationKey}
      opacity={opacity}
      onClick={handleClick}
    >
      {isInCorpus ? (
        <g
          data-testid={`${CELL_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-${node.data.id}`}
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
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
      />
    </StyledGroup>
  );
}
