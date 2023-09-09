import { MouseEventHandler } from "react";
import { HierarchyPointNode } from "@visx/hierarchy/lib/types";
import { useRouter } from "next/router";
import { ROUTES } from "src/common/constants/routes";
import { TreeNodeWithState } from "../../common/types";
import Text from "./components/Text";
import RectOrCircle from "./components/RectOrCircle";
import { StyledGroup } from "./style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL } from "src/views/CellGuide/components/common/OntologyDagView/components/Node/constants";

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
  cellTypesWithMarkerGenes: string[] | null;
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
  cellTypesWithMarkerGenes,
}: NodeProps) {
  const router = useRouter();

  // text labels should only collapse/expand node for dummy nodes
  const onClick = node.data.id.startsWith("dummy-child")
    ? handleClick
    : undefined;
  const textCursor = node.data.id.startsWith("dummy-child")
    ? "pointer"
    : "default";
  const inMarkerGeneMode = !!cellTypesWithMarkerGenes;
  const hasMarkerGene = !!cellTypesWithMarkerGenes?.includes(node.data.id);
  return (
    <StyledGroup top={top} left={left} key={animationKey} opacity={opacity}>
      {isInCorpus ? (
        <g
          data-testid={`${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-${node.data.id}`}
          onClick={() => {
            track(EVENTS.CG_TREE_CELL_TYPE_CLICKED, {
              cell_type: node.data.name,
            });
            router.push(
              `${ROUTES.CELL_GUIDE}/${node.data.id
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
        <g onClick={onClick}>
          <Text
            isInCorpus={isInCorpus}
            name={node.data.name}
            maxWidth={maxWidth}
            cursor={textCursor}
          />
        </g>
      )}
      <RectOrCircle
        animationKey={`${animationKey}-rect-or-circle`}
        node={node.data}
        handleClick={handleClick}
        isTargetNode={isTargetNode}
        handleMouseOver={handleMouseOver}
        handleMouseOut={handleMouseOut}
        inMarkerGeneMode={inMarkerGeneMode}
        hasMarkerGene={hasMarkerGene}
      />
    </StyledGroup>
  );
}
