import { MouseEventHandler } from "react";
import { TreeNodeWithState } from "../../../../common/types";
import {
  tertiaryColor,
  primaryColor,
  highlightColor,
  smallSize,
  largeSize,
  markerGeneModeColor,
} from "../../../../common/constants";
import { StyledRect, StyledCircle } from "./style";
import { CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID } from "src/views/CellGuide/components/common/OntologyDagView/components/Node/components/RectOrCircle/constants";

const DUMMY_CHILD = "dummy-child";

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
  inMarkerGeneMode: boolean;
  hasMarkerGene: boolean;
}

export default function RectOrCircle({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  animationKey,
  inMarkerGeneMode,
  hasMarkerGene,
}: RectOrCircleProps) {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
  if (isTargetNode) {
    color = highlightColor;
  }

  const size = node.n_cells === 0 ? smallSize : largeSize;
  if (node.id.startsWith(DUMMY_CHILD)) {
    color = primaryColor;
  }

  if (inMarkerGeneMode && hasMarkerGene) {
    color = markerGeneModeColor;
  } else if (inMarkerGeneMode) {
    color = tertiaryColor;
  }

  const onMouseOver = node.id.startsWith(DUMMY_CHILD)
    ? undefined
    : (event: React.MouseEvent<SVGElement>) => {
        handleMouseOver(event, node);
      };
  const onMouseOut = node.name.startsWith(DUMMY_CHILD)
    ? undefined
    : handleMouseOut;
  return node?.children?.length || node.id.startsWith(DUMMY_CHILD) ? (
    <StyledCircle
      data-testid={`${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${node.id}-has-children-isTargetNode=${isTargetNode}`}
      r={size}
      fill={color}
      key={animationKey}
      onClick={handleClick}
      strokeWidth={0.5}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
    />
  ) : (
    <StyledRect
      data-testid={`${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${node.id}-no-children-isTargetNode=${isTargetNode}`}
      height={size * 2}
      width={size * 2}
      y={-size}
      x={-size}
      key={animationKey}
      fill={color}
      onClick={handleClick}
      strokeWidth={0.5}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
    />
  );
}
