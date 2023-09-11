import { MouseEventHandler } from "react";
import { MarkerGeneStats, TreeNodeWithState } from "../../../../common/types";
import {
  tertiaryColor,
  primaryColor,
  highlightColor,
  smallSize,
  largeSize,
  markerGeneModeColor,
} from "../../../../common/constants";
import { StyledRect, StyledCircle, StyledPie } from "./style";
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
  cellTypeMarkerGeneStats: MarkerGeneStats | undefined;
  inMarkerGeneMode: boolean;
}

export default function RectOrCircle({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  animationKey,
  inMarkerGeneMode,
  cellTypeMarkerGeneStats,
}: RectOrCircleProps) {
  const hasMarkerGene = !!cellTypeMarkerGeneStats;
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
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
  if (isTargetNode) {
    color = highlightColor;
  }
  const onMouseOver = node.id.startsWith(DUMMY_CHILD)
    ? undefined
    : (event: React.MouseEvent<SVGElement>) => {
        handleMouseOver(event, node);
      };
  const onMouseOut = node.name.startsWith(DUMMY_CHILD)
    ? undefined
    : handleMouseOut;
  const sizeScaler = cellTypeMarkerGeneStats ? 0.8 : 1;
  return node?.children?.length ||
    node.id.startsWith(DUMMY_CHILD) ||
    cellTypeMarkerGeneStats ? (
    <g>
      <StyledCircle
        data-testid={`${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${node.id}-has-children-isTargetNode=${isTargetNode}`}
        r={size * sizeScaler}
        fillColor={
          cellTypeMarkerGeneStats ? cellTypeMarkerGeneStats.marker_score : color
        }
        key={animationKey}
        onClick={handleClick}
        strokeWidth={0.5}
        onMouseOver={onMouseOver}
        onMouseOut={onMouseOut}
      />
      {cellTypeMarkerGeneStats && (
        <foreignObject
          x={`-${size}`}
          y={`-${size}`}
          width={`${size * 2}`}
          height={`${size * 2}`}
        >
          <StyledPie
            data-testid={`${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${node.id}-has-children-isTargetNode=${isTargetNode}`}
            key={animationKey}
            onClick={
              handleClick as MouseEventHandler<HTMLDivElement> | undefined
            }
            onMouseOver={
              onMouseOver as MouseEventHandler<HTMLDivElement> | undefined
            }
            onMouseOut={onMouseOut}
            degree={cellTypeMarkerGeneStats.pc * 360}
            size={size * 2}
            fill={cellTypeMarkerGeneStats.marker_score}
          />
        </foreignObject>
      )}
    </g>
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
