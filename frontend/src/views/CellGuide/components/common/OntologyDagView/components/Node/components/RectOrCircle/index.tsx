import { MouseEventHandler } from "react";
import { MarkerGeneStats, TreeNodeWithState } from "../../../../common/types";
import {
  tertiaryColor,
  primaryColor,
  highlightColor,
  smallSize,
  largeSize,
  markerGeneModeColor,
  nodePieChartCircleScaler,
} from "../../../../common/constants";
import { StyledRect } from "./style";
import { StyledCircle, StyledPie } from "../../../../common/style";
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
const getColor = (
  node: TreeNodeWithState,
  inMarkerGeneMode: boolean,
  hasMarkerGene: boolean,
  isTargetNode: boolean
) => {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
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

  return color;
};

const getMouseOverHandler = (
  node: TreeNodeWithState,
  handleMouseOver: RectOrCircleProps["handleMouseOver"]
) => {
  return node.id.startsWith(DUMMY_CHILD)
    ? undefined
    : (event: React.MouseEvent<SVGElement>) => {
        handleMouseOver(event, node);
      };
};

const getMouseOutHandler = (
  node: TreeNodeWithState,
  handleMouseOut: RectOrCircleProps["handleMouseOut"]
) => {
  return node.name.startsWith(DUMMY_CHILD) ? undefined : handleMouseOut;
};

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
  const color = getColor(node, inMarkerGeneMode, hasMarkerGene, isTargetNode);
  const size = (node.n_cells === 0 ? smallSize : largeSize) * 2;
  const onMouseOver = getMouseOverHandler(node, handleMouseOver);
  const onMouseOut = getMouseOutHandler(node, handleMouseOut);
  const sizeScaler = cellTypeMarkerGeneStats ? nodePieChartCircleScaler : 1;
  const word = node?.children?.length ? "has" : "no";
  const dataTestId = `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${node.id}-${word}-children-isTargetNode=${isTargetNode}`;
  return node?.children?.length ||
    node.id.startsWith(DUMMY_CHILD) ||
    cellTypeMarkerGeneStats ? (
    <g>
      {cellTypeMarkerGeneStats ? (
        <foreignObject
          data-testid={dataTestId}
          x={-size / 2}
          y={-size / 2}
          width={size}
          height={size}
          onClick={handleClick}
          onMouseOver={onMouseOver}
          onMouseOut={onMouseOut}
        >
          <StyledPie
            key={animationKey}
            degree={cellTypeMarkerGeneStats.pc * 360}
            size={size}
            fill={cellTypeMarkerGeneStats.marker_score}
          >
            <StyledCircle
              size={size * sizeScaler}
              fill={
                cellTypeMarkerGeneStats
                  ? cellTypeMarkerGeneStats.marker_score
                  : color
              }
              key={animationKey}
              opacity={0.5}
            />
          </StyledPie>
        </foreignObject>
      ) : (
        <foreignObject
          data-testid={dataTestId}
          x={-size / 2}
          y={-size / 2}
          width={size}
          height={size}
          onClick={handleClick}
          onMouseOver={onMouseOver}
          onMouseOut={onMouseOut}
        >
          <StyledCircle
            size={size * sizeScaler}
            fill={color}
            key={animationKey}
          />
        </foreignObject>
      )}
    </g>
  ) : (
    <StyledRect
      data-testid={dataTestId}
      height={size}
      width={size}
      y={-size / 2}
      x={-size / 2}
      key={animationKey}
      fill={color}
      onClick={handleClick}
      strokeWidth={0.5}
      onMouseOver={onMouseOver}
      onMouseOut={onMouseOut}
    />
  );
}
