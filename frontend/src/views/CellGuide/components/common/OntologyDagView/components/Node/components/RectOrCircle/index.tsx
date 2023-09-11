import { MouseEventHandler } from "react";
import { MarkerGeneStats, TreeNodeWithState } from "../../../../common/types";
import {
  tertiaryColor,
  primaryColor,
  highlightColor,
  smallSize,
  largeSize,
  leftBorderSpacing,
  nodePieChartCircleScaler,
} from "../../../../common/constants";
import { StyledPie } from "../../../../common/style";
import { CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID } from "src/views/CellGuide/components/common/OntologyDagView/components/Node/components/RectOrCircle/constants";
import { Border, NodeWrapper } from "./style";

const DUMMY_CHILD = "dummy-child";

interface RectOrCircleProps {
  handleClick?: MouseEventHandler<HTMLDivElement>;
  isTargetNode: boolean;
  node: TreeNodeWithState;
  handleMouseOver: (
    event: React.MouseEvent<HTMLDivElement>,
    datum: TreeNodeWithState
  ) => void;
  handleMouseOut: () => void;
  animationKey: string;
  cellTypeMarkerGeneStats: MarkerGeneStats | undefined;
  inMarkerGeneMode: boolean;
}

const getColor = (
  node: TreeNodeWithState,
  isTargetNode: boolean,
  inMarkerGeneMode: boolean
) => {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
  if (node.id.startsWith(DUMMY_CHILD)) {
    color = primaryColor;
  }

  if (isTargetNode) {
    color = highlightColor;
  }
  if (inMarkerGeneMode) {
    color = tertiaryColor;
  }
  return color;
};

const getMouseOverHandler = (
  node: TreeNodeWithState,
  handleMouseOver: RectOrCircleProps["handleMouseOver"]
) => {
  return node.id.startsWith(DUMMY_CHILD)
    ? undefined
    : (event: React.MouseEvent<HTMLDivElement>) => {
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
  cellTypeMarkerGeneStats,
  inMarkerGeneMode,
}: RectOrCircleProps) {
  const color = getColor(node, isTargetNode, inMarkerGeneMode);
  const size = (node.id.startsWith(DUMMY_CHILD) ? smallSize : largeSize) * 2;
  const onMouseOver = getMouseOverHandler(node, handleMouseOver);
  const onMouseOut = getMouseOutHandler(node, handleMouseOut);
  const sizeScaler = cellTypeMarkerGeneStats ? nodePieChartCircleScaler : 1;
  const word = node?.children?.length ? "has" : "no";
  const dataTestId = `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${node.id}-${word}-children-isTargetNode=${isTargetNode}`;
  const hasDescendants =
    node?.children?.length || node.id.startsWith(DUMMY_CHILD);

  const xOffset = hasDescendants ? 0 : -leftBorderSpacing * 2;
  return (
    <foreignObject
      data-testid={dataTestId}
      x={-size / 2 + xOffset}
      y={-size / 2}
      width={size - xOffset}
      height={size}
      key={animationKey}
    >
      <div
        onMouseLeave={onMouseOut}
        onClick={handleClick}
        onMouseOver={onMouseOver}
      >
        <NodeWrapper columnGap={leftBorderSpacing}>
          {!hasDescendants && (
            <Border width={leftBorderSpacing} height={size} />
          )}
          {cellTypeMarkerGeneStats ? (
            <StyledPie
              degree={cellTypeMarkerGeneStats.pc * 360}
              size={size}
              fill={cellTypeMarkerGeneStats.marker_score}
            >
              <StyledPie
                size={size * sizeScaler}
                fill={
                  cellTypeMarkerGeneStats
                    ? cellTypeMarkerGeneStats.marker_score
                    : color
                }
                opacity={0.5}
                degree={360}
                center
              />
            </StyledPie>
          ) : (
            <StyledPie size={size * sizeScaler} fill={color} degree={360} />
          )}
        </NodeWrapper>
      </div>
    </foreignObject>
  );
}
