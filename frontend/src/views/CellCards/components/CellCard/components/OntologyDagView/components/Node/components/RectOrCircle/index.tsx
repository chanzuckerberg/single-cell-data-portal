import { MouseEventHandler } from "react";
import { TreeNodeWithState } from "../../../../common/types";
import {
  tertiaryColor,
  primaryColor,
  highlightColor,
  smallSize,
  largeSize,
} from "../../../../common/constants";
import { StyledRect, StyledCircle } from "./style";

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
}

export default function RectOrCircle({
  node,
  handleClick,
  isTargetNode,
  handleMouseOver,
  handleMouseOut,
  animationKey,
}: RectOrCircleProps) {
  let color = tertiaryColor;
  if (node.n_cells > 0) {
    color = primaryColor;
  }
  if (isTargetNode) {
    color = highlightColor;
  }
  const size = node.n_cells === 0 ? smallSize : largeSize;
  if (node.id.startsWith("dummy-child")) {
    color = primaryColor;
  }

  const onMouseOver = node.id.startsWith("dummy-child")
    ? undefined
    : (event: React.MouseEvent<SVGElement>) => {
        handleMouseOver(event, node);
      };
  const onMouseOut = node.name.startsWith("dummy-child")
    ? undefined
    : handleMouseOut;
  return node?.children?.length || node.id.startsWith("dummy-child") ? (
    <StyledCircle
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
