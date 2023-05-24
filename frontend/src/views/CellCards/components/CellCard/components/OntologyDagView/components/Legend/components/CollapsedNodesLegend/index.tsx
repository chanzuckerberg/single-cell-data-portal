import { LegendProps } from "../../../../common/types";
import { largeSize, backgroundColor } from "../../../../common/constants";
import { StyledLegendText } from "../common/style";

const CollapsedNodesLegend = ({ xPos, yPos }: LegendProps) => {
  return (
    <g>
      <rect
        x={xPos + 30 - largeSize + 2}
        y={yPos + 15 - largeSize + 1}
        fill={backgroundColor}
        stroke="#999999"
        strokeWidth={1}
        width={largeSize * 2 - 2}
        height={largeSize * 2 - 2}
      />
      <StyledLegendText x={xPos - 5} y={yPos}>
        Hidden terms
      </StyledLegendText>
    </g>
  );
};

export default CollapsedNodesLegend;
