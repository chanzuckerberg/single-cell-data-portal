import { LegendProps } from "../../../../common/types";
import { StyledLegendText } from "../common/style";
import { largeSize } from "../../../../common/constants";

const DescendantsLegend = ({ xPos, yPos }: LegendProps) => {
  return (
    <g>
      <defs>
        <pattern
          id="crosshatch"
          width="5"
          height="5"
          patternUnits="userSpaceOnUse"
        >
          <path d="M 5 0 L 0 5" stroke="#999999" strokeWidth="0.5" />
        </pattern>
      </defs>
      <StyledLegendText x={xPos + 2.5} y={yPos + 30}>
        Yes
      </StyledLegendText>
      <circle
        cx={xPos + 12.5}
        cy={yPos + 15}
        fill="url(#crosshatch)"
        stroke="#999999"
        strokeWidth={1}
        r={largeSize}
      />
      <StyledLegendText x={xPos - 5} y={yPos}>
        Descendants
      </StyledLegendText>
      <StyledLegendText x={xPos + 33.5} y={yPos + 30}>
        No
      </StyledLegendText>
      <rect
        x={xPos + 40 - largeSize + 2}
        y={yPos + 15 - largeSize + 1}
        fill="url(#crosshatch)"
        stroke="#999999"
        strokeWidth={1}
        width={largeSize * 2 - 2}
        height={largeSize * 2 - 2}
      />
    </g>
  );
};

export default DescendantsLegend;
