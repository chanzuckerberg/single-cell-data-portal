import { LegendProps } from "../../../../common/types";
import { StyledLegendText } from "../common/style";
import {
  largeSize,
  tertiaryColor,
  markerGeneModeColor,
} from "../../../../common/constants";

const InCorpusLegend = ({ xPos, yPos }: LegendProps) => {
  return (
    <g>
      <StyledLegendText x={xPos + 2.5} y={yPos + 30}>
        Yes
      </StyledLegendText>
      <circle
        cx={xPos + 12.5}
        cy={yPos + 15}
        fill={markerGeneModeColor}
        r={largeSize}
      />
      <StyledLegendText x={xPos - 20} y={yPos}>
        Has Marker Gene
      </StyledLegendText>
      <StyledLegendText x={xPos + 33} y={yPos + 30}>
        No
      </StyledLegendText>
      <circle
        cx={xPos + 40}
        cy={yPos + 15}
        fill={tertiaryColor}
        r={largeSize}
      />
    </g>
  );
};

export default InCorpusLegend;
