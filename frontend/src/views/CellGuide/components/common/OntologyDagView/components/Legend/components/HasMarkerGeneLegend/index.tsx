import { LegendProps } from "../../../../common/types";
import { StyledLegendText } from "../common/style";
import {
  largeSize,
  tertiaryColor,
  markerGeneModeColor,
} from "../../../../common/constants";

interface HasMarkerGeneLegendProps extends LegendProps {
  selectedGene: string;
}
const HasMarkerGeneLegend = ({
  xPos,
  yPos,
  selectedGene,
}: HasMarkerGeneLegendProps) => {
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
      <StyledLegendText x={xPos} y={yPos}>
        {selectedGene}
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

export default HasMarkerGeneLegend;
