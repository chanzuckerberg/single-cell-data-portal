import { LegendProps } from "../../../../common/types";
import { StyledLegendText } from "../common/style";
import {
  primaryColor,
  largeSize,
  tertiaryColor,
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
        fill={primaryColor}
        r={largeSize}
      />
      <StyledLegendText x={xPos} y={yPos}>
        In Corpus
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
