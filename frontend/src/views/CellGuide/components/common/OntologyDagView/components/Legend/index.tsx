import { backgroundColor } from "../../common/constants";
import InCorpusLegend from "./components/InCorpusLegend";
import DescendantsLegend from "./components/DescendantsLegend";
import HasMarkerGeneLegend from "./components/HasMarkerGeneLegend";

interface LegendProps {
  width: number;
  markerGeneMode: boolean;
}
export default function Legend({ width, markerGeneMode }: LegendProps) {
  return (
    <g>
      <rect
        x={width - 180}
        y={-10}
        width={180}
        height={60}
        fill={backgroundColor}
        rx={4}
      />
      {markerGeneMode ? (
        <HasMarkerGeneLegend xPos={width - 160} yPos={10} />
      ) : (
        <InCorpusLegend xPos={width - 160} yPos={10} />
      )}
      <DescendantsLegend xPos={width - 80} yPos={10} />
    </g>
  );
}
