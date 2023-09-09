import { backgroundColor } from "../../common/constants";
import InCorpusLegend from "./components/InCorpusLegend";
import DescendantsLegend from "./components/DescendantsLegend";
import HasMarkerGeneLegend from "./components/HasMarkerGeneLegend";

interface LegendProps {
  width: number;
  selectedGene: string | undefined;
}
export default function Legend({ width, selectedGene }: LegendProps) {
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
      {selectedGene ? (
        <HasMarkerGeneLegend
          selectedGene={selectedGene}
          xPos={width - 160}
          yPos={10}
        />
      ) : (
        <InCorpusLegend xPos={width - 160} yPos={10} />
      )}
      <DescendantsLegend xPos={width - 80} yPos={10} />
    </g>
  );
}
