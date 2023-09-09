import { backgroundColor } from "../../common/constants";
import InCorpusLegend from "./components/InCorpusLegend";
import DescendantsLegend from "./components/DescendantsLegend";

interface LegendProps {
  width: number;
}
export default function Legend({ width }: LegendProps) {
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
      <InCorpusLegend xPos={width - 160} yPos={10} />
      <DescendantsLegend xPos={width - 80} yPos={10} />
    </g>
  );
}
