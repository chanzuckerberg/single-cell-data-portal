import { backgroundColor } from "../../common/constants";
import InCorpusLegend from "./components/InCorpusLegend";
import DescendantsLegend from "./components/DescendantsLegend";
import CollapsedNodesLegend from "./components/CollapsedNodesLegend";

interface LegendProps {
  width: number;
}
export default function Legend({ width }: LegendProps) {
  return (
    <g>
      <rect
        x={width - 260}
        y={-10}
        width={260}
        height={60}
        fill={backgroundColor}
        rx={4}
      />
      <InCorpusLegend xPos={width - 240} yPos={10} />
      <DescendantsLegend xPos={width - 160} yPos={10} />
      <CollapsedNodesLegend xPos={width - 80} yPos={10} />
    </g>
  );
}
