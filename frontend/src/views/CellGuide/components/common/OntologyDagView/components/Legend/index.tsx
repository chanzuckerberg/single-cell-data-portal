import {
  CenterText,
  FlexColumn,
  HasDescendants,
  HasNoDescendants,
  IsInCorpus,
  IsNotInCorpus,
  IsTargetNode,
  LegendItem,
  LegendItemWrapper,
  LegendWrapper,
  MarkerScoreWrapper,
} from "./style";
import { StyledPie, StyledCircle } from "../../common/style";
import { CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID } from "./constants";
import { largeSize, nodePieChartCircleScaler } from "../../common/constants";
interface LegendProps {
  selectedGene: string | undefined;
}
export default function Legend({ selectedGene }: LegendProps) {
  const targetNodeLegendComponent = (
    <LegendItemWrapper>
      Current
      <IsTargetNode />
    </LegendItemWrapper>
  );
  const descendantsLegendComponent = (
    <LegendItemWrapper>
      Descendants
      <LegendItem>
        <FlexColumn>
          <HasDescendants />
          <CenterText>Yes</CenterText>
        </FlexColumn>
        <FlexColumn>
          <HasNoDescendants />
          <CenterText>No</CenterText>
        </FlexColumn>
      </LegendItem>
    </LegendItemWrapper>
  );
  const corpusLegendComponent = (
    <LegendItemWrapper>
      In Corpus
      <LegendItem>
        <FlexColumn>
          <IsInCorpus />
          <CenterText>Yes</CenterText>
        </FlexColumn>
        <FlexColumn>
          <IsNotInCorpus />
          <CenterText>No</CenterText>
        </FlexColumn>
      </LegendItem>
    </LegendItemWrapper>
  );
  const size = largeSize * 2;
  const numPies = 5;
  const fill = "#999999";
  const hasMarkerGeneLegend = (
    <LegendItemWrapper>
      Marker Score({selectedGene})
      <MarkerScoreWrapper>
        {Array.from({ length: numPies }).map((_, i) => (
          <MarkerScoreLegendItem
            key={i}
            degree={(360 / numPies) * i}
            fill={fill}
            size={size}
          />
        ))}
      </MarkerScoreWrapper>
    </LegendItemWrapper>
  );

  return (
    <LegendWrapper data-testid={CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID}>
      {!selectedGene && targetNodeLegendComponent}
      {selectedGene ? hasMarkerGeneLegend : corpusLegendComponent}
      {!selectedGene && descendantsLegendComponent}
    </LegendWrapper>
  );
}

interface MarkerScoreLegendItemProps {
  degree: number;
  fill: string | number;
  size: number;
}
const MarkerScoreLegendItem = ({
  degree,
  fill,
  size,
}: MarkerScoreLegendItemProps) => (
  <StyledPie degree={degree} fill={fill} size={size}>
    <StyledCircle fill={`${fill}80`} size={size * nodePieChartCircleScaler} />
  </StyledPie>
);
