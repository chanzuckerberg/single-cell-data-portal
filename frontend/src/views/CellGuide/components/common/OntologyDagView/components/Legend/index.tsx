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
import { StyledPie } from "../../common/style";
import { CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID } from "./constants";
import {
  largeSize,
  maxMarkerScore,
  nodePieChartCircleScaler,
} from "../../common/constants";
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
  const numPies = 6;
  const numMarkerScores = maxMarkerScore + 1;

  const fill = "#999999";
  const expressedInCellsPercLegendComponent = (
    <LegendItemWrapper>
      Expressed in Cells(%)
      <MarkerScoreWrapper>
        {Array.from({ length: numPies }).map((_, i) => (
          <FlexColumn key={i}>
            <ExpressedInCells
              degree={(360 / (numPies - 1)) * i}
              fill={fill}
              size={size}
            />
            {(i == 0 || i == numPies - 1) && (
              <CenterText>{(100 / (numPies - 1)) * i}</CenterText>
            )}
          </FlexColumn>
        ))}
      </MarkerScoreWrapper>
    </LegendItemWrapper>
  );

  const markerScoreLegendComponent = (
    <LegendItemWrapper>
      Marker Score ({selectedGene})
      <MarkerScoreWrapper>
        {Array.from({ length: numMarkerScores }).map((_, i) => (
          <FlexColumn key={i}>
            <MarkerScore fill={i} size={size} />
            {(i == 0 || i == numMarkerScores - 1) && (
              <CenterText>{i}</CenterText>
            )}
          </FlexColumn>
        ))}
      </MarkerScoreWrapper>
    </LegendItemWrapper>
  );
  const hasMarkerGeneLegend = (
    <>
      {markerScoreLegendComponent}
      {expressedInCellsPercLegendComponent}
    </>
  );
  return (
    <LegendWrapper data-testid={CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID}>
      {!selectedGene && targetNodeLegendComponent}
      {selectedGene ? hasMarkerGeneLegend : corpusLegendComponent}
      {!selectedGene && descendantsLegendComponent}
    </LegendWrapper>
  );
}

interface ExpressedInCellsProps {
  degree: number;
  fill: string;
  size: number;
}
const ExpressedInCells = ({ degree, fill, size }: ExpressedInCellsProps) => (
  <StyledPie degree={degree} fill={fill} size={size}>
    <StyledPie
      degree={360}
      fill={fill}
      size={size * nodePieChartCircleScaler}
      opacity={0.5}
      center
    />
  </StyledPie>
);

interface MarkerScoreProps {
  fill: number;
  size: number;
}
const MarkerScore = ({ fill, size }: MarkerScoreProps) => (
  <StyledPie degree={360} fill={fill} size={size} />
);
