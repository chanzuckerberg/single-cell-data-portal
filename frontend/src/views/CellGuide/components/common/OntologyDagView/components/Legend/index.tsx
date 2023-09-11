import {
  CenterText,
  FlexColumn,
  IsInCorpus,
  IsNotInCorpus,
  IsTargetNode,
  LegendItem,
  LegendItemWrapper,
  LegendWrapper,
  MarkerScoreWrapper,
  NoDescendantsLine,
} from "./style";
import { StyledPie } from "../../common/style";
import { CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID } from "./constants";
import {
  largeSize,
  leftBorderSpacing,
  maxMarkerScore,
  nodePieChartCircleScaler,
} from "../../common/constants";
import { Border, NodeWrapper } from "../Node/components/RectOrCircle/style";
interface LegendProps {
  selectedGene: string | undefined;
  isTissue?: boolean;
}
export default function Legend({ selectedGene, isTissue }: LegendProps) {
  const targetNodeLegendComponent = (
    <LegendItemWrapper>
      Current
      <IsTargetNode />
    </LegendItemWrapper>
  );
  const size = largeSize * 2;
  const descendantsLegendComponent = (
    <LegendItemWrapper>
      <FlexColumn>
        No Descendents
        <NodeWrapper columnGap={leftBorderSpacing}>
          <NoDescendantsLine size={size} />
          <Border width={leftBorderSpacing} height={size} />
          <IsNotInCorpus />
        </NodeWrapper>
      </FlexColumn>
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
      Marker Score
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
      {!selectedGene && !isTissue && targetNodeLegendComponent}
      {selectedGene ? hasMarkerGeneLegend : corpusLegendComponent}
      {descendantsLegendComponent}
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
