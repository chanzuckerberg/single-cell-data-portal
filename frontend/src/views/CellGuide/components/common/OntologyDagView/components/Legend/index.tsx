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
  StyledBorder,
  StyledNodeWrapper,
} from "./style";
import { StyledPie } from "../../common/style";
import { CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID } from "./constants";
import {
  largeSize,
  maxMarkerScore,
  nodePieChartCircleScaler,
  secondaryColor,
} from "../../common/constants";
interface LegendProps {
  selectedGene: string | undefined;
  isTissue?: boolean;
}
export default function Legend({ selectedGene, isTissue }: LegendProps) {
  const targetNodeLegendComponent = (
    <LegendItemWrapper>
      Cell Type
      <IsTargetNode />
      Selected
    </LegendItemWrapper>
  );
  const size = largeSize * 2;
  const descendantsLegendComponent = (
    <LegendItemWrapper>
      <FlexColumn>
        <StyledNodeWrapper>
          <NoDescendantsLine size={size} />
          <StyledBorder width={1} height={size * 2} />
        </StyledNodeWrapper>
        No Descendants
      </FlexColumn>
    </LegendItemWrapper>
  );
  const corpusLegendComponent = (
    <>
      <LegendItemWrapper>
        <LegendItem>
          <FlexColumn>
            <IsInCorpus />
            <CenterText>In Census</CenterText>
          </FlexColumn>
        </LegendItem>
      </LegendItemWrapper>
      <LegendItemWrapper>
        <LegendItem>
          <FlexColumn>
            <IsNotInCorpus />
            <CenterText>Not In Census</CenterText>
          </FlexColumn>
        </LegendItem>
      </LegendItemWrapper>
    </>
  );

  const numPies = 6;
  const numMarkerScores = maxMarkerScore + 1;

  const expressedInCellsPercLegendComponent = (
    <LegendItemWrapper>
      Expressed in Cells(%)
      <MarkerScoreWrapper>
        {Array.from({ length: numPies }).map((_, i) => (
          <FlexColumn key={i}>
            <ExpressedInCells
              degree={(360 / (numPies - 1)) * i}
              fill={secondaryColor}
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
  <StyledPie center degree={degree} fill={fill} size={size}>
    <StyledPie
      degree={360}
      fill={fill}
      size={size * nodePieChartCircleScaler}
      opacity={0.5}
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
