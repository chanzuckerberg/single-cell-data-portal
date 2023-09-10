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
  NoMarkerGene,
  YesMarkergene,
} from "./style";

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
  const hasMarkerGeneLegend = (
    <LegendItemWrapper>
      {selectedGene} Is Marker
      <LegendItem>
        <FlexColumn>
          <YesMarkergene />
          <CenterText>Yes</CenterText>
        </FlexColumn>
        <FlexColumn>
          <NoMarkerGene />
          <CenterText>No</CenterText>
        </FlexColumn>
      </LegendItem>
    </LegendItemWrapper>
  );

  return (
    <LegendWrapper>
      {targetNodeLegendComponent}
      {selectedGene ? hasMarkerGeneLegend : corpusLegendComponent}
      {descendantsLegendComponent}
    </LegendWrapper>
  );
}
