import cloneDeep from "lodash/cloneDeep";
import { memo, useContext, useMemo, useRef, useState } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  generateTermsByKey,
  OntologyTerm,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { HEATMAP_CONTAINER_ID } from "../../common/constants";
import { DispatchContext, State } from "../../common/store";
import { addCellInfoCellType } from "../../common/store/actions";
import {
  CellType,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "../../common/types";
import Loader from "../Loader";
import Chart from "./components/Chart";
import XAxisChart from "./components/XAxisChart";
import { CellCountLabel } from "./components/XAxisChart/style";
import YAxisChart from "./components/YAxisChart";
import { useSortedCellTypesByTissueName } from "./hooks/useSortedCellTypesByTissueName";
import {
  useSortedGeneNames,
  useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
} from "./hooks/useSortedGeneNames";
import { useTrackHeatMapLoaded } from "./hooks/useTrackHeatMapLoaded";
import {
  ChartWrapper,
  Container,
  ContainerWrapper,
  TopLeftCornerMask,
  YAxisWrapper,
  XAxisMask,
  XAxisWrapper,
} from "./style";

interface Props {
  className?: string;
  selectedTissues: string[];
  cellTypes: { [tissue: Tissue]: CellType[] };
  genes: State["selectedGenes"];
  tissuesWithDeletedCellTypes: string[];
  allTissueCellTypes: { [tissue: Tissue]: CellType[] };
  selectedGeneExpressionSummariesByTissueName: {
    [tissueName: string]: GeneExpressionSummary[];
  };
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isLoadingAPI: boolean;
  isScaled: boolean;
  cellTypeSortBy: SORT_BY;
  geneSortBy: SORT_BY;
  selectedOrganismId: string;
}

export default memo(function HeatMap({
  className,
  selectedTissues,
  cellTypes,
  genes,
  tissuesWithDeletedCellTypes,
  allTissueCellTypes,
  selectedGeneExpressionSummariesByTissueName,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isLoadingAPI,
  isScaled,
  cellTypeSortBy,
  geneSortBy,
  selectedOrganismId,
}: Props): JSX.Element {
  useTrackHeatMapLoaded({ selectedGenes: genes, selectedTissues });

  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);

  const dispatch = useContext(DispatchContext);

  const { data } = usePrimaryFilterDimensions();
  // Get tissueName to ID map for use in find marker genes
  const tissuesByName = useMemo(() => {
    let result: { [name: string]: OntologyTerm } = {};

    if (!data) return result;

    const { tissues } = data;

    result = generateTermsByKey(tissues, "name");

    return result;
  }, [data]);

  const generateMarkerGenes = (cellType: CellType, tissueID: string) => {
    if (!dispatch) return;
    dispatch(addCellInfoCellType({ cellType, tissueID }));
  };

  const tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap =
    useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
      selectedGeneExpressionSummariesByTissueName
    );

  const sortedGeneNames = useSortedGeneNames({
    geneSortBy,
    genes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const sortedCellTypesByTissueName = useSortedCellTypesByTissueName({
    cellTypeSortBy,
    genes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const geneNameToIndex = useMemo(() => {
    const result: { [key: string]: number } = {};

    for (const [index, gene] of Object.entries(sortedGeneNames)) {
      result[gene] = Number(index);
    }

    return result;
  }, [sortedGeneNames]);

  const orderedSelectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: { [tissueName: string]: GeneExpressionSummary[] } = {};

    for (const [tissueName, geneExpressionSummary] of Object.entries(
      selectedGeneExpressionSummariesByTissueName
    )) {
      // (thuang): sort() mutates the array, so we need to clone it
      result[tissueName] = cloneDeep(
        geneExpressionSummary.sort((a, b) => {
          if (!a || !b) return -1;

          return geneNameToIndex[a.name] - geneNameToIndex[b.name];
        })
      );
    }

    return result;
  }, [selectedGeneExpressionSummariesByTissueName, geneNameToIndex]);

  return (
    <ContainerWrapper>
      <TopLeftCornerMask>
        <CellCountLabel>Cell Count</CellCountLabel>
      </TopLeftCornerMask>
      <Container {...{ className }} id={HEATMAP_CONTAINER_ID}>
        {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}
        <XAxisWrapper>
          <XAxisMask />
          <XAxisChart geneNames={sortedGeneNames} />
        </XAxisWrapper>
        <YAxisWrapper>
          {selectedTissues.map((tissue) => {
            const tissueCellTypes = getTissueCellTypes({
              cellTypeSortBy,
              cellTypes,
              sortedCellTypesByTissueName,
              tissue,
            });

            return (
              <YAxisChart
                key={tissue}
                tissue={tissue}
                tissueID={tissuesByName[tissue].id}
                cellTypes={tissueCellTypes}
                hasDeletedCellTypes={tissuesWithDeletedCellTypes.includes(
                  tissue
                )}
                availableCellTypes={allTissueCellTypes[tissue]}
                generateMarkerGenes={generateMarkerGenes}
                selectedOrganismId={selectedOrganismId}
              />
            );
          })}
        </YAxisWrapper>
        <ChartWrapper ref={chartWrapperRef}>
          {selectedTissues.map((tissue) => {
            const tissueCellTypes = getTissueCellTypes({
              cellTypeSortBy,
              cellTypes,
              sortedCellTypesByTissueName,
              tissue,
            });

            return (
              <Chart
                isScaled={isScaled}
                key={tissue}
                tissue={tissue}
                cellTypes={tissueCellTypes}
                selectedGeneData={
                  orderedSelectedGeneExpressionSummariesByTissueName[tissue]
                }
                setIsLoading={setIsLoading}
                scaledMeanExpressionMax={scaledMeanExpressionMax}
                scaledMeanExpressionMin={scaledMeanExpressionMin}
              />
            );
          })}
        </ChartWrapper>
      </Container>
    </ContainerWrapper>
  );
});

function getTissueCellTypes({
  cellTypes,
  sortedCellTypesByTissueName,
  tissue,
  cellTypeSortBy,
}: {
  cellTypes: { [tissue: Tissue]: CellType[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellType[] };
  tissue: Tissue;
  cellTypeSortBy: SORT_BY;
}) {
  const tissueCellTypes = cellTypes[tissue];
  const sortedTissueCellTypes = sortedCellTypesByTissueName[tissue];

  return (
    (cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? tissueCellTypes
      : sortedTissueCellTypes) || EMPTY_ARRAY
  );
}

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
