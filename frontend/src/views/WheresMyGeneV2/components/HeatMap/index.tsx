import cloneDeep from "lodash/cloneDeep";
import {
  Dispatch,
  memo,
  SetStateAction,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";

import {
  EMPTY_ARRAY,
  EMPTY_OBJECT,
  EMPTY_SET,
} from "src/common/constants/utils";
import {
  CellTypeRow,
  generateTermsByKey,
  OntologyTerm,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import {
  DispatchContext,
  State,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  addCellInfoCellType,
  setFilteredCellTypes,
} from "src/views/WheresMyGene/common/store/actions";
import {
  CellType,
  ChartProps,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGene/common/types";
import YAxisChart from "./components/YAxisChart";
import { useTrackHeatMapLoaded } from "./hooks/useTrackHeatMapLoaded";

import {
  useSortedGeneNames,
  useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
} from "src/views/WheresMyGene/components/HeatMap/hooks/useSortedGeneNames";
import { useSortedCellTypesByTissueName } from "src/views/WheresMyGene/components/HeatMap/hooks/useSortedCellTypesByTissueName";
import {
  CellTypeFilterContainer,
  CellTypeTagContainer,
  ChartWrapper,
  Container,
  ContainerWrapper,
  StyledTag,
  XAxisMask,
  YAxisWrapper,
} from "src/views/WheresMyGene/components/HeatMap/style";
import { CellCountLabel } from "src/views/WheresMyGene/components/HeatMap/components/XAxisChart/style";
import {
  HEATMAP_CONTAINER_ID,
  MARGIN_BETWEEN_HEATMAPS,
} from "src/views/WheresMyGene/common/constants";
import Loader from "src/views/WheresMyGene/components/Loader";
import XAxisChart from "src/views/WheresMyGene/components/HeatMap/components/XAxisChart";
import Chart from "src/views/WheresMyGene/components/HeatMap/components/Chart";
import { hyphenize } from "src/views/WheresMyGene/components/HeatMap/utils";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";
import { Autocomplete, DefaultAutocompleteOption } from "@czi-sds/components";
import { Divider, TopLeftCornerMask, XAxisWrapper } from "./style";

interface Props {
  className?: string;
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  genes: State["selectedGenes"];
  selectedGeneExpressionSummariesByTissueName: {
    [tissueName: string]: GeneExpressionSummary[];
  };
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isLoadingAPI: boolean;
  isScaled: boolean;
  cellTypeSortBy: SORT_BY;
  geneSortBy: SORT_BY;
  echartsRendererMode: "svg" | "canvas";
  setAllChartProps: Dispatch<
    SetStateAction<{
      [tissue: string]: ChartProps;
    }>
  >;
  allChartProps: { [tissue: string]: ChartProps };
  tissuesByName: { [name: string]: OntologyTerm };
  setTissuesByName: Dispatch<
    SetStateAction<{
      [name: string]: OntologyTerm;
    }>
  >;
  expandedTissues: Set<string>;
  setExpandedTissues: Dispatch<SetStateAction<Set<string>>>;
}

export default memo(function HeatMap({
  className,
  cellTypes,
  genes,
  selectedGeneExpressionSummariesByTissueName,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isLoadingAPI,
  isScaled,
  cellTypeSortBy,
  geneSortBy,
  echartsRendererMode,
  allChartProps,
  setAllChartProps,
  tissuesByName,
  setTissuesByName,
  expandedTissues,
  setExpandedTissues,
}: Props): JSX.Element {
  const {
    xAxisHeight,
    selectedFilters: { tissues: filteredTissueIds },
    filteredCellTypes,
  } = useContext(StateContext);

  const selectedCellTypeOptions = useMemo(() => {
    return filteredCellTypes.map((cellType) => ({
      name: cellType,
    }));
  }, [filteredCellTypes]);

  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);
  const dispatch = useContext(DispatchContext);

  const { data } = usePrimaryFilterDimensions(2); //temp explicit version

  // Get tissueName to ID map for use in find marker genes
  useEffect(() => {
    let result: { [name: string]: OntologyTerm } = EMPTY_OBJECT;

    if (data) {
      const { tissues } = data;

      result = generateTermsByKey(tissues, "name");
    }

    setTissuesByName(result);
  }, [data, setTissuesByName]);

  const cellTypesByName = useMemo(() => {
    const result: { [name: string]: CellType } = {};

    Object.values(cellTypes).forEach((cellTypes) => {
      cellTypes.forEach((cellType) => {
        result[cellType.cellTypeName] = cellType;
      });
    });

    return result;
  }, [cellTypes]);

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

  const initialDisplayedCellTypeIds = useMemo(() => {
    return Object.values(tissuesByName).reduce((acc, { id }) => {
      if (
        (filteredTissueIds.length > 0 && filteredTissueIds.includes(id)) ||
        filteredTissueIds.length === 0
      ) {
        acc.add(id + id);
      }
      return acc;
    }, new Set<string>());
  }, [tissuesByName, filteredTissueIds]);

  // set of tissue names that are visible and set of cell types that are visible
  // presence is visible
  const [displayedCellTypes, setDisplayedCellTypes] = useState<Set<string>>(
    initialDisplayedCellTypeIds
  );

  useEffect(() => {
    setDisplayedCellTypes(initialDisplayedCellTypeIds);
  }, [initialDisplayedCellTypeIds]);

  const handleExpandCollapse = useCallback(
    (tissueID: string, tissueName: Tissue) => {
      const newDisplayedCellTypes = new Set<string>(displayedCellTypes);
      const newExpandedTissues = new Set<string>(expandedTissues);
      let addedTissue = false;

      if (expandedTissues.has(tissueID)) {
        newExpandedTissues.delete(tissueID);
      } else {
        newExpandedTissues.add(tissueID);
        addedTissue = true;
        track(EVENTS.WMG_TISSUE_EXPAND, { tissue: tissueName });
      }
      if (addedTissue) {
        sortedCellTypesByTissueName[tissueName].forEach((cellType) => {
          if (
            filteredCellTypes.length == 0 ||
            (filteredCellTypes.length > 0 &&
              filteredCellTypes.includes(cellType.cellTypeName))
          )
            newDisplayedCellTypes.add(tissueID + cellType.cellTypeName);
        });
      } else {
        [...newDisplayedCellTypes].forEach((cellType) => {
          if (
            cellType.includes(tissueID) &&
            cellType !== `${tissueID}${tissueID}`
          ) {
            newDisplayedCellTypes.delete(cellType);
          }
        });
      }

      setDisplayedCellTypes(newDisplayedCellTypes);
      setExpandedTissues(newExpandedTissues);
    },
    [
      displayedCellTypes,
      expandedTissues,
      setExpandedTissues,
      sortedCellTypesByTissueName,
      filteredCellTypes,
    ]
  );

  const uniqueCellTypes = useMemo(() => {
    const result: Set<string> = new Set<string>();
    Object.values(sortedCellTypesByTissueName).forEach((cellTypes) => {
      cellTypes.forEach((cellType) => {
        result.add(cellType.cellTypeName);
      });
    });
    return [...result].sort().map((cellType) => ({ name: cellType }));
  }, [sortedCellTypesByTissueName]);

  const handleFilteredCellTypesChange = (
    _: unknown,
    rawNewFilteredCellTypes: DefaultAutocompleteOption[]
  ) => {
    if (!dispatch) return;
    dispatch(
      setFilteredCellTypes(
        rawNewFilteredCellTypes.map((cellType) => cellType.name)
      )
    );
  };
  useEffect(() => {
    if (filteredCellTypes.length === 0) {
      setDisplayedCellTypes(initialDisplayedCellTypeIds);
      setExpandedTissues(EMPTY_SET as Set<string>);
    }
  }, [
    filteredCellTypes.length,
    initialDisplayedCellTypeIds,
    setExpandedTissues,
  ]);

  // Reset `displayedCellTypes` and `expandedTissues` when the user clears `filteredCellTypes`
  useEffect(() => {
    if (filteredCellTypes.length === 0) {
      // This is handled in the above useEffect, but we need to return early here so we don't do the work below
      return;
    }

    const newDisplayedCellTypes = new Set<string>();
    const newExpandedTissues = new Set<string>();

    Object.entries(sortedCellTypesByTissueName).forEach(
      ([tissue, cellTypes]) => {
        if (
          filteredTissueIds.length > 0 &&
          !filteredTissueIds.includes(tissuesByName[tissue].id)
        ) {
          return;
        }

        cellTypes.forEach((cellType) => {
          if (filteredCellTypes.includes(cellType.name)) {
            newDisplayedCellTypes.add(
              tissuesByName[tissue].id + tissuesByName[tissue].id
            );
            newDisplayedCellTypes.add(tissuesByName[tissue].id + cellType.name);
            newExpandedTissues.add(tissuesByName[tissue].id);
          }
        });
      }
    );
    const filteredCellTypeIds = filteredCellTypes.map(
      (cellType) => cellTypesByName[cellType].id
    );
    track(EVENTS.WMG_SELECT_CELL_TYPE, {
      cell_types: filteredCellTypeIds,
    });

    setDisplayedCellTypes(newDisplayedCellTypes);
    setExpandedTissues(newExpandedTissues);
  }, [
    cellTypesByName,
    filteredCellTypes,
    filteredTissueIds,
    initialDisplayedCellTypeIds,
    setExpandedTissues,
    sortedCellTypesByTissueName,
    tissuesByName,
  ]);

  const handleCellTypeDelete = (cellTypeToDelete: string) => () => {
    if (!dispatch) return;
    const newValue = filteredCellTypes.filter(
      (cellType) => !(cellTypeToDelete === cellType)
    );
    dispatch(setFilteredCellTypes(newValue));
  };

  useTrackHeatMapLoaded({
    selectedGenes: genes,
    displayedCellTypes,
    selectedCellTypes: filteredCellTypes,
  });

  /**
   * All tissue cell types to render in YAxisCharts
   */
  const allTissueCellTypes = useMemo(() => {
    return Object.values(tissuesByName)
      .sort((a, b) => {
        // sort tissues alphabetically
        return a.name.localeCompare(b.name);
      })
      .map((tissue: OntologyTerm) => {
        const { id, name } = tissue;

        return {
          tissueId: id,
          tissueName: name,
          tissueCellTypes: getTissueCellTypes({
            cellTypeSortBy,
            cellTypes,
            sortedCellTypesByTissueName,
            tissue: name,
            tissueID: id,
            displayedCellTypes,
          }),
        };
      });
  }, [
    tissuesByName,
    cellTypeSortBy,
    cellTypes,
    sortedCellTypesByTissueName,
    displayedCellTypes,
  ]);

  return (
    <>
      <ContainerWrapper>
        <TopLeftCornerMask height={xAxisHeight}>
          <CellTypeFilterContainer id="celltype-filter-container">
            <Autocomplete
              className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
              multiple
              label="Search cell types"
              value={selectedCellTypeOptions}
              search
              onChange={handleFilteredCellTypesChange}
              options={uniqueCellTypes}
            />
            <CellTypeTagContainer>
              {filteredCellTypes.map((cellType) => (
                <StyledTag
                  label={cellType}
                  key={cellType}
                  onDelete={handleCellTypeDelete(cellType)}
                />
              ))}
            </CellTypeTagContainer>
          </CellTypeFilterContainer>
          <CellCountLabel>Cell Count</CellCountLabel>
        </TopLeftCornerMask>
        <Container {...{ className }} id={HEATMAP_CONTAINER_ID}>
          {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}
          <XAxisWrapper id="x-axis-wrapper">
            <XAxisMask data-testid="x-axis-mask" height={xAxisHeight} />
            <XAxisChart geneNames={sortedGeneNames} />
          </XAxisWrapper>
          <YAxisWrapper top={0}>
            {allTissueCellTypes.map(
              ({ tissueId, tissueName, tissueCellTypes }) => {
                return (
                  tissueCellTypes.length > 0 && (
                    <YAxisChart
                      key={tissueName}
                      tissue={tissueName}
                      tissueID={tissueId}
                      cellTypes={tissueCellTypes}
                      generateMarkerGenes={generateMarkerGenes}
                      expandedTissues={expandedTissues}
                      handleExpandCollapse={handleExpandCollapse}
                    />
                  )
                );
              }
            )}
          </YAxisWrapper>
          <ChartWrapper ref={chartWrapperRef} top={xAxisHeight}>
            {allTissueCellTypes.map(({ tissueName, tissueCellTypes }) => {
              // Don't render anything if tissue has no cell types for some reason
              if (!tissueCellTypes.length) return;

              const selectedGeneData =
                orderedSelectedGeneExpressionSummariesByTissueName[tissueName];

              /**
               * (thuang): If there is no selected gene data, we don't want to render
               * the chart, because it will cause the chart to render with 0 width,
               * which is an error for echarts
               */
              if (!selectedGeneData?.length) {
                const height =
                  document.getElementById(`${hyphenize(tissueName)}-y-axis`)
                    ?.clientHeight ?? 0;
                return (
                  <div
                    id={`no-chart-data-${hyphenize(tissueName)}`} // Not used, just to make it stand out
                    key={`${tissueName}-${echartsRendererMode}`}
                    style={{
                      height: `${height + MARGIN_BETWEEN_HEATMAPS}px`,
                    }}
                  />
                );
              }

              return (
                <Chart
                  isScaled={isScaled}
                  /**
                   * (thuang): We use `key` to force re-render the HeatMap component
                   * when the renderer mode changes, so echarts can create new instances
                   */
                  key={`${tissueName}-${echartsRendererMode}`}
                  tissue={tissueName}
                  cellTypes={tissueCellTypes}
                  selectedGeneData={
                    orderedSelectedGeneExpressionSummariesByTissueName[
                      tissueName
                    ]
                  }
                  setIsLoading={setIsLoading}
                  scaledMeanExpressionMax={scaledMeanExpressionMax}
                  scaledMeanExpressionMin={scaledMeanExpressionMin}
                  echartsRendererMode={echartsRendererMode}
                  setAllChartProps={setAllChartProps}
                  chartProps={allChartProps[tissueName]}
                />
              );
            })}
          </ChartWrapper>
        </Container>
        <Divider className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME} />
      </ContainerWrapper>
    </>
  );
});

function getTissueCellTypes({
  cellTypes,
  sortedCellTypesByTissueName,
  tissue,
  tissueID,
  cellTypeSortBy,
  displayedCellTypes,
}: {
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellTypeRow[] };
  tissue: Tissue;
  tissueID: string;
  cellTypeSortBy: SORT_BY;
  displayedCellTypes: Set<string>;
}) {
  const tissueCellTypes = cellTypes[tissue];

  if (!tissueCellTypes || tissueCellTypes.length === 0) return EMPTY_ARRAY;

  const sortedTissueCellTypes = sortedCellTypesByTissueName[tissue];

  let ret =
    (cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? tissueCellTypes
      : sortedTissueCellTypes) || EMPTY_ARRAY;

  ret = ret.filter((cellType) =>
    displayedCellTypes.has(tissueID + cellType.cellTypeName)
  );

  return ret;
}

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
