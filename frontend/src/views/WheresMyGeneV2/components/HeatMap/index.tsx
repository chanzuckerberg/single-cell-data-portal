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

import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
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
} from "src/views/WheresMyGeneV2/common/store";
import {
  addCellInfoCellType,
  setFilteredCellTypes,
  toggleExpandedTissueId,
} from "src/views/WheresMyGeneV2/common/store/actions";
import {
  CellType,
  ChartProps,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGeneV2/common/types";
import YAxisChart from "./components/YAxisChart";

import {
  useSortedGeneNames,
  useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
} from "src/views/WheresMyGene/components/HeatMap/hooks/useSortedGeneNames";
import { useSortedCellTypesByTissueName } from "src/views/WheresMyGene/components/HeatMap/hooks/useSortedCellTypesByTissueName";
import {
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
} from "src/views/WheresMyGeneV2/common/constants";
import Loader from "src/views/WheresMyGene/components/Loader";
import XAxisChart from "src/views/WheresMyGene/components/HeatMap/components/XAxisChart";
import Chart from "src/views/WheresMyGeneV2/components/HeatMap/components/Chart";
import { hyphenize } from "src/views/WheresMyGene/components/HeatMap/utils";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";
import { Autocomplete, DefaultAutocompleteOption } from "@czi-sds/components";
import {
  CellTypeFilterContainer,
  Divider,
  TopLeftCornerMask,
  XAxisWrapper,
} from "./style";
import {
  useHandleExpandedTissueIds,
  useTrackHeatMapLoaded,
} from "src/views/WheresMyGeneV2/components/HeatMap/hooks";

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
  sidebarWidth: number;
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
  sidebarWidth,
}: Props): JSX.Element {
  const {
    xAxisHeight,
    selectedFilters: { tissues: filteredTissueIds },
    filteredCellTypes,
    filteredCellTypeIds,
    expandedTissueIds,
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
    dispatch?.(addCellInfoCellType({ cellType, tissueID }));
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

  /**
   * (thuang): Tissues to display after applying filters
   */
  const displayedTissues = useMemo(() => {
    return Object.values(tissuesByName)
      .filter(({ id }) => {
        return !filteredTissueIds.length || filteredTissueIds.includes(id);
      })
      .filter(({ name }) => {
        if (!filteredCellTypes.length) return true;

        const tissueCellTypes = sortedCellTypesByTissueName[name];

        return tissueCellTypes?.some((cellType) => {
          return filteredCellTypes.includes(cellType.cellTypeName);
        });
      });
  }, [
    filteredTissueIds,
    filteredCellTypes,
    sortedCellTypesByTissueName,
    tissuesByName,
  ]);

  const displayedTissueIds = useMemo(() => {
    return displayedTissues.map(({ id }) => id);
  }, [displayedTissues]);

  /**
   * (thuang): Derive displayed cell types from `displayedTissues`,
   * `expandedTissueIds`, and `filteredCellTypes`
   */
  const displayedCellTypes = useMemo(() => {
    const result = new Set<string>();

    displayedTissues.forEach(({ id, name }) => {
      result.add(id + id);

      if (expandedTissueIds.includes(id)) {
        const tissueCellTypes = sortedCellTypesByTissueName[name];

        tissueCellTypes?.forEach((cellType) => {
          if (
            !filteredCellTypes.length ||
            filteredCellTypes.includes(cellType.cellTypeName)
          ) {
            result.add(id + cellType.cellTypeName);
          }
        });
      }
    });

    return result;
  }, [
    displayedTissues,
    expandedTissueIds,
    filteredCellTypes,
    sortedCellTypesByTissueName,
  ]);

  const handleExpandCollapse = useCallback(
    (tissueId: string, tissueName: Tissue) => {
      dispatch?.(toggleExpandedTissueId({ tissueId, tissueName }));
    },
    [dispatch]
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
    const cellTypeNames = rawNewFilteredCellTypes.map(
      (cellType) => cellType.name
    );
    const cellTypeIds = cellTypeNames.map((name) => cellTypesByName[name].id);

    dispatch?.(
      setFilteredCellTypes({
        filteredCellTypes: cellTypeNames,
        filteredCellTypeIds: cellTypeIds,
        displayedTissueIds,
      })
    );
  };

  const handleCellTypeDelete = (cellTypeNameToDelete: string) => () => {
    const cellTypeIdToDelete = cellTypesByName[cellTypeNameToDelete].id;
    const newCellTypeNames = filteredCellTypes.filter(
      (cellType) => !(cellTypeNameToDelete === cellType)
    );
    const newCellTypeIds = filteredCellTypeIds.filter(
      (cellTypeId) => !(cellTypeIdToDelete === cellTypeId)
    );

    dispatch?.(
      setFilteredCellTypes({
        filteredCellTypes: newCellTypeNames,
        filteredCellTypeIds: newCellTypeIds,
        displayedTissueIds,
      })
    );
  };

  useTrackHeatMapLoaded({
    selectedGenes: genes,
    displayedCellTypes,
    selectedCellTypes: filteredCellTypes,
  });

  useHandleExpandedTissueIds({
    filteredCellTypeIds,
    filteredTissueIds,
    displayedTissueIds,
    dispatch,
  });

  /**
   * All tissue cell types to render in YAxisCharts
   */
  const allTissueCellTypes = useMemo(() => {
    return displayedTissues
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
    cellTypeSortBy,
    cellTypes,
    sortedCellTypesByTissueName,
    displayedCellTypes,
    displayedTissues,
  ]);

  return (
    <>
      <ContainerWrapper>
        <TopLeftCornerMask height={xAxisHeight}>
          <CellTypeFilterContainer
            id="celltype-filter-container"
            data-testid="celltype-filter"
          >
            <Autocomplete
              className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
              data-testid="cell-type-search"
              multiple
              label="Search cell types"
              value={selectedCellTypeOptions}
              search
              onChange={handleFilteredCellTypesChange}
              options={uniqueCellTypes}
            />
            <CellTypeTagContainer>
              {filteredCellTypes.map((cellType) => (
                <span key={cellType} data-testid={`cell-type-tag-${cellType}`}>
                  <StyledTag
                    label={cellType}
                    onDelete={handleCellTypeDelete(cellType)}
                  />
                </span>
              ))}
            </CellTypeTagContainer>
          </CellTypeFilterContainer>
          <CellCountLabel>Cell Count</CellCountLabel>
        </TopLeftCornerMask>
        <Container {...{ className }} id={HEATMAP_CONTAINER_ID}>
          {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}
          <XAxisWrapper id="x-axis-wrapper">
            <XAxisMask data-testid="x-axis-mask" height={xAxisHeight} />
            <XAxisChart
              geneNames={sortedGeneNames}
              sidebarWidth={sidebarWidth}
            />
          </XAxisWrapper>
          <YAxisWrapper top={0}>
            {allTissueCellTypes.map(
              ({ tissueId, tissueName, tissueCellTypes }) => {
                return (
                  <YAxisChart
                    key={tissueName}
                    tissue={tissueName}
                    tissueID={tissueId}
                    cellTypes={tissueCellTypes}
                    generateMarkerGenes={generateMarkerGenes}
                    expandedTissueIds={expandedTissueIds}
                    handleExpandCollapse={handleExpandCollapse}
                  />
                );
              }
            )}
          </YAxisWrapper>
          <ChartWrapper ref={chartWrapperRef} top={xAxisHeight}>
            {allTissueCellTypes.map(({ tissueName, tissueCellTypes }) => {
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
                  selectedGeneData={selectedGeneData}
                  setIsLoading={setIsLoading}
                  scaledMeanExpressionMax={scaledMeanExpressionMax}
                  scaledMeanExpressionMin={scaledMeanExpressionMin}
                  echartsRendererMode={echartsRendererMode}
                  setAllChartProps={setAllChartProps}
                  chartProps={allChartProps[tissueName]}
                  maxExpression={scaledMeanExpressionMax}
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
