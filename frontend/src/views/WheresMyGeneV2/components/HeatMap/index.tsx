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
import SearchIcon from "@mui/icons-material/Search";
import { EMPTY_ARRAY } from "src/common/constants/utils";
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
import { addCellInfoCellType } from "src/views/WheresMyGene/common/store/actions";
import {
  CellType,
  ChartProps,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGene/common/types";
import YAxisChart from "./components/YAxisChart";
import { useTrackHeatMapLoaded } from "./hooks/useTrackHeatMapLoaded";
import { memoize } from "lodash";
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
  StyledAutocomplete,
  StyledTag,
  TopLeftCornerMask,
  XAxisMask,
  XAxisWrapper,
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
import { InputAdornment, TextField } from "@mui/material";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";

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
  selectedOrganismId: string;
  echartsRendererMode: "svg" | "canvas";
  setAllChartProps: Dispatch<
    SetStateAction<{
      [tissue: string]: ChartProps;
    }>
  >;
  allChartProps: { [tissue: string]: ChartProps };
  setTissuesByName: Dispatch<
    SetStateAction<{
      [name: string]: OntologyTerm;
    }>
  >;
  expandedTissues: Set<string>;
  setExpandedTissues: Dispatch<SetStateAction<Set<string>>>;
  filteredCellTypes: string[];
  setFilteredCellTypes: Dispatch<SetStateAction<string[]>>;
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
  selectedOrganismId,
  echartsRendererMode,
  allChartProps,
  setAllChartProps,
  setTissuesByName,
  expandedTissues,
  setExpandedTissues,
  filteredCellTypes,
  setFilteredCellTypes,
}: Props): JSX.Element {
  const {
    xAxisHeight,
    selectedFilters: { tissues: filteredTissues },
  } = useContext(StateContext);
  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);
  const dispatch = useContext(DispatchContext);

  const { data } = usePrimaryFilterDimensions(2); //temp explicit version

  // Get tissueName to ID map for use in find marker genes
  const tissuesByName = useMemo(() => {
    let result: { [name: string]: OntologyTerm } = {};

    if (!data) return result;

    const { tissues } = data;

    result = generateTermsByKey(tissues, "name");

    setTissuesByName(result);

    return result;
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

  const initialDisplayedCellTypes = useMemo(
    () =>
      Object.entries(tissuesByName).reduce((acc, [_, { id }]) => {
        if (
          (filteredTissues.length > 0 && filteredTissues.includes(id)) ||
          filteredTissues.length === 0
        ) {
          acc.add(id + id);
        }
        return acc;
      }, new Set<string>()),
    [tissuesByName, filteredTissues]
  );

  // set of tissue names that are visible and set of cell types that are visible
  // presence is visible
  const [displayedCellTypes, setDisplayedCellTypes] = useState<Set<string>>(
    initialDisplayedCellTypes
  );

  useEffect(() => {
    setDisplayedCellTypes(initialDisplayedCellTypes);
    setExpandedTissues(new Set<string>());
  }, [initialDisplayedCellTypes, setExpandedTissues, selectedOrganismId]);

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
    return [...result].sort();
  }, [sortedCellTypesByTissueName]);

  // update displayedCellTypes and expandedTissues
  const handleFilteredCellTypesChange = (
    _: unknown,
    rawNewFilteredCellTypes: unknown
  ) => {
    const newFilteredCellTypes = rawNewFilteredCellTypes as string[];
    if (newFilteredCellTypes.length === 0) {
      setDisplayedCellTypes(initialDisplayedCellTypes);
      setExpandedTissues(new Set<string>());
      setFilteredCellTypes(newFilteredCellTypes);
      return;
    }
    const newDisplayedCellTypes = new Set<string>();
    const newExpandedTissues = new Set<string>();
    Object.entries(sortedCellTypesByTissueName).forEach(
      ([tissue, cellTypes]) => {
        if (
          filteredTissues.length > 0 &&
          !filteredTissues.includes(tissuesByName[tissue].id)
        )
          return;
        cellTypes.forEach((cellType) => {
          if (newFilteredCellTypes.includes(cellType.name)) {
            newDisplayedCellTypes.add(
              tissuesByName[tissue].id + tissuesByName[tissue].id
            );
            newDisplayedCellTypes.add(tissuesByName[tissue].id + cellType.name);
            newExpandedTissues.add(tissuesByName[tissue].id);
          }
        });
      }
    );
    if (newFilteredCellTypes.length > filteredCellTypes.length) {
      const filteredCellTypeIDs = newFilteredCellTypes.map(
        (cellType) => cellTypesByName[cellType].id
      );
      track(EVENTS.WMG_SELECT_CELL_TYPE, {
        cell_types: filteredCellTypeIDs,
      });
    }

    setDisplayedCellTypes(newDisplayedCellTypes);
    setExpandedTissues(newExpandedTissues);
    setFilteredCellTypes(newFilteredCellTypes);
  };

  const handleCellTypeDelete = (cellTypeToDelete: string) => () => {
    const newValue = filteredCellTypes.filter(
      (cellType) => !(cellTypeToDelete === cellType)
    );
    handleFilteredCellTypesChange(null, newValue);
  };

  useTrackHeatMapLoaded({
    selectedGenes: genes,
    displayedCellTypes,
    selectedCellTypes: filteredCellTypes,
  });

  return (
    <>
      <ContainerWrapper>
        <TopLeftCornerMask height={xAxisHeight}>
          <CellTypeFilterContainer
            id="celltype-filter-container"
            className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
          >
            <StyledAutocomplete
              multiple
              value={filteredCellTypes}
              onChange={handleFilteredCellTypesChange}
              renderInput={(params) => (
                <TextField
                  {...params}
                  InputProps={{
                    ...params.InputProps,
                    startAdornment: (
                      <InputAdornment position="start">
                        <SearchIcon />
                      </InputAdornment>
                    ),
                    endAdornment: undefined,
                  }}
                  placeholder="Search cell types"
                ></TextField>
              )}
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
            {Object.values(tissuesByName)
              .sort((a, b) =>
                // sort tissues alphabetically
                a.name.localeCompare(b.name)
              )
              .map((tissue: OntologyTerm) => {
                const tissueCellTypes = memoizedGetTissueCellTypes({
                  cellTypeSortBy,
                  cellTypes,
                  sortedCellTypesByTissueName,
                  tissue: tissue.name,
                  tissueID: tissue.id,
                  displayedCellTypes,
                });

                return (
                  tissueCellTypes.length > 0 && (
                    <YAxisChart
                      key={tissue.name}
                      tissue={tissue.name}
                      tissueID={tissue.id}
                      cellTypes={tissueCellTypes}
                      generateMarkerGenes={generateMarkerGenes}
                      selectedOrganismId={selectedOrganismId}
                      expandedTissues={expandedTissues}
                      handleExpandCollapse={handleExpandCollapse}
                    />
                  )
                );
              })}
          </YAxisWrapper>
          <ChartWrapper ref={chartWrapperRef} top={xAxisHeight}>
            {Object.values(tissuesByName)
              .sort((a, b) =>
                // sort tissues alphabetically
                a.name.localeCompare(b.name)
              )
              .map((tissue: OntologyTerm) => {
                const tissueCellTypes = memoizedGetTissueCellTypes({
                  cellTypeSortBy,
                  cellTypes,
                  sortedCellTypesByTissueName,
                  tissue: tissue.name,
                  tissueID: tissue.id,
                  displayedCellTypes: displayedCellTypes,
                });

                // Don't render anything if tissue has no cell types for some reason
                if (!tissueCellTypes.length) return;

                const selectedGeneData =
                  orderedSelectedGeneExpressionSummariesByTissueName[
                    tissue.name
                  ];

                /**
                 * (thuang): If there is no selected gene data, we don't want to render
                 * the chart, because it will cause the chart to render with 0 width,
                 * which is an error for echarts
                 */
                if (!selectedGeneData?.length) {
                  const height =
                    document.getElementById(`${hyphenize(tissue.name)}-y-axis`)
                      ?.clientHeight ?? 0;
                  return (
                    <div
                      id={`no-chart-data-${hyphenize(tissue.name)}`} // Not used, just to make it stand out
                      key={`${tissue.name}-${echartsRendererMode}`}
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
                    key={`${tissue.name}-${echartsRendererMode}`}
                    tissue={tissue.name}
                    cellTypes={tissueCellTypes}
                    selectedGeneData={
                      orderedSelectedGeneExpressionSummariesByTissueName[
                        tissue.name
                      ]
                    }
                    setIsLoading={setIsLoading}
                    scaledMeanExpressionMax={scaledMeanExpressionMax}
                    scaledMeanExpressionMin={scaledMeanExpressionMin}
                    echartsRendererMode={echartsRendererMode}
                    setAllChartProps={setAllChartProps}
                    chartProps={allChartProps[tissue.name]}
                  />
                );
              })}
          </ChartWrapper>
        </Container>
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

  if (!tissueCellTypes || tissueCellTypes.length === 0) return [];
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

const memoizedGetTissueCellTypes = memoize(
  getTissueCellTypes,
  ({
    cellTypes,
    sortedCellTypesByTissueName,
    tissue,
    tissueID,
    cellTypeSortBy,
    displayedCellTypes,
  }) => {
    return `${tissue}-${tissueID}-${cellTypeSortBy}-${[
      ...displayedCellTypes,
    ]?.join("")}-${cellTypes[tissue]
      ?.map((cellType) => cellType.viewId)
      .join("")}-${sortedCellTypesByTissueName[tissue]?.join("")}`;
  }
);

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
