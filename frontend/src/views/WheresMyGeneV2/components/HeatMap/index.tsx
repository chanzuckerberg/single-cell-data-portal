import { memo } from "react";

import { Tissue } from "src/views/WheresMyGeneV2/common/types";
import YAxisChart from "./components/YAxisChart";

import {
  CellTypeTagContainer,
  ChartWrapper,
  Container,
  ContainerWrapper,
  SkeletonContainer,
  SkeletonWrapper,
  StyledTag,
  XAxisMask,
  YAxisWrapper,
} from "src/views/WheresMyGeneV2/components/HeatMap/style";
import { CellCountLabel } from "src/views/WheresMyGeneV2/components/HeatMap/components/XAxisChart/style";
import {
  HEATMAP_CONTAINER_ID,
  MARGIN_BETWEEN_HEATMAPS,
} from "src/views/WheresMyGeneV2/common/constants";
import Loader from "src/views/WheresMyGeneV2/components/Loader";
import XAxisChart from "src/views/WheresMyGeneV2/components/HeatMap/components/XAxisChart";
import Chart from "src/views/WheresMyGeneV2/components/HeatMap/components/Chart";
import { hyphenize } from "src/views/WheresMyGeneV2/components/HeatMap/utils";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";
import { Autocomplete } from "@czi-sds/components";
import {
  CellTypeFilterContainer,
  Divider,
  TopLeftCornerMask,
  XAxisWrapper,
  StyledSkeleton,
} from "./style";

import { useConnect } from "src/views/WheresMyGeneV2/components/HeatMap/connect";

import { Props } from "./types";

export default memo(function HeatMap(props: Props): JSX.Element {
  const {
    allChartProps,
    className,
    echartsRendererMode,
    isLoadingAPI,
    isScaled,
    scaledMeanExpressionMax,
    scaledMeanExpressionMin,
    setAllChartProps,
    sidebarWidth,
  } = props;

  const {
    allTissueCellTypes,
    chartWrapperRef,
    expandedTissueIds,
    filteredCellTypes,
    generateMarkerGenes,
    handleCellTypeDelete,
    handleExpandCollapse,
    handleFilteredCellTypesChange,
    isLoading,
    orderedSelectedGeneExpressionSummariesByTissueName,
    selectedCellTypeOptions,
    setIsLoading,
    sortedGeneNames,
    totalElementsCount,
    uniqueCellTypes,
    xAxisHeight,
  } = useConnect(props);

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
          <YAxisWrapper top={xAxisHeight}>
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
          {isLoadingAPI || isAnyTissueLoading(isLoading) ? (
            <ChartWrapper
              top={xAxisHeight}
              hidden={!isAnyTissueLoading(isLoading)}
            >
              <SkeletonContainer>
                {[...Array(totalElementsCount)].map((i) => (
                  <SkeletonWrapper key={i}>
                    {[...Array(sortedGeneNames.length)].map((i) => (
                      <StyledSkeleton
                        variant="circular"
                        width={19}
                        height={19}
                        key={i}
                      />
                    ))}
                  </SkeletonWrapper>
                ))}
              </SkeletonContainer>
            </ChartWrapper>
          ) : null}
          <ChartWrapper
            ref={chartWrapperRef}
            top={xAxisHeight}
            hidden={isAnyTissueLoading(isLoading)}
          >
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

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}
