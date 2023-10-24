import { DrawerSize } from "@blueprintjs/core";
import Head from "next/head";
import CellInfoSideBar from "src/views/WheresMyGeneV2/components/CellInfoSideBar";
import Filters from "src/views/WheresMyGeneV2/components/Filters";
import GeneInfoSideBar from "src/components/GeneInfoSideBar";
import SideBar from "src/components/common/SideBar";
import InfoPanel from "src/views/WheresMyGene/components/InfoPanel";
import Legend from "src/views/WheresMyGene/components/InfoPanel/components/Legend";
import Loader from "src/views/WheresMyGene/components/Loader";
import {
  StyledBannerContainer,
  StyledSidebarDrawer,
} from "src/views/WheresMyGene/components/Main/style";
import ScreenTint from "src/views/WheresMyGene/components/ScreenTint";
import {
  SideBarPositioner,
  SideBarWrapper,
  Top,
  Wrapper,
} from "src/views/WheresMyGene/style";
import { View } from "src/views/globalStyle";
import HeatMap from "../HeatMap";
import BottomBanner from "src/components/BottomBanner";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { UnderlyingDataChangeBanner } from "../GeneSearchBar/components/SaveExport/ExportBanner";
import { GENE_EXPRESSION_BANNER_SURVEY_LINK } from "src/common/constants/airtableLinks";
import { StyledRightSideBar } from "./style";
import { useConnect } from "./connect";
import { GENGE_EXPRESSION_PAGE_TITLE } from "./constants";

export default function WheresMyGene(): JSX.Element {
  const {
    set,
    handle,
    selected,
    echartsRendererMode,
    check,
    sortBy,
    geneInfoGene,
    cellInfoCellType,
    filteredCellTypes,
    expandedTissueIds,
    tissuesByID,
    availableFilters,
    sidebarWidth,
    cellTypesByTissueName,
    downloadStatus,
    allChartProps,
    tissuesByName,
    scaledMeanExpressionMax,
    scaledMeanExpressionMin,
  } = useConnect();

  return (
    <>
      <Head>
        <title>{GENGE_EXPRESSION_PAGE_TITLE}</title>
      </Head>
      <SideBar
        label="Filters"
        SideBarWrapperComponent={SideBarWrapper}
        SideBarPositionerComponent={SideBarPositioner}
        testId="filters-panel"
        wmgSideBar
        onWidthChange={set.sidebarWidth}
      >
        <Filters
          isLoading={check.loading}
          availableFilters={availableFilters}
          setAvailableFilters={set.availableFilters}
          setIsScaled={set.isScaled}
        />
      </SideBar>
      {cellInfoCellType && tissuesByID ? (
        <StyledRightSideBar width={CELL_INFO_SIDEBAR_WIDTH_PX}>
          <CellInfoSideBar
            generateGeneInfo={handle.generateGeneInfo}
            cellInfoCellType={cellInfoCellType}
            tissueInfo={tissuesByID[cellInfoCellType.tissueID]}
            handleClose={handle.closeRightSideBar}
            title={`${cellInfoCellType.cellType.name}`}
          />
          {geneInfoGene && (
            <GeneInfoSideBar
              geneInfoGene={geneInfoGene}
              handleClose={handle.closeGeneInfoSideBar}
              title={`${geneInfoGene}`}
            />
          )}
        </StyledRightSideBar>
      ) : (
        geneInfoGene && (
          <StyledRightSideBar>
            <GeneInfoSideBar
              geneInfoGene={geneInfoGene}
              handleClose={handle.closeGeneInfoSideBar}
              title={geneInfoGene}
            />
          </StyledRightSideBar>
        )
      )}
      <View id="view" overflow="hidden">
        <Wrapper>
          {check.loading && <Loader />}

          {downloadStatus.isLoading && (
            <StyledBannerContainer>
              <UnderlyingDataChangeBanner />
            </StyledBannerContainer>
          )}

          <Top id="top-legend">
            <Legend
              selectedCellTypes={cellTypesByTissueName}
              selectedGenes={selected.genes}
              isScaled={check.scaled}
              handleRightSidebarButtonClick={handle.sourceDatasetButtonClick}
              setDownloadStatus={set.downloadStatus}
              setEchartsRendererMode={set.echartsRendererMode}
              allChartProps={allChartProps}
              availableFilters={availableFilters}
              tissues={tissuesByName}
              expandedTissueIds={expandedTissueIds}
              filteredCellTypes={filteredCellTypes}
            />
          </Top>

          <StyledSidebarDrawer
            position="right"
            isOpen={check.sourceDatasetSidebarOpen}
            title="Source Data"
            canEscapeKeyClose={true}
            canOutsideClickClose={true}
            onClose={handle.sourceDatasetButtonClick}
            size={DrawerSize.SMALL}
          >
            <InfoPanel />
          </StyledSidebarDrawer>

          <ScreenTint isDownloading={downloadStatus} />

          <HeatMap
            /**
             * (thuang): Use the selected organism ID to reset the heatmap state
             */
            key={selected.organismId}
            echartsRendererMode={echartsRendererMode}
            cellTypeSortBy={sortBy.cellTypes}
            geneSortBy={sortBy.genes}
            isScaled={check.scaled}
            isLoadingAPI={check.loading}
            cellTypes={cellTypesByTissueName}
            genes={selected.genes}
            selectedGeneExpressionSummariesByTissueName={
              selected.geneExpressionSummariesByTissueName
            }
            scaledMeanExpressionMax={scaledMeanExpressionMax}
            scaledMeanExpressionMin={scaledMeanExpressionMin}
            allChartProps={allChartProps}
            setAllChartProps={set.allChartProps}
            setTissuesByName={set.tissuesByName}
            tissuesByName={tissuesByName}
            /**
             * (thuang): This is needed to reposition gene search bar when the
             * sidebar width expands/collapses
             */
            sidebarWidth={sidebarWidth}
          />
        </Wrapper>
      </View>

      <BottomBanner
        airtableLink={GENE_EXPRESSION_BANNER_SURVEY_LINK}
        includeSurveyLink
      />
    </>
  );
}
