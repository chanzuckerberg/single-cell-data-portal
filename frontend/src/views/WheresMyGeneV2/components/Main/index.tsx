import { DrawerSize } from "@blueprintjs/core";
import Head from "next/head";

import SideBar from "src/components/common/SideBar";

import {
  addGeneInfoGene,
  clearGeneInfoGene,
  closeRightSidebar,
} from "src/views/WheresMyGeneV2/common/store/actions";

import CellInfoSideBar from "src/views/WheresMyGeneV2/components/CellInfoSideBar";
import Filters from "src/views/WheresMyGeneV2/components/Filters";
import GeneInfoSideBar from "src/components/GeneInfoSideBar";

import InfoPanel from "src/views/WheresMyGeneV2/components/InfoPanel";
import Legend from "src/views/WheresMyGeneV2/components/InfoPanel/components/Legend";
import Loader from "src/views/WheresMyGeneV2/components/Loader";
import {
  StyledBannerContainer,
  StyledSidebarDrawer,
} from "src/views/WheresMyGeneV2/components/Main/style";
import ScreenTint from "src/views/WheresMyGeneV2/components/ScreenTint";
import {
  SideBarPositioner,
  SideBarWrapper,
  Top,
  Wrapper,
} from "src/views/WheresMyGeneV2/style";
import { View } from "src/views/globalStyle";
import HeatMap from "../HeatMap";
import BottomBanner from "src/components/BottomBanner";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { UnderlyingDataChangeBanner } from "../GeneSearchBar/components/SaveExport/ExportBanner";
import { GENE_EXPRESSION_BANNER_SURVEY_LINK } from "src/common/constants/airtableLinks";
import { StyledRightSideBar } from "./style";
import { useConnect } from "src/views/WheresMyGeneV2/components/Main/connect";
import Notification from "src/components/Notification";

export const INFO_PANEL_WIDTH_PX = 320;

export default function WheresMyGene(): JSX.Element {
  const {
    dispatch,
    sidebarWidth,
    setSidebarWidth,
    isLoading,
    availableFilters,
    setAvailableFilters,
    isScaled,
    setIsScaled,
    handleSourceDatasetButtonClick,
    downloadStatus,
    setDownloadStatus,
    echartsRendererMode,
    setEchartsRendererMode,
    allChartProps,
    setAllChartProps,
    tissuesByName,
    setTissuesByName,
    sortBy,
    geneInfoGene,
    cellInfoCellType,
    tissuesByID,
    cellTypesByTissueName,
    selectedGenes,
    expandedTissueIds,
    filteredCellTypes,
    scaledMeanExpressionMax,
    scaledMeanExpressionMin,
    isSourceDatasetSidebarOpen,
    selectedOrganismId,
    selectedGeneExpressionSummariesByTissueName,
  } = useConnect();

  const handleCloseRightSideBar = () => {
    if (!dispatch) return;
    dispatch(closeRightSidebar());
  };

  const handleCloseGeneInfoSideBar = () => {
    if (!dispatch) return;
    dispatch(clearGeneInfoGene());
  };

  const generateGeneInfo = (gene: string) => {
    if (!dispatch) return;
    dispatch(addGeneInfoGene(gene));
  };

  return (
    <>
      <Head>
        <title>Gene Expression - CZ CELLxGENE Discover</title>
      </Head>

      <SideBar
        label="Filters"
        SideBarWrapperComponent={SideBarWrapper}
        SideBarPositionerComponent={SideBarPositioner}
        testId="filters-panel"
        wmgSideBar
        onWidthChange={setSidebarWidth}
      >
        <Filters
          isLoading={isLoading}
          availableFilters={availableFilters}
          setAvailableFilters={setAvailableFilters}
          setIsScaled={setIsScaled}
        />
      </SideBar>
      {cellInfoCellType && tissuesByID ? (
        <StyledRightSideBar width={CELL_INFO_SIDEBAR_WIDTH_PX}>
          <CellInfoSideBar
            generateGeneInfo={generateGeneInfo}
            cellInfoCellType={cellInfoCellType}
            tissueInfo={tissuesByID[cellInfoCellType.tissueID]}
            handleClose={handleCloseRightSideBar}
            title={`${cellInfoCellType.cellType.name}`}
          />

          {
            // Split right sidebar view if fmg AND gene info is populated
            !!geneInfoGene && (
              <GeneInfoSideBar
                geneInfoGene={geneInfoGene}
                handleClose={handleCloseGeneInfoSideBar}
                title={`${geneInfoGene}`}
              />
            )
          }
        </StyledRightSideBar>
      ) : (
        // Gene info full right sidebar length
        geneInfoGene && (
          <StyledRightSideBar>
            <GeneInfoSideBar
              geneInfoGene={geneInfoGene}
              handleClose={handleCloseGeneInfoSideBar}
              title={geneInfoGene}
            />
          </StyledRightSideBar>
        )
      )}

      <View id="view" overflow="hidden">
        <Wrapper>
          {isLoading && <Loader />}
          {/* Used for PNG and SVG exports to render message banner to render in output */}
          {downloadStatus.isLoading && (
            <StyledBannerContainer>
              <UnderlyingDataChangeBanner />
            </StyledBannerContainer>
          )}
          <Top id="top-legend">
            <Notification />
            <Legend
              selectedCellTypes={cellTypesByTissueName}
              selectedGenes={selectedGenes}
              isScaled={isScaled}
              handleRightSidebarButtonClick={handleSourceDatasetButtonClick}
              setDownloadStatus={setDownloadStatus}
              setEchartsRendererMode={setEchartsRendererMode}
              allChartProps={allChartProps}
              availableFilters={availableFilters}
              tissues={tissuesByName}
              expandedTissueIds={expandedTissueIds}
              filteredCellTypes={filteredCellTypes}
              maxExpression={scaledMeanExpressionMax}
            />
          </Top>
          <StyledSidebarDrawer
            position="right"
            isOpen={isSourceDatasetSidebarOpen}
            title="Source Data"
            canEscapeKeyClose={true}
            canOutsideClickClose={true}
            onClose={handleSourceDatasetButtonClick}
            size={DrawerSize.SMALL}
          >
            <InfoPanel />
          </StyledSidebarDrawer>

          <ScreenTint isDownloading={downloadStatus} />

          <HeatMap
            /**
             * (thuang): Use the selected organism ID to reset the heatmap state
             */
            key={selectedOrganismId}
            echartsRendererMode={echartsRendererMode}
            cellTypeSortBy={sortBy.cellTypes}
            geneSortBy={sortBy.genes}
            isScaled={isScaled}
            isLoadingAPI={isLoading}
            cellTypes={cellTypesByTissueName}
            genes={selectedGenes}
            selectedGeneExpressionSummariesByTissueName={
              selectedGeneExpressionSummariesByTissueName
            }
            scaledMeanExpressionMax={scaledMeanExpressionMax}
            scaledMeanExpressionMin={scaledMeanExpressionMin}
            allChartProps={allChartProps}
            setAllChartProps={setAllChartProps}
            setTissuesByName={setTissuesByName}
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
