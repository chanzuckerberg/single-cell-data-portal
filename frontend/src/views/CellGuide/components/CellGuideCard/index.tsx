import React, { useMemo } from "react";
import { Global } from "@emotion/react";
import {
  Wrapper,
  CellGuideCardName,
  CellGuideCardHeader,
  StyledTag,
  CellGuideView,
  CellGuideCardHeaderInnerWrapper,
  StyledRightSideBar,
  MobileTooltipTitle,
  MobileTooltipWrapper,
  MobileTooltipHeader,
  NavBarDropdownWrapper,
  CellGuideWrapper,
  StyledCellTagSideBar,
  StyledGeneTagSideBar,
  StyledAutocomplete,
} from "./style";
import Description from "./components/Description";
import MarkerGeneTables from "./components/MarkerGeneTables";
import OntologyDagView from "../common/OntologyDagView";
import FullScreenProvider from "../common/FullScreenProvider";
import SourceDataTable from "./components/SourceDataTable";
import CellGuideCardSidebar from "./components/CellGuideCardSidebar";
import CellGuideMobileHeader from "../CellGuideMobileHeader";
import GeneInfoSideBar from "src/components/GeneInfoSideBar";
import { titleize } from "src/common/utils/string";
import Head from "next/head";
import CellGuideBottomBanner from "../CellGuideBottomBanner";
import { StickySidebarStyle } from "./components/CellGuideCardSidebar/style";
import {
  CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_GLOBAL_MARKER_GENE_DROPDOWN,
  CELL_GUIDE_CARD_HEADER_NAME,
  CELL_GUIDE_CARD_HEADER_TAG,
  RIGHT_SIDEBAR_WIDTH_PX,
  SELECT_A_GENE,
} from "src/views/CellGuide/components/CellGuideCard/constants";
import {
  ALL_TISSUES,
  TISSUE_AGNOSTIC,
} from "./components/MarkerGeneTables/constants";
import {
  Dropdown,
  InputDropdownProps,
  Button,
  DefaultAutocompleteOption,
} from "@czi-sds/components";
import { useComponentWidth } from "./components/common/hooks/useIsComponentPastBreakpoint";
import { DEFAULT_ONTOLOGY_HEIGHT } from "../common/OntologyDagView/common/constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import CellGuideInfoSideBar from "../CellGuideInfoSideBar";
import { useComputationalMarkerGenesTableRowsAndFilters } from "./components/MarkerGeneTables/hooks/computational_markers";
import { useConnect } from "./connect";
import { SDSOrgan } from "src/views/CellGuide/components/CellGuideCard/types";
import { getCellTypeLink } from "src/views/CellGuide/common/utils";
import { AutocompleteValue, TextField } from "@mui/material";

export const SDS_INPUT_DROPDOWN_PROPS: InputDropdownProps = {
  sdsStyle: "square",
} as InputDropdownProps;

interface Props {
  name: string;
  seoDescription: string;
  synonyms?: string[];
}

export default function CellGuideCard({
  // From getServerSideProps
  name,
  // From getServerSideProps
  seoDescription: rawSeoDescription,
  // From getServerSideProps
  synonyms,
}: Props): JSX.Element {
  const {
    router,
    pageNavIsOpen,
    setPageNavIsOpen,
    selectedGene,
    headerRef,
    topPadding,
    sectionRef0,
    sectionRef1,
    sectionRef2,
    sectionRef3,
    selectGene,
    setSelectedGene,
    skinnyMode,
    tooltipContent,
    setTooltipContent,
    queryCellTypeId,
    cellTypeId,
    geneInfoGene,
    setGeneInfoGene,
    cellInfoCellType,
    setCellInfoCellType,
    sdsOrganismsList,
    sdsOrgansList,
    selectedOrgan,
    selectedOrganId,
    selectedOrganism,
    setSelectedOrganism,
  } = useConnect();

  const tissueName = selectedOrgan?.name || "";

  const cellGuideSideBar = useMemo(() => {
    return (
      <CellGuideCardSidebar
        sectionClickHandler={() => setPageNavIsOpen(false)}
        skinnyMode={skinnyMode}
        items={[
          { elementRef: sectionRef0, title: "Intro" },
          { elementRef: sectionRef1, title: "Cell Ontology" },
          { elementRef: sectionRef2, title: "Marker Genes" },
          { elementRef: sectionRef3, title: "Data" },
        ]}
      />
    );
  }, [
    skinnyMode,
    sectionRef0,
    sectionRef1,
    sectionRef2,
    sectionRef3,
    setPageNavIsOpen,
  ]);

  const cellTypeName = name || "";
  const titleizedCellTypeName = titleize(cellTypeName);

  const handleChangeOrgan = (
    _: React.SyntheticEvent, // event
    option: AutocompleteValue<DefaultAutocompleteOption, false, false, false>
  ) => {
    const { name, id } = (option || {}) as SDSOrgan;

    if (!option || !name || name === tissueName) return;

    const optionName = name === TISSUE_AGNOSTIC ? ALL_TISSUES : name;

    // Continue tracking the analytics event as All Tissues
    track(EVENTS.CG_SELECT_TISSUE, { tissue: optionName });

    const url = getCellTypeLink({ tissueId: id, cellTypeId });

    /**
     * (thuang): Product requirement to keep the scroll position
     */
    router.push(url, url, { scroll: false });
  };

  const handleChangeOrganism = (
    _: React.SyntheticEvent,
    option: AutocompleteValue<DefaultAutocompleteOption, false, false, false>
  ) => {
    if (!option || option.name === selectedOrganism.name) return;
    setSelectedOrganism(option);
    track(EVENTS.CG_SELECT_ORGANISM, { organism: option.name });
  };

  function handleCloseGeneInfoSideBar() {
    setGeneInfoGene(null);
  }
  function handleCloseCellGuideInfoSideBar() {
    setGeneInfoGene(null);
    setCellInfoCellType(null);
  }

  const cellTypePrefix =
    tissueName === TISSUE_AGNOSTIC ? "" : `${tissueName} specific `;

  const { computationalMarkerGeneTableData } =
    useComputationalMarkerGenesTableRowsAndFilters({
      cellTypeId,
      organismName: selectedOrganism.name,
      organId: selectedOrganId,
    });

  const genesList = useMemo(() => {
    return computationalMarkerGeneTableData.map((gene) => gene.symbol);
  }, [computationalMarkerGeneTableData]);

  const handleChangeGene = (option: string) => {
    if (!option) {
      setSelectedGene(undefined);
    } else if (option) {
      setSelectedGene(option);
    }
  };

  const geneDropdownComponent = (
    <StyledAutocomplete
      data-testid={CELL_GUIDE_CARD_GLOBAL_MARKER_GENE_DROPDOWN}
      options={genesList}
      renderInput={(params) => (
        <TextField {...params} placeholder={SELECT_A_GENE} />
      )}
      onChange={(_, value) => {
        handleChangeGene(value as string);
        track(EVENTS.CG_SELECT_MARKER_GENE, { marker_gene: value });
      }}
      value={selectedGene}
    />
  );
  const title = `${titleize(
    cellTypePrefix
  )}${titleizedCellTypeName} Cell Types - CZ CELLxGENE CellGuide`;

  const seoDescription = `Find comprehensive information about ${cellTypePrefix}"${cellTypeName}" cell types (synonyms: ${
    synonyms?.join(", ") || "N/A"
  }). ${rawSeoDescription}`;

  const OrganismSelectorDropdown = (
    <Dropdown<DefaultAutocompleteOption, false, false, false>
      InputDropdownProps={SDS_INPUT_DROPDOWN_PROPS}
      search
      label={selectedOrganism?.name}
      onChange={handleChangeOrganism}
      options={sdsOrganismsList}
      value={selectedOrganism}
      data-testid={CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN}
    />
  );
  const dropdownComponents = (
    <CellGuideCardHeaderInnerWrapper>
      {OrganismSelectorDropdown}
      <Dropdown<DefaultAutocompleteOption, false, false, false>
        InputDropdownProps={SDS_INPUT_DROPDOWN_PROPS}
        search
        label={tissueName}
        onChange={handleChangeOrgan}
        options={sdsOrgansList}
        value={selectedOrgan}
        data-testid={CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN}
      />
    </CellGuideCardHeaderInnerWrapper>
  );
  const pageNav = (
    <div>
      {cellGuideSideBar}
      <NavBarDropdownWrapper>{dropdownComponents}</NavBarDropdownWrapper>
    </div>
  );
  const { width, containerRef } = useComponentWidth();

  const geneInfoGeneTitle = (
    <span>
      {geneInfoGene}{" "}
      <StyledGeneTagSideBar
        label="Gene"
        sdsType="secondary"
        sdsStyle="square"
        color="info"
        tagColor="info"
      />
    </span>
  );

  return (
    <>
      {/* This is a fix that overrides a global overflow css prop to get sticky elements to work */}
      <Global styles={StickySidebarStyle} />

      <Head>
        <title>{title}</title>
        <meta property="title" key="title" content={title} />
        <meta property="og:title" key="og:title" content={title} />
        <meta property="twitter:title" key="twitter:title" content={title} />

        <meta name="description" key="description" content={seoDescription} />
        <meta
          property="og:description"
          key="og:description"
          content={seoDescription}
        />
        <meta
          property="twitter:description"
          key="twitter:description"
          content={seoDescription}
        />

        {/* This prevents auto zooming on the input box on mobile */}
        <meta
          name="viewport"
          content="width=device-width, initial-scale=1, minimum-scale=1"
        />
      </Head>

      {/* Cell Guide Mobile Navigation */}
      {skinnyMode && (
        <CellGuideMobileHeader
          title={titleizedCellTypeName}
          pageNav={pageNav}
          pageNavIsOpen={pageNavIsOpen}
          setPageNavIsOpen={setPageNavIsOpen}
        />
      )}
      {/* when tooltipContent is null, remove the tooltip */}
      {/* setTooltipContent sets the title and content element */}
      {skinnyMode && tooltipContent && (
        <MobileTooltipWrapper>
          <MobileTooltipHeader>
            <MobileTooltipTitle>{tooltipContent.title}</MobileTooltipTitle>
            <Button
              onClick={() => {
                setTooltipContent(null);
              }}
              icon="XMark"
              sdsStyle="icon"
              sdsSize="medium"
            />
          </MobileTooltipHeader>
          <div>{tooltipContent.element}</div>
        </MobileTooltipWrapper>
      )}
      <div>
        {/* Intro section */}
        <div ref={sectionRef0} id="section-0" data-testid="section-0" />
        {/* Don't show title of the cell card if we're on mobile, since the title is already in the header nav */}
        {!skinnyMode && (
          <CellGuideCardHeader ref={headerRef} width={width}>
            <CellGuideCardHeaderInnerWrapper>
              <CellGuideCardName data-testid={CELL_GUIDE_CARD_HEADER_NAME}>
                {titleizedCellTypeName}
              </CellGuideCardName>
              <a
                href={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${queryCellTypeId}`}
                target="_blank"
                rel="noreferrer noopener"
              >
                <StyledTag
                  data-testid={CELL_GUIDE_CARD_HEADER_TAG}
                  label={cellTypeId}
                  sdsType="secondary"
                  sdsStyle="square"
                  color="neutral"
                  hover
                />
              </a>
            </CellGuideCardHeaderInnerWrapper>
            {dropdownComponents}
          </CellGuideCardHeader>
        )}
      </div>
      <CellGuideWrapper skinnyMode={skinnyMode}>
        <CellGuideView skinnyMode={skinnyMode}>
          {/* Flex item left */}
          <Wrapper
            topMargin={topPadding}
            skinnyMode={skinnyMode}
            ref={containerRef}
          >
            {/* (thuang): Somehow we need a parent to prevent error:
              NotFoundError: Failed to execute 'insertBefore' on 'Node'
            */}
            <Description
              cellTypeId={cellTypeId}
              cellTypeName={cellTypeName}
              skinnyMode={skinnyMode}
              setTooltipContent={setTooltipContent}
              synonyms={synonyms}
            />

            {/* Cell Ontology section */}
            <div ref={sectionRef1} id="section-1" data-testid="section-1" />
            {/* (thuang): Somehow we need a parent <div /> to prevent error:
          NotFoundError: Failed to execute 'insertBefore' on 'Node'
         */}
            <div>
              <FullScreenProvider cellInfoSideBarDisplayed={!!cellInfoCellType}>
                <OntologyDagView
                  key={`${cellTypeId}-${selectedOrganId}-${selectedGene}`}
                  cellTypeId={cellTypeId}
                  cellTypeName={cellTypeName}
                  tissueName={tissueName}
                  tissueId={selectedOrganId}
                  inputWidth={width}
                  setCellInfoCellType={setCellInfoCellType}
                  inputHeight={DEFAULT_ONTOLOGY_HEIGHT}
                  selectedOrganism={selectedOrganism.name}
                  selectedGene={selectedGene}
                  geneDropdownComponent={geneDropdownComponent}
                />
              </FullScreenProvider>
            </div>

            {/* Marker Genes section */}
            <div ref={sectionRef2} id="section-2" data-testid="section-2" />
            <MarkerGeneTables
              setTooltipContent={setTooltipContent}
              key={cellTypeId}
              cellTypeId={cellTypeId}
              setGeneInfoGene={setGeneInfoGene}
              cellTypeName={cellTypeName}
              skinnyMode={skinnyMode}
              organName={tissueName}
              organId={selectedOrganId}
              organismName={selectedOrganism.name}
              selectedGene={selectedGene}
              selectGene={selectGene}
            />

            {/* Source Data section */}
            <div ref={sectionRef3} id="section-3" data-testid="section-3" />
            <SourceDataTable
              cellTypeId={cellTypeId}
              organId={selectedOrganId}
              organismName={selectedOrganism.name}
              skinnyMode={skinnyMode}
              setTooltipContent={setTooltipContent}
            />
          </Wrapper>

          {/* Side bar */}
          {!skinnyMode && cellGuideSideBar}
        </CellGuideView>
      </CellGuideWrapper>
      {cellInfoCellType ? (
        <StyledRightSideBar
          width={RIGHT_SIDEBAR_WIDTH_PX}
          skinnyMode={skinnyMode}
        >
          <CellGuideInfoSideBar
            cellInfoCellType={cellInfoCellType}
            setCellInfoCellType={setCellInfoCellType}
            handleClose={handleCloseCellGuideInfoSideBar}
            title={
              <span>
                {titleize(cellInfoCellType.cellTypeName)}{" "}
                <StyledCellTagSideBar
                  label="Cell Type"
                  sdsType="secondary"
                  sdsStyle="square"
                  color="info"
                  tagColor="info"
                />
              </span>
            }
            setGeneInfoGene={setGeneInfoGene}
            selectedOrganName={tissueName}
            selectedOrganId={selectedOrganId}
            organismName={selectedOrganism.name}
            selectedGene={selectedGene}
            selectGene={selectGene}
            setTooltipContent={setTooltipContent}
            skinnyMode={skinnyMode}
          />

          {
            // Split right sidebar view if fmg AND gene info is populated
            !!geneInfoGene && (
              <GeneInfoSideBar
                geneInfoGene={geneInfoGene}
                handleClose={handleCloseGeneInfoSideBar}
                title={geneInfoGeneTitle}
              />
            )
          }
        </StyledRightSideBar>
      ) : (
        // Gene info full right sidebar length
        geneInfoGene && (
          <StyledRightSideBar
            width={RIGHT_SIDEBAR_WIDTH_PX}
            skinnyMode={skinnyMode}
          >
            <GeneInfoSideBar
              geneInfoGene={geneInfoGene}
              handleClose={handleCloseGeneInfoSideBar}
              title={geneInfoGeneTitle}
            />
          </StyledRightSideBar>
        )
      )}
      {/* dont include long survey link text if in mobile view */}
      <CellGuideBottomBanner includeSurveyLink={!skinnyMode} />
    </>
  );
}
