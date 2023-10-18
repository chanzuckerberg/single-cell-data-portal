import React, { useCallback, useEffect, useMemo, useState } from "react";
import { useRouter } from "next/router";
import { Global } from "@emotion/react";
import {
  Wrapper,
  CellGuideCardName,
  CellGuideCardHeader,
  StyledTag,
  CellGuideView,
  CellGuideCardHeaderInnerWrapper,
  LEFT_RIGHT_PADDING_PX_XXL,
  StyledRightSideBar,
  MobileTooltipTitle,
  MobileTooltipWrapper,
  MobileTooltipHeader,
  NavBarDropdownWrapper,
  CellGuideWrapper,
  StyledCellTagSideBar,
  StyledGeneTagSideBar,
} from "./style";
import Description from "./components/Description";
import MarkerGeneTables from "./components/MarkerGeneTables";
import OntologyDagView from "../common/OntologyDagView";
import FullScreenProvider from "../common/FullScreenProvider";
import SourceDataTable from "./components/SourceDataTable";
import CellGuideCardSidebar from "./components/CellGuideCardSidebar";
import CellGuideMobileHeader from "../CellGuideMobileHeader";
import { Gene } from "src/views/WheresMyGene/common/types";
import { throttle } from "lodash";
import GeneInfoSideBar from "src/components/GeneInfoSideBar";
import { titleize } from "src/common/utils/string";
import Head from "next/head";
import CellGuideBottomBanner from "../CellGuideBottomBanner";
import { StickySidebarStyle } from "./components/CellGuideCardSidebar/style";
import {
  SKINNY_MODE_BREAKPOINT_WIDTH,
  CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_HEADER_NAME,
  CELL_GUIDE_CARD_HEADER_TAG,
  RIGHT_SIDEBAR_WIDTH_PX,
} from "src/views/CellGuide/components/CellGuideCard/constants";
import { useOrganAndOrganismFilterListForCellType } from "./components/MarkerGeneTables/hooks/common";
import {
  ALL_TISSUES,
  NO_ORGAN_ID,
  TISSUE_AGNOSTIC,
} from "./components/MarkerGeneTables/constants";
import {
  DefaultDropdownMenuOption,
  Dropdown,
  InputDropdownProps,
  ButtonIcon,
} from "@czi-sds/components";
import { useComponentWidth } from "./components/common/hooks/useIsComponentPastBreakpoint";
import { DEFAULT_ONTOLOGY_HEIGHT } from "../common/OntologyDagView/common/constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import CellGuideInfoSideBar from "../CellGuideInfoSideBar";
import { CellType } from "../common/OntologyDagView/common/types";

const SDS_INPUT_DROPDOWN_PROPS: InputDropdownProps = {
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
  const router = useRouter();

  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);
  const [selectedGene, setSelectedGene] = useState<string | undefined>(
    undefined
  );

  // Navigation
  const sectionRef0 = React.useRef<HTMLDivElement>(null);
  const sectionRef1 = React.useRef<HTMLDivElement>(null);
  const sectionRef2 = React.useRef<HTMLDivElement>(null);
  const sectionRef3 = React.useRef<HTMLDivElement>(null);

  const selectGene = (gene: string) => {
    if (gene === selectedGene) {
      setSelectedGene(undefined);
    } else {
      setSelectedGene(gene);
      if (sectionRef1.current) {
        window.scrollTo({
          top:
            sectionRef1.current.getBoundingClientRect().top +
            window.scrollY -
            50,
          behavior: "smooth",
        });
      }
    }
  };
  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);

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
  }, [skinnyMode]);

  // Set the mobile tooltip view content
  const [tooltipContent, setTooltipContent] = useState<{
    title: string;
    element: JSX.Element;
  } | null>(null);

  // cell type id
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypeName = name || "";
  const titleizedCellTypeName = titleize(cellTypeName);

  const handleResize = useCallback(() => {
    setSkinnyMode(
      window.innerWidth <
        SKINNY_MODE_BREAKPOINT_WIDTH + 2 * LEFT_RIGHT_PADDING_PX_XXL
    );
  }, []);

  const throttledHandleResize = useMemo(() => {
    return throttle(handleResize, 100);
  }, [handleResize]);

  useEffect(() => {
    throttledHandleResize();
    window.addEventListener("resize", throttledHandleResize);

    return () => window.removeEventListener("resize", throttledHandleResize);
  }, [throttledHandleResize]);

  const [geneInfoGene, setGeneInfoGene] = useState<Gene["name"] | null>(null);
  const [cellInfoCellType, setCellInfoCellType] = useState<CellType | null>(
    null
  );

  const { organismsList, organsMap } =
    useOrganAndOrganismFilterListForCellType(cellTypeId);

  const sdsOrganismsList = useMemo(
    () =>
      organismsList.map((organism) => ({
        name: organism,
      })),
    [organismsList]
  );

  const sdsOrgansList = useMemo(
    () =>
      Array.from(organsMap.keys()).map((organ) => ({
        name: organ === ALL_TISSUES ? TISSUE_AGNOSTIC : organ,
      })),
    [organsMap]
  );

  const [selectedOrgan, setSelectedOrgan] = useState<DefaultDropdownMenuOption>(
    sdsOrgansList.find(
      (organ) => organ.name === TISSUE_AGNOSTIC
    ) as DefaultDropdownMenuOption
  );

  const [selectedOrganId, setSelectedOrganId] = useState(NO_ORGAN_ID);

  const handleChangeOrgan = (option: DefaultDropdownMenuOption | null) => {
    if (!option || option.name === selectedOrgan.name) return;
    setSelectedOrgan(option);
    const optionName =
      option.name === TISSUE_AGNOSTIC ? ALL_TISSUES : option.name;
    setSelectedOrganId(organsMap.get(optionName) ?? "");
    // Continue tracking the analytics event as All Tissues
    track(EVENTS.CG_SELECT_TISSUE, { tissue: optionName });
  };

  const [selectedOrganism, setSelectedOrganism] =
    useState<DefaultDropdownMenuOption>(sdsOrganismsList[0]);

  const handleChangeOrganism = (option: DefaultDropdownMenuOption | null) => {
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

  useEffect(() => {
    setSelectedGene(undefined);
  }, [selectedOrgan, selectedOrganism, setSelectedGene]);

  const title = `${titleizedCellTypeName} Cell Types - CZ CELLxGENE CellGuide`;
  const seoDescription = `Find comprehensive information about "${cellTypeName}" cell types (synonyms: ${
    synonyms?.join(", ") || "N/A"
  }). ${rawSeoDescription}`;

  const dropdownComponents = (
    <CellGuideCardHeaderInnerWrapper>
      <Dropdown
        InputDropdownProps={SDS_INPUT_DROPDOWN_PROPS}
        search
        label={selectedOrganism?.name}
        onChange={handleChangeOrganism}
        options={sdsOrganismsList}
        value={selectedOrganism}
        data-testid={CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN}
      />
      <Dropdown
        InputDropdownProps={SDS_INPUT_DROPDOWN_PROPS}
        search
        label={selectedOrgan?.name}
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
          content="width=device-width, initial-scale=1, minimum-scale=1, maximum-scale=1"
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
            <ButtonIcon
              onClick={() => {
                setTooltipContent(null);
              }}
              sdsIcon="xMark"
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
          <CellGuideCardHeader width={width}>
            <CellGuideCardHeaderInnerWrapper>
              <CellGuideCardName data-testid={CELL_GUIDE_CARD_HEADER_NAME}>
                {titleizedCellTypeName}
              </CellGuideCardName>
              <a
                href={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
                target="_blank"
                rel="noreferrer noopener"
              >
                <StyledTag
                  data-testid={CELL_GUIDE_CARD_HEADER_TAG}
                  label={cellTypeId}
                  sdsType="secondary"
                  sdsStyle="square"
                  color="gray"
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
          <Wrapper skinnyMode={skinnyMode} ref={containerRef}>
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
                  key={`${cellTypeId}-${selectedOrganId}`}
                  cellTypeId={cellTypeId}
                  tissueName={selectedOrgan.name}
                  tissueId={selectedOrganId}
                  inputWidth={width}
                  setCellInfoCellType={setCellInfoCellType}
                  inputHeight={DEFAULT_ONTOLOGY_HEIGHT}
                  selectedOrganism={selectedOrganism.name}
                  selectedGene={selectedGene}
                  selectGene={selectGene}
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
              organName={selectedOrgan.name}
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
            selectedOrganName={selectedOrgan.name}
            selectedOrganId={selectedOrganId}
            organismName={selectedOrganism.name}
            selectedGene={selectedGene}
            selectGene={selectGene}
            setTooltipContent={setTooltipContent}
            skinnyMode={skinnyMode}
          />

          {
            // Split right sidebar view if fmg AND gene info is populated
            geneInfoGene && (
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
