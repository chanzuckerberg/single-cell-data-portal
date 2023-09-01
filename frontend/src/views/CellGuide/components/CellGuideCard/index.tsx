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
  LEFT_RIGHT_PADDING_PX,
  StyledRightSideBar,
  StyledSynonyms,
  FlexContainer,
  StyledOntologyId,
  MobileTooltipTitle,
  MobileTooltipWrapper,
  MobileTooltipHeader,
  NavBarDropdownWrapper,
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
import { useCellTypeMetadata } from "src/common/queries/cellGuide";
import {
  SKINNY_MODE_BREAKPOINT_WIDTH,
  CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_HEADER_NAME,
  CELL_GUIDE_CARD_HEADER_TAG,
  CELL_GUIDE_CARD_SYNONYMS,
} from "src/views/CellGuide/components/CellGuideCard/constants";
import { useOrganAndOrganismFilterListForCelltype } from "./components/MarkerGeneTables/hooks/common";
import {
  ALL_TISSUES,
  NO_ORGAN_ID,
} from "./components/MarkerGeneTables/constants";
import {
  DefaultDropdownMenuOption,
  Dropdown,
  InputDropdownProps,
  ButtonIcon,
} from "@czi-sds/components";

const RIGHT_SIDEBAR_WIDTH_PX = 400;

const SDS_INPUT_DROPDOWN_PROPS: InputDropdownProps = {
  sdsStyle: "square",
} as InputDropdownProps;

interface Props {
  name: string;
  seoDescription: string;
}

export default function CellGuideCard({
  // From getServerSideProps
  name,
  // From getServerSideProps
  seoDescription: rawSeoDescription,
}: Props): JSX.Element {
  const router = useRouter();

  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);

  // Navigation
  const sectionRef0 = React.useRef(null);
  const sectionRef1 = React.useRef(null);
  const sectionRef2 = React.useRef(null);
  const sectionRef3 = React.useRef(null);

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

  const { data: cellTypesById } = useCellTypeMetadata();

  const cellType = cellTypesById && cellTypesById[cellTypeId];

  const { synonyms } = cellType || {};

  const handleResize = useCallback(() => {
    setSkinnyMode(
      window.innerWidth <
        SKINNY_MODE_BREAKPOINT_WIDTH + 2 * LEFT_RIGHT_PADDING_PX
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

  const { organismsList, organsMap } =
    useOrganAndOrganismFilterListForCelltype(cellTypeId);

  const sdsOrganismsList = organismsList.map((organism) => ({
    name: organism,
  }));

  const sdsOrgansList = Array.from(organsMap.keys()).map((organ) => ({
    name: organ,
  }));

  const [selectedOrgan, setSelectedOrgan] = useState<DefaultDropdownMenuOption>(
    {
      name: ALL_TISSUES,
    } as DefaultDropdownMenuOption
  );
  const [selectedOrganId, setSelectedOrganId] = useState(NO_ORGAN_ID);

  const handleChangeOrgan = (option: DefaultDropdownMenuOption | null) => {
    if (!option) return;
    setSelectedOrgan(option);
    setSelectedOrganId(organsMap.get(option.name) ?? "");
  };

  const [selectedOrganism, setSelectedOrganism] =
    useState<DefaultDropdownMenuOption>(sdsOrganismsList[0]);

  const handleChangeOrganism = (option: DefaultDropdownMenuOption | null) => {
    if (!option) return;
    setSelectedOrganism(option);
  };

  function handleCloseGeneInfoSideBar() {
    setGeneInfoGene(null);
  }

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

      <CellGuideView skinnyMode={skinnyMode}>
        {/* Flex item left */}
        <Wrapper>
          {/* (thuang): Somehow we need a parent to prevent error:
          NotFoundError: Failed to execute 'insertBefore' on 'Node'
         */}
          {/* Intro section */}
          <div ref={sectionRef0} id="section-0" data-testid="section-0" />
          {!skinnyMode && (
            <CellGuideCardHeader>
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
          <Description
            cellTypeId={cellTypeId}
            cellTypeName={cellTypeName}
            skinnyMode={skinnyMode}
            setTooltipContent={setTooltipContent}
          />

          <FlexContainer>
            {skinnyMode && (
              <StyledOntologyId
                url={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
                ontologyId={cellTypeId}
              />
            )}
            <StyledSynonyms
              synonyms={synonyms}
              data-testid={CELL_GUIDE_CARD_SYNONYMS}
            />
          </FlexContainer>

          {/* Cell Ontology section */}
          <div ref={sectionRef1} id="section-1" data-testid="section-1" />
          {/* (thuang): Somehow we need a parent <div /> to prevent error:
          NotFoundError: Failed to execute 'insertBefore' on 'Node'
         */}
          <div>
            <FullScreenProvider>
              <OntologyDagView
                key={`${cellTypeId}-${selectedOrganId}`}
                cellTypeId={cellTypeId}
                tissueName={selectedOrgan.name}
                tissueId={selectedOrganId}
                skinnyMode={skinnyMode}
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
          />

          {/* Source Data section */}
          <div ref={sectionRef3} id="section-3" data-testid="section-3" />
          <SourceDataTable
            cellTypeId={cellTypeId}
            organName={selectedOrgan.name}
            organId={selectedOrganId}
            organismName={selectedOrganism.name}
            skinnyMode={skinnyMode}
            setTooltipContent={setTooltipContent}
          />
        </Wrapper>

        {/* Side bar */}
        {!skinnyMode && cellGuideSideBar}
      </CellGuideView>
      <StyledRightSideBar
        width={RIGHT_SIDEBAR_WIDTH_PX}
        skinnyMode={skinnyMode}
      >
        {geneInfoGene && (
          <GeneInfoSideBar
            geneInfoGene={geneInfoGene}
            handleClose={handleCloseGeneInfoSideBar}
            title={geneInfoGene}
          />
        )}
      </StyledRightSideBar>
      <CellGuideBottomBanner includeSurveyLink={!skinnyMode} />
    </>
  );
}
