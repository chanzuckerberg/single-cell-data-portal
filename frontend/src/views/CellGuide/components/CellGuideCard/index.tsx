import React, { useCallback, useEffect, useMemo, useState } from "react";
import { useRouter } from "next/router";
import {
  Wrapper,
  CellGuideCardName,
  CellGuideCardHeader,
  StyledTag,
  CellGuideView,
  CellGuideCardHeaderInnerWrapper,
  SearchBarWrapper,
  LEFT_RIGHT_PADDING_PX,
  SearchBarPositioner,
  StyledRightSideBar,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellGuide";
import Description from "./components/Description";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import MarkerGeneTables from "./components/MarkerGeneTables";
import OntologyDagView from "../common/OntologyDagView";
import FullScreenProvider from "../common/FullScreenProvider";
import SourceDataTable from "./components/SourceDataTable";
import CellGuideCardSidebar from "./components/CellGuideCardSidebar";
import { Gene } from "src/views/WheresMyGene/common/types";
import { throttle } from "lodash";
import GeneInfoSideBar from "src/components/GeneInfoSideBar";

export const CELL_GUIDE_CARD_HEADER_NAME = "cell-guide-card-header-name";
export const CELL_GUIDE_CARD_HEADER_TAG = "cell-guide-card-header-tag";
const RIGHT_SIDEBAR_WIDTH_PX = 400;

// This is the desired width of the CellGuideCard components right after the sidebar is hidden.
const BREAKPOINT_WIDTH = 960;

export default function CellGuideCard(): JSX.Element {
  const router = useRouter();

  // Navigation
  const sectionRef0 = React.useRef(null);
  const sectionRef1 = React.useRef(null);
  const sectionRef2 = React.useRef(null);
  const sectionRef3 = React.useRef(null);

  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);
  // cell type id
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById() ?? {};
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

  const handleResize = useCallback(() => {
    setSkinnyMode(
      window.innerWidth < BREAKPOINT_WIDTH + 2 * LEFT_RIGHT_PADDING_PX
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

  function handleCloseGeneInfoSideBar() {
    setGeneInfoGene(null);
  }

  return (
    <>
      <CellGuideView skinnyMode={skinnyMode}>
        {/* Flex item left */}
        <Wrapper>
          {/* (thuang): Somehow we need a parent to prevent error:
          NotFoundError: Failed to execute 'insertBefore' on 'Node'
         */}
          <div>
            {skinnyMode && (
              <SearchBarPositioner>
                <SearchBarWrapper>
                  <CellGuideCardSearchBar />
                </SearchBarWrapper>
              </SearchBarPositioner>
            )}
          </div>
          {/* Intro section */}
          <div ref={sectionRef0} id="section-0" data-testid="section-0" />
          <CellGuideCardHeader>
            <CellGuideCardHeaderInnerWrapper>
              <CellGuideCardName data-testid={CELL_GUIDE_CARD_HEADER_NAME}>
                {cellTypeName}
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
          </CellGuideCardHeader>
          <Description cellTypeId={cellTypeId} cellTypeName={cellTypeName} />

          {/* Cell Ontology section */}
          <div ref={sectionRef1} id="section-1" data-testid="section-1" />
          <FullScreenProvider>
            <OntologyDagView cellTypeId={cellTypeId} skinnyMode={skinnyMode} />
          </FullScreenProvider>

          {/* Marker Genes section */}
          <div ref={sectionRef2} id="section-2" data-testid="section-2" />
          <MarkerGeneTables
            cellTypeId={cellTypeId}
            setGeneInfoGene={setGeneInfoGene}
            cellTypeName={cellTypeName}
          />

          {/* Source Data section */}
          <div ref={sectionRef3} id="section-3" data-testid="section-3" />
          <SourceDataTable cellTypeId={cellTypeId} />
        </Wrapper>
        {!skinnyMode && (
          <CellGuideCardSidebar
            items={[
              { elementRef: sectionRef0, title: "Intro" },
              { elementRef: sectionRef1, title: "Cell Ontology" },
              { elementRef: sectionRef2, title: "Marker Genes" },
              { elementRef: sectionRef3, title: "Source Data" },
            ]}
          />
        )}
      </CellGuideView>
      <StyledRightSideBar width={RIGHT_SIDEBAR_WIDTH_PX}>
        {geneInfoGene && (
          <GeneInfoSideBar
            geneInfoGene={geneInfoGene}
            handleClose={handleCloseGeneInfoSideBar}
            title={geneInfoGene}
          />
        )}
      </StyledRightSideBar>
    </>
  );
}
