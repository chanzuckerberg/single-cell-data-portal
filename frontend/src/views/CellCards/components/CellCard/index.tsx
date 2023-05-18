import React, { useEffect, useState } from "react";
import { useRouter } from "next/router";
import {
  Wrapper,
  CellCardName,
  CellCardHeader,
  StyledTag,
  CellCardsView,
  CellCardHeaderInnerWrapper,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import Description from "./components/Description";
import CanonicalMarkerGeneTable from "./components/CanonicalMarkerGeneTable";
import EnrichedGenesTable from "./components/EnrichedGenesTable";
import SourceDataTable from "./components/SourceDataTable";
import CellCardSidebar, {
  INTRO_SECTION_ID,
} from "./components/CellCardSidebar";
import FullScreenProvider from "./components/FullScreenProvider";
import OntologyDagView from "./components/OntologyDagView";
import CellCardSearchBar from "../CellCardSearchBar";
import { SearchBarWrapper } from "./components/CellCardSidebar/style";

export default function CellCard(): JSX.Element {
  const router = useRouter();

  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);
  // cell type id
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById() ?? {};
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

  useEffect(() => {
    const element = document.getElementById("global-layout-wrapper");
    if (element) {
      element.scrollTo(0, 0);
    }
  }, [cellTypeId]);

  useEffect(() => {
    const width = window.innerWidth;
    if (width < 1040) {
      setSkinnyMode(true);
    } else if (width >= 1040) {
      setSkinnyMode(false);
    }

    const handleResize = () => {
      const width = window.innerWidth;
      if (width < 1040) {
        setSkinnyMode(true);
      } else if (width >= 1040) {
        setSkinnyMode(false);
      }
    };
    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  }, []);

  return (
    <>
      <style>
        {`
          /* Hack because main has a global overflow CSS prop which interferes with sticky sidebar and scroll listener */
          main {
            overflow: unset !important;
          }
          html {
            height: unset !important;
          }
        `}
      </style>
      <CellCardsView>
        {/* Flex item left */}
        <Wrapper>
          <CellCardHeader id={INTRO_SECTION_ID}>
            <CellCardHeaderInnerWrapper>
              <CellCardName>
                {cellTypeName.charAt(0).toUpperCase() + cellTypeName.slice(1)}
              </CellCardName>
              <StyledTag
                label={cellTypeId}
                sdsType="secondary"
                sdsStyle="square"
                color="gray"
                hover={false}
              />
            </CellCardHeaderInnerWrapper>
            {skinnyMode && (
              <SearchBarWrapper>
                <CellCardSearchBar />
              </SearchBarWrapper>
            )}
          </CellCardHeader>
          <Description cellTypeName={cellTypeName} />
          <FullScreenProvider>
            <OntologyDagView
              skinnyMode={skinnyMode}
              cellTypeId={cellTypeId.replace(":", "_")}
            />
          </FullScreenProvider>
          <CanonicalMarkerGeneTable cellTypeId={cellTypeId} />
          <EnrichedGenesTable cellTypeId={cellTypeId} />
          <SourceDataTable cellTypeId={cellTypeId} />
        </Wrapper>

        {/* Flex item right */}
        {!skinnyMode && <CellCardSidebar />}
      </CellCardsView>
    </>
  );
}
