import React, { useEffect, useState } from "react";
import { useRouter } from "next/router";
import {
  Wrapper,
  CellCardName,
  CellCardHeader,
  StyledTag,
  CellCardsView,
  CellCardHeaderInnerWrapper,
  SearchBarWrapper,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import Description from "./components/Description";
import CellCardSearchBar from "../CellCardSearchBar";
import CanonicalMarkerGeneTable from "./components/CanonicalMarkerGeneTable";
import EnrichedGenesTable from "./components/EnrichedGenesTable";
import OntologyDagView from "./components/OntologyDagView";
import FullScreenProvider from "./components/FullScreenProvider";
import SourceDataTable from "./components/SourceDataTable";
import CellCardSidebar from "./components/CellCardSidebar";

export const CELL_CARD_HEADER_NAME = "cell-card-header-name";
export const CELL_CARD_HEADER_TAG = "cell-card-header-tag";

export default function CellCard(): JSX.Element {
  const router = useRouter();

  // Navigation
  const sectionRef0 = React.useRef(null);
  const sectionRef1 = React.useRef(null);
  const sectionRef2 = React.useRef(null);
  const sectionRef3 = React.useRef(null);
  const sectionRef4 = React.useRef(null);

  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);
  // cell type id
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById() ?? {};
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

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
    <CellCardsView>
      {/* Flex item left */}
      <Wrapper>
        <div ref={sectionRef0} id="section-0">
          <CellCardHeader>
            <CellCardHeaderInnerWrapper>
              <CellCardName data-testid={CELL_CARD_HEADER_NAME}>
                {cellTypeName.charAt(0).toUpperCase() + cellTypeName.slice(1)}
              </CellCardName>
              <a
                href={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
                target="_blank"
              >
                <StyledTag
                  data-testid={CELL_CARD_HEADER_TAG}
                  label={cellTypeId}
                  sdsType="secondary"
                  sdsStyle="square"
                  color="gray"
                  hover
                />
              </a>
            </CellCardHeaderInnerWrapper>
            {skinnyMode && (
              <SearchBarWrapper>
                <CellCardSearchBar />
              </SearchBarWrapper>
            )}
          </CellCardHeader>
          <Description cellTypeId={cellTypeId} cellTypeName={cellTypeName} />
        </div>

        <div ref={sectionRef1} id="section-1">
          <FullScreenProvider>
            <OntologyDagView cellTypeId={cellTypeId} skinnyMode={skinnyMode} />
          </FullScreenProvider>
        </div>
        <div ref={sectionRef2} id="section-2">
          <CanonicalMarkerGeneTable cellTypeId={cellTypeId} />
        </div>
        <div ref={sectionRef3} id="section-3">
          <EnrichedGenesTable cellTypeId={cellTypeId} />
        </div>
        <div ref={sectionRef4} id="section-4">
          <SourceDataTable cellTypeId={cellTypeId} />
        </div>
      </Wrapper>

      {/* Flex item right */}
      {!skinnyMode && (
        <CellCardSidebar
          items={[
            { elementRef: sectionRef0, title: "Intro" },
            { elementRef: sectionRef1, title: "Cell Ontology" },
            { elementRef: sectionRef2, title: "Marker Genes" },
            { elementRef: sectionRef3, title: "Highly Expressed Genes" },
            { elementRef: sectionRef4, title: "Source Data" },
          ]}
        />
      )}
    </CellCardsView>
  );
}
