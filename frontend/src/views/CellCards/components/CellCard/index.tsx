import React, { useEffect } from "react";
import { useRouter } from "next/router";
import {
  Wrapper,
  CellCardName,
  CellCardHeader,
  StyledTag,
  Divider,
  CellCardsView,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import Description from "./components/Description";
import CanonicalMarkerGeneTable from "./components/CanonicalMarkerGeneTable";
import EnrichedGenesTable from "./components/EnrichedGenesTable";
import SourceDataTable from "./components/SourceDataTable";
import CellCardSidebar, {
  INTRO_SECTION_ID,
} from "./components/CellCardSidebar";
import OntologyDagView from "./components/OntologyDagView";

// enum of available descriptions

export default function CellCard(): JSX.Element {
  const router = useRouter();

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
          </CellCardHeader>
          <Divider />
          <Description />
          <OntologyDagView
            cellTypeId={cellTypeId.replace(":", "_")}
            width={1040}
            height={500}
          />
          <CanonicalMarkerGeneTable cellTypeId={cellTypeId} />
          <EnrichedGenesTable cellTypeId={cellTypeId} />
          <SourceDataTable cellTypeId={cellTypeId} />
        </Wrapper>

        {/* Flex item right */}
        <CellCardSidebar />
      </CellCardsView>
    </>
  );
}
