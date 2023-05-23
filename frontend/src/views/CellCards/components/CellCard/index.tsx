import React from "react";
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

export const CELL_CARD_HEADER_NAME = "cell-card-header-name";
export const CELL_CARD_HEADER_TAG = "cell-card-header-tag";

export default function CellCard(): JSX.Element {
  const router = useRouter();

  // cell type id
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById() ?? {};
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

  return (
    <CellCardsView>
      <Wrapper>
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
          <SearchBarWrapper>
            <CellCardSearchBar />
          </SearchBarWrapper>
        </CellCardHeader>
        <Description cellTypeId={cellTypeId} cellTypeName={cellTypeName} />
        <CanonicalMarkerGeneTable cellTypeId={cellTypeId} />
      </Wrapper>
    </CellCardsView>
  );
}
