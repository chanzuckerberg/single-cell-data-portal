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
            <CellCardName>
              {cellTypeName.charAt(0).toUpperCase() + cellTypeName.slice(1)}
            </CellCardName>
            <a
              href={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
              target="_blank"
            >
              <StyledTag
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
        <Description cellTypeName={cellTypeName} />
      </Wrapper>
    </CellCardsView>
  );
}
