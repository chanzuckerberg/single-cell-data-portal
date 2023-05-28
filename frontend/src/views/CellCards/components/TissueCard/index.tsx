import { useRouter } from "next/router";
import { useTissueCards } from "src/common/queries/cellCards";
import { TISSUE_CARD_MAX_WIDTH, TissueCardHeader, TissueCardHeaderInnerWrapper, TissueCardName, StyledTag, Wrapper, SearchBarWrapper } from "./style";
import CellCardSearchBar from "../CellCardSearchBar";
import OntologyDagView from "../common/OntologyDagView";
import FullScreenProvider from "../common/FullScreenProvider";

export const TISSUE_CARD_HEADER_NAME = "tissue-card-header-name";
export const TISSUE_CARD_HEADER_TAG = "tissue-card-header-tag";

export default function TissueCard(): JSX.Element {
  const router = useRouter();

  // cell type id
  const { tissueId: tissueIdRaw } = router.query;
  const tissueId = (tissueIdRaw as string)?.replace("_", ":") ?? "";
  const { data: tissueCardsData } = useTissueCards();
  const tissueName = tissueCardsData?.[tissueIdRaw as keyof typeof tissueCardsData]?.name ?? tissueId;
  
  return (
    <Wrapper>
      <TissueCardHeader>
        <TissueCardHeaderInnerWrapper>
          <TissueCardName data-testid={TISSUE_CARD_HEADER_NAME}>
            {tissueName.charAt(0).toUpperCase() + tissueName.slice(1)}
          </TissueCardName>
          <a
            href={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${tissueIdRaw}`}
            target="_blank"
          >
            <StyledTag
              data-testid={TISSUE_CARD_HEADER_TAG}
              label={tissueId}
              sdsType="secondary"
              sdsStyle="square"
              color="gray"
              hover
            />
          </a>
        </TissueCardHeaderInnerWrapper>
        <SearchBarWrapper>
          <CellCardSearchBar />
        </SearchBarWrapper>
      </TissueCardHeader>
      <FullScreenProvider>
        <OntologyDagView tissueId={tissueId} skinnyMode={false} initialWidth={TISSUE_CARD_MAX_WIDTH} initialHeight={560}/>
      </FullScreenProvider>
    </Wrapper>
  );
}
