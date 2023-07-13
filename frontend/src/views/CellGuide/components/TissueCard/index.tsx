import { useEffect, useState } from "react";
import { useRouter } from "next/router";
import {
  useTissuesById,
  useUberonDescription,
} from "src/common/queries/cellGuide";
import {
  TISSUE_CARD_MAX_WIDTH,
  TissueCardHeader,
  TissueCardHeaderInnerWrapper,
  TissueCardName,
  StyledTag,
  Wrapper,
  SearchBarWrapper,
  DescriptionWrapper,
  LEFT_RIGHT_PADDING_PX,
} from "./style";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import OntologyDagView from "../common/OntologyDagView";
import FullScreenProvider from "../common/FullScreenProvider";
import {
  CellGuideCardDescription,
  Source,
  SourceLink,
} from "../CellGuideCard/components/Description/style";
import Link from "../CellGuideCard/components/common/Link";

export const TISSUE_CARD_HEADER_NAME = "tissue-card-header-name";
export const TISSUE_CARD_HEADER_TAG = "tissue-card-header-tag";
export const TISSUE_CARD_UBERON_DESCRIPTION = "tissue-card-uberon-description";

export default function TissueCard(): JSX.Element {
  const router = useRouter();

  // cell type id
  const { tissueId: tissueIdRaw } = router.query;
  const tissueId = (tissueIdRaw as string)?.replace("_", ":") ?? "";
  const tissuesById = useTissuesById();
  const tissueName =
    tissuesById?.[tissueId as keyof typeof tissuesById]?.label ?? tissueId;

  // get current height of viewport
  const [height, setHeight] = useState(1000);
  useEffect(() => {
    const onResize = () => setHeight(window.innerHeight - 200);
    onResize();
    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, []);

  const [descriptionUberon, setDescriptionUberon] = useState<string>("");
  const { data: rawDescriptionUberon } = useUberonDescription(tissueId);

  useEffect(() => {
    if (rawDescriptionUberon) setDescriptionUberon(rawDescriptionUberon);
    else setDescriptionUberon("");
  }, [rawDescriptionUberon]);

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
      </TissueCardHeader>
      <SearchBarWrapper>
        <CellGuideCardSearchBar />
      </SearchBarWrapper>
      <DescriptionWrapper>
        <CellGuideCardDescription data-testid={TISSUE_CARD_UBERON_DESCRIPTION}>
          {descriptionUberon}
          <Source>
            <SourceLink>
              {"Source: "}
              <Link
                url={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${tissueIdRaw}`}
                label={"UBERON Ontology"}
              />
            </SourceLink>
          </Source>
        </CellGuideCardDescription>
      </DescriptionWrapper>
      <FullScreenProvider>
        <OntologyDagView
          tissueId={tissueId}
          tissueName={tissueName}
          skinnyMode={false}
          initialWidth={TISSUE_CARD_MAX_WIDTH - LEFT_RIGHT_PADDING_PX * 2}
          initialHeight={height}
        />
      </FullScreenProvider>
    </Wrapper>
  );
}
