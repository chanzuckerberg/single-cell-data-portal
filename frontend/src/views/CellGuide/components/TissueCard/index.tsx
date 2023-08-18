import { useEffect, useState } from "react";
import { useRouter } from "next/router";
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
import Head from "next/head";
import { titleize } from "src/common/utils/string";
import CellGuideBottomBanner from "../CellGuideBottomBanner";
import {
  TISSUE_CARD_HEADER_NAME,
  TISSUE_CARD_HEADER_TAG,
  TISSUE_CARD_UBERON_DESCRIPTION,
} from "src/views/CellGuide/components/TissueCard/constants";

interface Props {
  // From getServerSideProps
  description: string;
  // From getServerSideProps
  name: string;
}

export default function TissueCard({ description, name }: Props): JSX.Element {
  const router = useRouter();

  // cell type id
  const { tissueId: tissueIdRaw } = router.query;
  const tissueId = (tissueIdRaw as string)?.replace("_", ":") ?? "";
  const tissueName = name || tissueId;
  const titleizedName = titleize(tissueName);

  // get current height of viewport
  const [height, setHeight] = useState(1000);
  useEffect(() => {
    const onResize = () => setHeight(window.innerHeight - 200);
    onResize();
    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, []);

  const title = `${titleizedName} Tissue - CZ CELLxGENE CellGuide`;
  const seoDescription = `Find comprehensive information about ${tissueName} tissue: ${description}`;

  return (
    <>
      <Wrapper>
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
        </Head>
        <TissueCardHeader>
          <TissueCardHeaderInnerWrapper>
            <TissueCardName data-testid={TISSUE_CARD_HEADER_NAME}>
              {titleizedName}
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
          <CellGuideCardDescription
            data-testid={TISSUE_CARD_UBERON_DESCRIPTION}
          >
            {description}
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
            key={tissueId}
            tissueId={tissueId}
            tissueName={tissueName}
            skinnyMode={false}
            initialWidth={TISSUE_CARD_MAX_WIDTH - LEFT_RIGHT_PADDING_PX * 2}
            initialHeight={height}
          />
        </FullScreenProvider>
      </Wrapper>
      <CellGuideBottomBanner />
    </>
  );
}
