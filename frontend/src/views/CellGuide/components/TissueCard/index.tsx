import { useCallback, useEffect, useMemo, useState } from "react";
import { useRouter } from "next/router";
import {
  TissueCardHeader,
  TissueCardHeaderInnerWrapper,
  TissueCardName,
  StyledTag,
  Wrapper,
  SearchBarWrapper,
  DescriptionWrapper,
  LEFT_RIGHT_PADDING_PX_XXL,
  TissueCardView,
  TissueCardWrapper,
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
import { Global } from "@emotion/react";
import { StickySidebarStyle } from "../CellGuideCard/components/CellGuideCardSidebar/style";
import CellGuideMobileHeader from "../CellGuideMobileHeader";
import { SKINNY_MODE_BREAKPOINT_WIDTH } from "../CellGuideCard/constants";
import { throttle } from "lodash";
import CellGuideCardSidebar from "../CellGuideCard/components/CellGuideCardSidebar";
import React from "react";
import { StyledOntologyId } from "../CellGuideCard/components/Description/style";
import {
  TableTitleWrapper,
  TableTitle,
} from "../CellGuideCard/components/common/style";
import { useComponentWidth } from "../CellGuideCard/components/common/hooks/useIsComponentPastBreakpoint";

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

  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);
  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);

  const handleResize = useCallback(() => {
    setSkinnyMode(
      window.innerWidth <
        SKINNY_MODE_BREAKPOINT_WIDTH + 2 * LEFT_RIGHT_PADDING_PX_XXL
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

  // Navigation
  const sectionRef0 = React.useRef(null);
  const sectionRef1 = React.useRef(null);

  const { width, containerRef } = useComponentWidth();
  return (
    <>
      {/* Intro section */}
      <div ref={sectionRef0} id="section-0" data-testid="section-0" />
      {skinnyMode && (
        <CellGuideMobileHeader
          title={titleizedName}
          pageNav={
            <CellGuideCardSidebar
              skinnyMode={skinnyMode}
              items={[
                { elementRef: sectionRef0, title: "Intro" },
                { elementRef: sectionRef1, title: "Ontology" },
              ]}
            />
          }
          pageNavIsOpen={pageNavIsOpen}
          setPageNavIsOpen={setPageNavIsOpen}
        />
      )}
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
      <TissueCardWrapper skinnyMode={skinnyMode}>
        <TissueCardView skinnyMode={skinnyMode}>
          <Wrapper skinnyMode={skinnyMode} ref={containerRef}>
            {!skinnyMode && (
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
            )}

            {/* Don't show search on page for mobile view */}
            {!skinnyMode && (
              <SearchBarWrapper>
                <CellGuideCardSearchBar />
              </SearchBarWrapper>
            )}

            <DescriptionWrapper>
              <CellGuideCardDescription
                data-testid={TISSUE_CARD_UBERON_DESCRIPTION}
              >
                {description}
                <Source>
                  <SourceLink>
                    {"Source: "}
                    <Link
                      url={`https://www.ebi.ac.uk/ols/ontologies/uberon/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F${tissueIdRaw}`}
                      label={"UBERON Ontology"}
                    />
                  </SourceLink>
                </Source>
              </CellGuideCardDescription>
            </DescriptionWrapper>

            {skinnyMode && (
              <StyledOntologyId
                url={`https://www.ebi.ac.uk/ols/ontologies/uberon/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F${tissueIdRaw}`}
                ontologyId={tissueId}
              />
            )}

            {/* Ontology section */}
            <div ref={sectionRef1} id="section-1" data-testid="section-1">
              <TableTitleWrapper>
                <TableTitle>Ontology</TableTitle>
              </TableTitleWrapper>
              <FullScreenProvider>
                <OntologyDagView
                  key={tissueId}
                  tissueId={tissueId}
                  tissueName={tissueName}
                  inputWidth={width}
                  inputHeight={height}
                />
              </FullScreenProvider>
            </div>
          </Wrapper>
        </TissueCardView>
      </TissueCardWrapper>
      <CellGuideBottomBanner includeSurveyLink={!skinnyMode} />
    </>
  );
}
