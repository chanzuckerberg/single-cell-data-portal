import Head from "next/head";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import { StyledHeader, Wrapper } from "./style";
import CellGuideBottomBanner from "../CellGuideBottomBanner";
import { LANDING_PAGE_HEADER } from "src/views/CellGuide/components/LandingPage/constants";
import CellGuideMobileHeader from "../CellGuideMobileHeader";
import { useCallback, useEffect, useMemo, useState } from "react";
import { throttle } from "lodash";
import { MAX_WIDTH_BREAKPOINT_PX } from "src/components/LandingHeader/style";

const TITLE = "CellGuide Cell Types and Cell Tissues - CZ CELLxGENE";

const DESCRIPTION =
  "Explore single-cell transcriptomics data in CellGuide, a comprehensive resource that empowers researchers with deep insights into the intricacies of cell types";

export default function LandingPage(): JSX.Element {
  // This is for mobile only
  const [openMobileSearch, setOpenMobileSearch] = useState(false);

  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);

  const handleResize = useCallback(() => {
    setSkinnyMode(
      window.innerWidth < MAX_WIDTH_BREAKPOINT_PX // This is the value the global header nav condenses for mobile
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

  return (
    <div
      onBlur={() => {
        if (openMobileSearch) {
          setOpenMobileSearch(false);
        }
      }}
    >
      <Head>
        <title>{TITLE}</title>
        <meta property="title" key="title" content={TITLE} />
        <meta property="og:title" key="og:title" content={TITLE} />
        <meta property="twitter:title" key="twitter:title" content={TITLE} />

        <meta name="description" key="description" content={DESCRIPTION} />
        <meta
          property="og:description"
          key="og:description"
          content={DESCRIPTION}
        />
        <meta
          property="twitter:description"
          key="twitter:description"
          content={DESCRIPTION}
        />

        {/* This prevents auto zooming on the input box on mobile */}
        <meta
          name="viewport"
          content="width=device-width, initial-scale=1, minimum-scale=1, maximum-scale=1"
        />
      </Head>

      {/* When clicking landing page input box, open the search at the top */}
      {skinnyMode && openMobileSearch && (
        <CellGuideMobileHeader
          title=""
          pageNav={null}
          openSearch={true}
          top={0}
        />
      )}

      <Wrapper searchBarOpen={openMobileSearch}>
        <StyledHeader data-testid={LANDING_PAGE_HEADER}>
          CellGuide is a <br />
          comprehensive resource for <br />
          knowledge about cell types
        </StyledHeader>
        <div
          onClick={() => {
            if (skinnyMode) {
              setOpenMobileSearch(true);
            }
          }}
        >
          {/* Search will open at top of page for mobile */}
          {!openMobileSearch && (
            <CellGuideCardSearchBar
              skinnyModeBreakpointWidth={MAX_WIDTH_BREAKPOINT_PX}
            />
          )}
        </div>
      </Wrapper>
      <CellGuideBottomBanner includeSurveyLink={!skinnyMode} />
    </div>
  );
}
