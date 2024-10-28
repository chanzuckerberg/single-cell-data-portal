import { useRouter } from "next/router";
import { ReactNode } from "react";
import { ROUTES } from "src/common/constants/routes";
import Header from "../Header";
import LandingFooter from "../LandingFooter";
import LandingHeader from "../MobileFriendlyHeader";
import { Wrapper } from "./style";
import BottomBanner from "../BottomBanner";
import { BANNER_FEEDBACK_SURVEY_LINK } from "src/common/constants/airtableLinks";

interface Props {
  children: ReactNode;
}

const Layout = ({ children }: Props) => {
  const { pathname } = useRouter();

  if (pathname === ROUTES.HOMEPAGE || pathname === ROUTES.SITEMAP) {
    return (
      <>
        <LandingHeader
          data-testid="landing-footer"
          banner={true}
          bannerOptions={{
            content: (
              <p>
                New feature: Visualize and analyze spatial transcriptomics
                datasets in Explorer!&nbsp;
                <a href="https://cellxgene.cziscience.com/datasets?utm_campaign=spatial_jul_2024?utm_source=homepage_banner&utm_medium=product">
                  Browse datasets
                </a>
                &nbsp;using the Visium or Slide-seqV2 assay filter.
              </p>
            ),
          }}
        />
        {children}
        <LandingFooter />
        <BottomBanner surveyLink={BANNER_FEEDBACK_SURVEY_LINK} />
      </>
    );
  } else if (pathname === ROUTES.CELL_GUIDE) {
    return (
      <Wrapper data-testid="global-layout-wrapper">
        <LandingHeader title="CellGuide" />
        {children}
      </Wrapper>
    );
  } else if (pathname.startsWith(ROUTES.CELL_GUIDE)) {
    return (
      <>
        <LandingHeader
          title="CellGuide"
          homeUrl={ROUTES.HOMEPAGE}
          labelUrl={ROUTES.CELL_GUIDE}
        />
        {children}
      </>
    );
  } else {
    return (
      <Wrapper data-testid="global-layout-wrapper">
        <Header />
        {children}
        {/* <Footer /> */}
      </Wrapper>
    );
  }
};

export default Layout;
