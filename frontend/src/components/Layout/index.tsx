import { useRouter } from "next/router";
import { ReactNode } from "react";
import { ROUTES } from "src/common/constants/routes";
import Header from "../Header";
import LandingFooter from "../LandingFooter";
import LandingHeader from "../MobileFriendlyHeader";
import { Wrapper } from "./style";

interface Props {
  children: ReactNode;
}

const Layout = ({ children }: Props) => {
  const { pathname } = useRouter();

  const showBanner = () => {
    const now = new Date();
    const startDate = new Date("2025-09-25");
    const endDate = new Date("2026-09-25");
    return now >= startDate && now <= endDate;
  };

  if (pathname === ROUTES.HOMEPAGE || pathname === ROUTES.SITEMAP) {
    return (
      <>
        <LandingHeader
          data-testid="landing-footer"
          banner={showBanner()}
          bannerOptions={{
            content: (
              <p>
                New! Explore our{" "}
                <a
                  href="https://virtualcellmodels.cziscience.com/benchmarks"
                  target="_blank"
                >
                  benchmarking
                </a>{" "}
                suite to evaluate single-cell foundation model performance on
                tasks like cell clustering and cell type classification.
              </p>
            ),
          }}
        />
        {children}
        <LandingFooter />
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
