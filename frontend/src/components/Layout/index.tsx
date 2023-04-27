import { useRouter } from "next/router";
import { ReactNode } from "react";
import { ROUTES } from "src/common/constants/routes";
import Header from "../Header";
import LandingFooter from "../LandingFooter";
import LandingHeader from "../LandingHeader";
import { Wrapper } from "./style";

interface Props {
  children: ReactNode;
}

const Layout = ({ children }: Props) => {
  const { pathname } = useRouter();

  if (pathname === ROUTES.HOMEPAGE || pathname === ROUTES.SITEMAP) {
    return (
      <>
        <LandingHeader data-testid="landing-footer" />
        {children}
        <LandingFooter />
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
