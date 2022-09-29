import { useRouter } from "next/router";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import Header from "../Header";
import LandingFooter from "../LandingFooter";
import LandingHeader from "../LandingHeader";
import { Wrapper } from "./style";

const Layout: FC = ({ children }) => {
  const { pathname } = useRouter();

  if (pathname === ROUTES.HOMEPAGE) {
    return (
      <>
        <LandingHeader />
        {children}
        <LandingFooter />
      </>
    );
  } else {
    return (
      <Wrapper>
        <Header />
        {children}
        {/* <Footer /> */}
      </Wrapper>
    );
  }
};

export default Layout;
