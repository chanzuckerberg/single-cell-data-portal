import { useRouter } from "next/router";
import { FC } from "react";
import Footer from "../Footer";
import Header from "../Header";
import LandingFooter from "../LandingFooter";
import LandingHeader from "../LandingHeader";
import { Wrapper } from "./style";

const Layout: FC = ({ children }) => {
  const { pathname } = useRouter();

  // CHANGE TO "/" FOR PROD
  if (pathname === "/landing-page") {
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
        <Footer />
      </Wrapper>
    );
  }
};

export default Layout;
