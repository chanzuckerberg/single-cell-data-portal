import { useRouter } from "next/router";
import { ReactNode, useReducer } from "react";
import { ROUTES } from "src/common/constants/routes";
import Header from "../Header";
import LandingFooter from "../LandingFooter";
import LandingHeader from "../LandingHeader";
import { Wrapper } from "./style";
import CellGuideMobileHeader from "src/views/CellGuide/components/CellGuideMobileHeader";
import {
  DispatchContext,
  INITIAL_STATE,
  StateContext,
  reducer,
} from "src/views/CellGuide/common/store";

interface Props {
  children: ReactNode;
}

const Layout = ({ children }: Props) => {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);
  const { pathname } = useRouter();

  if (pathname === ROUTES.HOMEPAGE || pathname === ROUTES.SITEMAP) {
    return (
      <>
        <LandingHeader data-testid="landing-footer" />
        {children}
        <LandingFooter />
      </>
    );
  } else if (pathname === ROUTES.CELL_GUIDE) {
    return (
      <Wrapper data-testid="global-layout-wrapper">
        <LandingHeader title="CellGuide" />
        {children}
        {/* <Footer /> */}
      </Wrapper>
    );
  } else if (pathname === ROUTES.CELL_GUIDE + "/[cellTypeId]") {
    return (
      <DispatchContext.Provider value={dispatch}>
        <StateContext.Provider value={state}>
          <Wrapper data-testid="global-layout-wrapper">
            <LandingHeader title="CellGuide" homeUrl={ROUTES.CELL_GUIDE} />
            <CellGuideMobileHeader />
            {children}
          </Wrapper>
        </StateContext.Provider>
      </DispatchContext.Provider>
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
