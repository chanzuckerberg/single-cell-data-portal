import "@blueprintjs/icons/lib/css/blueprint-icons.css";
import React, { FC } from "react";
import "src/global.scss";
import { ThemeProvider } from "theme-ui";
import Footer from "../Footer";
import Header from "../Header";
import "../layout.css";
import theme from "../theme";
import { MainWrapper, Wrapper } from "./style";

const Layout: FC = ({ children }) => {
  return (
    <ThemeProvider theme={theme}>
      <Wrapper>
        <Header />
        <MainWrapper>
          <main>{children}</main>
        </MainWrapper>
        <Footer />
      </Wrapper>
    </ThemeProvider>
  );
};

export default Layout;
