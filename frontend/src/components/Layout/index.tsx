import React, { FC } from "react";
import Footer from "../Footer";
import Header from "../Header";
import { MainWrapper, Wrapper } from "./style";

const Layout: FC = ({ children }) => {
  return (
    <Wrapper>
      <Header />
      <MainWrapper>
        <main>{children}</main>
      </MainWrapper>
      <Footer />
    </Wrapper>
  );
};

export default Layout;
