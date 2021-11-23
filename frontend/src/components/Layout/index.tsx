import { FC } from "react";
import Footer from "../Footer";
import Header from "../Header";
import { MainWrapper, Wrapper } from "./style";

const Layout: FC = ({ children }) => {
  return (
    <Wrapper>
      <Header />
      <MainWrapper>{children}</MainWrapper>
      <Footer />
    </Wrapper>
  );
};

export default Layout;
