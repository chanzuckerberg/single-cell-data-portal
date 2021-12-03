import { FC } from "react";
import Footer from "../Footer";
import Header from "../Header";
import { Wrapper } from "./style";

const Layout: FC = ({ children }) => {
  return (
    <Wrapper>
      <Header />
      {children}
      <Footer />
    </Wrapper>
  );
};

export default Layout;
