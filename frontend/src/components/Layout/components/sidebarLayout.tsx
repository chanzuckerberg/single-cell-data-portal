import { FC } from "react";
import Layout from "src/components/Layout";
import { SidebarMainWrapper } from "../style";

const SidebarLayout: FC = ({ children }): JSX.Element => {
  return (
    <Layout>
      <SidebarMainWrapper>
        <main>{children}</main>
      </SidebarMainWrapper>
    </Layout>
  );
};

export default SidebarLayout;
