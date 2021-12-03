import { FC } from "react";
import Layout from "src/components/Layout";
import { DefaultMainWrapper } from "../style";

const DefaultLayout: FC = ({ children }): JSX.Element => {
  return (
    <Layout>
      <DefaultMainWrapper>
        <main>{children}</main>
      </DefaultMainWrapper>
    </Layout>
  );
};

export default DefaultLayout;
