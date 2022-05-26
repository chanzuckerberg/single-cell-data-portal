import { useRouter } from "next/router";
import { FC } from "react";
import Layout from "src/components/Layout";
import { DefaultMainWrapper } from "../style";

const DefaultLayout: FC = ({ children }): JSX.Element => {
  const { pathname } = useRouter();

  // REMOVES DEFAULT MAIN WRAPPER ON LANDING PAGE WHICH DISABLES STICKY NAV
  // CHANGE TO "/" FOR PROD
  if (pathname === "/landing-page") {
    return (
      <Layout>
        <main>{children}</main>
      </Layout>
    );
  } else {
    return (
      <Layout>
        <DefaultMainWrapper>
          <main>{children}</main>
        </DefaultMainWrapper>
      </Layout>
    );
  }
};

export default DefaultLayout;
