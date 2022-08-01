import { useRouter } from "next/router";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import Layout from "src/components/Layout";
import { DefaultMainWrapper } from "../style";

const DefaultLayout: FC = ({ children }): JSX.Element => {
  const { pathname } = useRouter();

  // REMOVES DEFAULT MAIN WRAPPER ON LANDING PAGE WHICH DISABLES STICKY NAV
  if (pathname === ROUTES.HOMEPAGE) {
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
